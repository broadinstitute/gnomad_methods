"""
Functions to group variants for constraint.

Covers consequence groups, VEP grouping fields, and genetic ancestry downsamplings.
"""

import functools
import logging
import operator
from typing import Dict, List, Optional, Tuple, Union

import hail as hl

from gnomad.assessment.summary_stats import generate_filter_combinations
from gnomad.utils.constraint.constants import CLASSIC_LOF_ANNOTATIONS
from gnomad.utils.vep import explode_by_vep_annotation, process_consequences

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_utils")
logger.setLevel(logging.INFO)


def build_constraint_consequence_groups(
    csq_expr: hl.expr.StringExpression,
    lof_modifier_expr: hl.expr.StringExpression,
    classic_lof_annotations: Tuple[str, ...] = CLASSIC_LOF_ANNOTATIONS,
    additional_groupings: Optional[
        Dict[str, Dict[str, hl.expr.BooleanExpression]]
    ] = None,
    additional_grouping_combinations: Optional[List[List[str]]] = None,
) -> Tuple[List[hl.expr.BooleanExpression], List[Dict[str, str]]]:
    """
    Build constraint consequence groups.

    Builds constraint groups based on the consequence expression and LoF
    modifier expression. By default, the following groups are built:

        - csq_set: synonymous_variant, missense_variant
        - lof: classic, hc_lc, hc

    The resulting meta and corresponding constraint group filters are:

        - {"csq_set": "syn"}: synonymous_variant
        - {"csq_set": "mis"}: missense_variant
        - {"lof": "classic"}: classic LoF annotations (stop_gained,
          splice_donor_variant, splice_acceptor_variant)
        - {"lof": "hc_lc"}: LOFTEE HC or LC modifier
        - {"lof": "hc"}: LOFTEE HC modifier only

    Additional groupings can be added via ``additional_groupings``, and
    grouping combinations via ``additional_grouping_combinations``.

    .. note::

        A combination of more than one grouping builds the **union** of the
        combined groups, not their intersection. For example, with an
        ``additional_groupings`` entry ``{"am": {"am_per_99": ...}}``, the
        combination ``["lof", "am"]`` produces groups such as
        ``{"lof": "hc", "am": "am_per_99"}`` that contain variants that are HC
        LoF **or** in the top 1% of AlphaMissense scores. This differs from
        frequency metadata (e.g. ``freq_meta``), where multiple keys mean all of
        them apply.

    :param csq_expr: VEP most severe consequence expression (e.g.,
        ``ht.most_severe_consequence``).
    :param lof_modifier_expr: LOFTEE modifier expression (e.g., ``ht.modifier``).
    :param classic_lof_annotations: Classic LoF annotations. Default is
        ``("stop_gained", "splice_donor_variant", "splice_acceptor_variant")``.
    :param additional_groupings: Additional groupings to add to the constraint
        groups. Default is None.
    :param additional_grouping_combinations: Additional grouping combinations to
        add to the constraint groups. Groupings combined in one entry are joined as a
        union (see the note above). Default is None.
    :return: Tuple of (constraint group filter expressions, meta dicts).
    """
    lof_classic_expr = hl.literal(set(classic_lof_annotations)).contains(csq_expr)
    lof_hc_expr = lof_modifier_expr == "HC"
    lof_hc_lc_expr = lof_hc_expr | (lof_modifier_expr == "LC")
    annotation_dict = {
        "csq_set": {
            "syn": csq_expr == "synonymous_variant",
            "mis": csq_expr == "missense_variant",
        },
        "lof": {
            "classic": lof_classic_expr,
            "hc_lc": lof_hc_lc_expr,
            "hc": lof_hc_expr,
        },
    }

    annotation_dict.update(additional_groupings or {})

    grouping_combinations = [["csq_set"], ["lof"]]
    grouping_combinations.extend(additional_grouping_combinations or [])

    meta = generate_filter_combinations(
        grouping_combinations,
        {k: list(v.keys()) for k, v in annotation_dict.items()},
    )
    # Multi-grouping combinations are unions, e.g. LoF or high-scoring missense.
    constraint_group_filters = [
        functools.reduce(operator.ior, [annotation_dict[k][v] for k, v in m.items()])
        for m in meta
    ]

    return constraint_group_filters, meta


def explode_downsamplings_oe(
    ht: hl.Table,
    downsampling_meta: Dict[str, List[str]],
    metrics: List[str] = ["syn", "lof", "mis"],
) -> hl.Table:
    """
    Explode observed and expected downsampling counts for each genetic ancestry group and metric.

    The input `ht` must contain struct of downsampling information for genetic ancestry
    groups under each metric name. For example: 'lof': struct {gen_anc_exp: struct
    {global: array<float64>}.

    :param ht: Input Table.
    :param metrics: List of metrics to explode. Default is '['syn', 'lof', 'mis']'.
    :param downsampling_meta: Dictionary containing downsampling metadata. Keys are the
        genetic ancestry group names and values are the list of downsamplings for that
        genetic ancestry group. Example: {'global': ['5000', '10000'], 'afr': ['5000',
        '10000']}.
    :return: Table with downsampling counts exploded so that observed and expected
        metric counts for each pair of genetic ancestry groups and downsamplings is a
        separate row.
    """
    ht = ht.select(
        _data=[
            hl.struct(
                gen_anc=gen_anc,
                downsampling=downsampling,
                **{
                    f"{metric}.{oe}": ht[metric][f"gen_anc_{oe}"][gen_anc][i]
                    for oe in ["obs", "exp"]
                    for metric in metrics
                },
            )
            for gen_anc, downsamplings in downsampling_meta.items()
            for i, downsampling in enumerate(downsamplings)
        ]
    )
    ht = ht.explode("_data")
    ht = ht.transmute(**ht._data)
    return ht


def get_all_gen_anc_lengths(
    ht: hl.Table,
    gen_ancs: Tuple[str, ...],
    obs_expr: hl.expr.StructExpression,
) -> List[Tuple[int, str]]:
    """
    Get the minimum length of observed variant counts array for each genetic ancestry group downsampling.

    The observed variant counts for each genetic ancestry group in `gen_ancs` are specified by
    annotations on the `obs_expr` expression.

    The function also performs a check that arrays of variant counts within genetic ancestry group
    downsamplings all have the same lengths.

    :param ht: Input Table containing `obs_expr`.
    :param gen_ancs: Genetic ancestry groups used to categorize observed variant counts in downsamplings.
    :param obs_expr: Expression for the genetic ancestry group observed variant counts. Should be a
        struct containing an array for each gen_anc in `gen_ancs`.
    :return: List of ``(length, gen_anc)`` tuples, one per group in ``gen_ancs``,
        where ``length`` is the minimum array length for that group.
    """
    # TODO: This function will be converted into doing just the length check if there
    #  is no usage of gen_anc_lengths in the constraint pipeline.
    # Get minimum length of downsamplings for each genetic ancestry group.
    gen_anc_downsampling_lengths = ht.aggregate(
        [hl.agg.min(hl.len(obs_expr[gen_anc])) for gen_anc in gen_ancs]
    )

    # Zip genetic ancestry group name with their downsampling length.
    gen_anc_lengths = list(zip(gen_anc_downsampling_lengths, gen_ancs))
    logger.info("Found: %s", "".join(map(str, gen_anc_lengths)))

    assert ht.all(
        hl.all(
            lambda f: f,
            [
                hl.len(obs_expr[gen_anc]) == length
                for length, gen_anc in gen_anc_lengths
            ],
        )
    ), (
        "The arrays of variant counts within genetic ancestry group downsamplings have different"
        " lengths!"
    )

    return gen_anc_lengths


def get_constraint_grouping_expr(
    vep_annotation_expr: hl.expr.StructExpression,
    coverage_expr: Optional[hl.expr.Int32Expression] = None,
    include_transcript_group: bool = True,
    include_canonical_group: bool = True,
    include_mane_select_group: bool = False,
) -> Dict[
    str,
    Union[hl.expr.StringExpression, hl.expr.Int32Expression, hl.expr.BooleanExpression],
]:
    """
    Collect annotations used for constraint groupings.

    Function collects the following annotations:

        - annotation - most_severe_consequence from ``vep_annotation_expr``
        - modifier - first non-missing of lof or polyphen_prediction
          from ``vep_annotation_expr``, or the literal "None"
        - gene - gene_symbol from ``vep_annotation_expr``
        - gene_id - gene_id from ``vep_annotation_expr``
        - coverage - exome coverage if ``coverage_expr`` is specified
        - transcript - transcript_id from ``vep_annotation_expr`` (added when
          ``include_transcript_group`` is True)
        - canonical - from ``vep_annotation_expr`` (added when
          ``include_canonical_group`` is True)
        - mane_select - from ``vep_annotation_expr`` (added when
          ``include_mane_select_group`` is True)

    .. note::

        This function expects that the following fields are present in
        ``vep_annotation_expr``:

            - lof
            - most_severe_consequence
            - gene_symbol
            - gene_id
            - polyphen_prediction (optional; missing used if absent)
            - transcript_id (if ``include_transcript_group`` is True)
            - canonical (if ``include_canonical_group`` is True)
            - mane_select (if ``include_mane_select_group`` is True)

    :param vep_annotation_expr: StructExpression of VEP annotation.
    :param coverage_expr: Int32Expression of exome coverage. Default is None.
    :param include_transcript_group: Whether to include the transcript annotation in the
        groupings. Default is True.
    :param include_canonical_group: Whether to include canonical annotation in the
        groupings. Default is True.
    :param include_mane_select_group: Whether to include mane_select annotation in the
        groupings. Default is False.
    :return: Dict mapping annotation names to Hail expressions.
    """
    lof_expr = vep_annotation_expr.lof
    if "polyphen_prediction" in vep_annotation_expr:
        polyphen_prediction_expr = vep_annotation_expr.polyphen_prediction
    else:
        polyphen_prediction_expr = hl.missing(hl.tstr)

    # Create constraint annotations to be used for groupings.
    groupings = {
        "annotation": vep_annotation_expr.most_severe_consequence,
        "modifier": hl.coalesce(lof_expr, polyphen_prediction_expr, "None"),
        "gene": vep_annotation_expr.gene_symbol,
        "gene_id": vep_annotation_expr.gene_id,
    }
    if coverage_expr is not None:
        groupings["coverage"] = coverage_expr

    # Add 'transcript' and 'canonical' annotation if requested.
    if include_transcript_group:
        groupings["transcript"] = vep_annotation_expr.transcript_id
    if include_canonical_group:
        groupings["canonical"] = hl.or_else(vep_annotation_expr.canonical == 1, False)
    if include_mane_select_group:
        groupings["mane_select"] = hl.is_defined(vep_annotation_expr.mane_select)

    return groupings


def annotate_exploded_vep_for_constraint_groupings(
    ht: hl.Table,
    coverage_expr: Optional[hl.expr.Int32Expression] = None,
    vep_annotation: str = "transcript_consequences",
    include_canonical_group: bool = True,
    include_mane_select_group: bool = False,
) -> Tuple[hl.Table, Tuple[str, ...]]:
    """
    Explode a VEP annotation and add constraint grouping fields.

    Explodes the specified VEP annotation (``vep_annotation``) and adds the
    following annotations via ``get_constraint_grouping_expr``:

        - annotation - most_severe_consequence from ``vep_annotation``
        - modifier - first non-missing of lof or polyphen_prediction
          from ``vep_annotation``, or the literal "None"
        - gene - gene_symbol from ``vep_annotation``
        - gene_id - gene_id from ``vep_annotation``
        - coverage - exome coverage if ``coverage_expr`` is specified
        - transcript - transcript_id from ``vep_annotation`` (added when
          ``vep_annotation`` is "transcript_consequences")
        - canonical - from ``vep_annotation`` (added when
          ``include_canonical_group`` is True)
        - mane_select - from ``vep_annotation`` (added when
          ``include_mane_select_group`` is True)

    .. note::

        This function expects that a ``vep`` annotation is present in ``ht``.

    .. note::

        Only the grouping fields are kept from the VEP annotation. The returned
        Table no longer has the ``vep`` annotation or the exploded
        ``vep_annotation`` struct, which keeps per-variant tables built from it
        small. All other fields of ``ht`` are kept.

    :param ht: Input Table.
    :param coverage_expr: Expression that defines the coverage metric. Default
        is None.
    :param vep_annotation: Name of annotation in vep (one of
        ``"transcript_consequences"`` and ``"worst_csq_by_gene"``) that will be
        used for obtaining constraint annotations. Default is
        ``"transcript_consequences"``.
    :param include_canonical_group: Whether to include canonical annotation
        in the groupings. Default is True. Ignored unless ``vep_annotation`` is
        ``"transcript_consequences"``.
    :param include_mane_select_group: Whether to include mane_select
        annotation in the groupings. Default is False. Ignored unless
        ``vep_annotation`` is ``"transcript_consequences"``.
    :return: Tuple of (Table with one row per element of ``vep_annotation`` and the
        grouping fields added in place of the VEP annotation, names of the added
        grouping fields).
    """
    # Annotate ht with coverage_expr set as a temporary annotation '_coverage_metric'
    # before modifying ht.
    if coverage_expr is not None:
        ht = ht.annotate(_coverage_metric=coverage_expr)

    if vep_annotation == "transcript_consequences":
        if not include_canonical_group and not include_mane_select_group:
            raise ValueError(
                "If 'vep_annotation' is 'transcript_consequences', one of either"
                " 'include_canonical_group' or 'include_mane_select_group' must be set!"
            )
        include_transcript_group = True
    else:
        logger.warning(
            "Setting both 'include_canonical_group' and 'include_mane_select_group' to"
            " False (options cannot be used unless 'vep_annotation' is"
            " 'transcript_consequences')."
        )
        include_transcript_group = False
        include_canonical_group = False
        include_mane_select_group = False

    # Annotate 'worst_csq_by_gene' to `ht` if it's specified for `vep_annotation`.
    if vep_annotation == "worst_csq_by_gene":
        ht = process_consequences(ht)

    # Explode the specified VEP annotation.
    ht = explode_by_vep_annotation(ht, vep_annotation)

    # Collect the annotations used for groupings.
    groupings = get_constraint_grouping_expr(
        ht[vep_annotation],
        coverage_expr=None if coverage_expr is None else ht._coverage_metric,
        include_transcript_group=include_transcript_group,
        include_canonical_group=include_canonical_group,
        include_mane_select_group=include_mane_select_group,
    )

    return ht.transmute(**groupings), tuple(groupings.keys())
