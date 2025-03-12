"""Script containing generic constraint functions that may be used in the constraint pipeline."""

import copy
import logging
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import hail as hl
from hail.utils.misc import divide_null, new_temp_file

from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.vep import (
    add_most_severe_consequence_to_consequence,
    explode_by_vep_annotation,
    process_consequences,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_utils")
logger.setLevel(logging.INFO)

COVERAGE_CUTOFF = 30
"""
Minimum median exome coverage differentiating high coverage sites from low coverage sites.

Low coverage sites require an extra calibration when computing the proportion of expected variation.
"""


def get_mu_annotation_expr(
    ht: hl.Table,
    mutation_ht: hl.Table,
    mu_expr: Union[str, hl.expr.Float64Expression] = "mu_snp",
) -> hl.expr.Float64Expression:
    """
    Get mutation rate annotation expression from the mutation rate Table.

    .. note::

        Function expects that `ht` includes `mutation_ht`'s key fields. Note that these
        annotations don't need to be the keys of `ht`. The standard keys used are:
        'context', 'ref', 'alt', and 'methylation_level'.

    This function converts the mutation rate Table into a dictionary and uses the
    mutation rate HT key fields in the input Table to get the mutation rate annotation
    expression. The dictionary is used instead of joining the mutation rate Table to
    the input Table to avoid unnecessary shuffling given the small size of the
    mutation rate Table.

    :param ht: Input Table that will be annotated with the returned mutation rate
        expression.
    :param mutation_ht: Mutation rate Table.
    :param mu_expr: Mutation rate expression or annotation name in `mutation_ht`.
        Default is 'mu_snp'.
    :return: Mutation rate annotation expression.
    """
    if isinstance(mu_expr, str):
        mu_expr = mutation_ht[mu_expr]

    mutation_dict = mutation_ht.aggregate(
        hl.agg.group_by(mutation_ht.key, hl.agg.collect(mu_expr)), _localize=False
    ).map_values(lambda x: x[0])

    return mutation_dict.get(hl.struct(**{k: ht[k] for k in mutation_ht.key}))


def annotate_with_mu(
    ht: hl.Table,
    mutation_ht: hl.Table,
    mu_annotation: str = "mu_snp",
) -> hl.Table:
    """
    Annotate SNP mutation rate for the input Table.

    .. note::

        Function expects that `ht` includes `mutation_ht`'s key fields. Note that these
        annotations don't need to be the keys of `ht`. The standard keys used are:
        'context', 'ref', 'alt', and 'methylation_level'.

    :param ht: Input Table to annotate.
    :param mutation_ht: Mutation rate Table.
    :param mu_annotation: The name of mutation rate annotation in `mutation_ht`.
        Default is 'mu_snp'.
    :return: Table with mutational rate annotation added.
    """
    mu = get_mu_annotation_expr(ht, mutation_ht, mutation_ht[mu_annotation])
    return ht.annotate(
        **{mu_annotation: hl.case().when(hl.is_defined(mu), mu).or_error("Missing mu")}
    )


def _get_expression_for_annotation(
    t: Optional[Union[hl.Table, hl.MatrixTable]] = None,
    annotation_name: Optional[str] = None,
    expr: Optional[hl.expr.Expression] = None,
    expr_param_name: Optional[str] = None,
) -> hl.expr.Expression:
    """
    Get an annotation from a Table or MatrixTable.

    Either `t` and `annotation_name` or `expr` must be provided. If `expr` is not None,
    it will be returned. Otherwise, the annotation with the name `annotation_name` will
    be returned from the input Table or MatrixTable.

    `expr_param_name` is an optional parameter that is only used in a logger message or
    error message if `expr` is None.

    :param t: Optional Input Table or MatrixTable.
    :param annotation_name: Optional Name of the annotation to get.
    :param expr: Optional expression to return. Default is None.
    :param expr_param_name: Optional name of the parameter that `expr` is associated
        with. Default is None.
    :return: Expression of the annotation.
    """
    if (t is None or annotation_name) and expr is None:
        raise ValueError("Either `t` and `annotation_name` or `expr` must be provided.")

    expr_param_name = expr_param_name or "expr"
    if expr is None and hasattr(t, annotation_name):
        logger.warning(
            f"{expr_param_name} was not provided, using '{annotation_name}'."
        )
        expr = t[annotation_name]
    elif expr is None:
        raise ValueError(
            f"{expr_param_name} was not provided and '{annotation_name}' annotation is "
            f"not present in the input Table."
        )

    return expr


def single_variant_count_expr(
    freq_expr: Optional[hl.expr.StructExpression] = None,
    ht: Optional[hl.Table] = None,
    singleton: bool = False,
    max_af: Optional[float] = None,
    count_missing: bool = False,
) -> Union[hl.Table, Any]:
    """
    Get the expression for a single variant count based on specified criteria.

    One of `ht` or `freq_expr` must be provided. If `freq_expr` is not provided, the
    function will use `ht.freq` if it exists.

    A single variant count can be 1 (variant meets criteria) or 0 (variant does not meet
    criteria). The criteria are specified by `singleton` and `max_af`.

    The variant count will be 1 if:
        - `singleton` is True and the variant is a singleton (`freq_expr.AC == 1`).
        - `max_af` is not None and the variant has an allele frequency less than or
          equal to `max_af` and has an allele count greater than 0 (`freq_expr.AF <=
          max_af` and `freq_expr.AC > 0`).
        - If `max_af` is None and `singleton` is False, and `freq_expr` is not None, the
          variant count will be 1 if the variant has an allele count greater than 0
          (`freq_expr.AC > 0`).
        - If `max_af` is None and `singleton` is False, and `freq_expr` is None, the
          variant count will be 1.

    :param freq_expr: Optional StructExpression with 'AC' and 'AF' annotations. Default
        is None.
    :param ht: Optional Input Table. Default is None.
    :param singleton: Whether to count singletons. Default is False.
    :param max_af: Maximum variant allele frequency to keep. By default, no cutoff is
        applied.
    :param count_missing: Whether to count missing frequency values. Default is False.
    :return: Expression for a single variant count.
    """
    if ht is None and freq_expr is None:
        raise ValueError("Either `ht` or `freq_expr` must be provided.")

    if max_af or singleton:
        freq_expr = _get_expression_for_annotation(ht, "freq", freq_expr, "freq_expr")
        if isinstance(freq_expr, hl.expr.ArrayExpression):
            freq_expr = freq_expr[0]

    if singleton:
        count_var = freq_expr.AC == 1, count_missing
    elif max_af:
        count_var = (freq_expr.AC > 0) & (freq_expr.AF <= max_af)
    elif freq_expr is not None:
        count_var = freq_expr.AC > 0
    else:
        count_var = True

    return hl.int(hl.or_else(count_var, count_missing))


def single_variant_observed_and_possible_expr(
    freq_expr: hl.ArrayExpression,
    max_af: Optional[float] = None,
    use_possible_adj: bool = True,
) -> hl.expr.StructExpression:
    """
    Get the observed and possible variant annotations.

    The function returns a struct indicating whether the variant should be included in
    the observed and possible variant counts. It uses `single_variant_count_expr`
    to determine the observed and possible variant counts for the frequency array.

    :param freq_expr: Frequency array expression.
    :param max_af: Optional maximum allele frequency to consider a variant as observed.
    :param use_possible_adj: Whether to use the full dataset high-quality genotypes
        (adj) frequency for the possible count (first element in the array). Otherwise,
        the possible count is calculated for all entries in the frequency array.
        Default is True.
    :return: Struct containing the observed and possible variant annotations.
    """
    pos_expr = hl.array([freq_expr[0]]) if use_possible_adj else freq_expr
    pos_expr = pos_expr.map(
        lambda x: single_variant_count_expr(
            freq_expr=x, max_af=max_af, count_missing=True
        )
    )
    pos_expr = pos_expr[0] if use_possible_adj else pos_expr
    return hl.struct(
        observed_variants=freq_expr.map(
            lambda x: single_variant_count_expr(freq_expr=x, max_af=max_af)
        ),
        possible_variants=pos_expr,
    )


def weighted_agg_sum_expr(
    expr: Union[hl.ArrayExpression, hl.NumericExpression],
    weight_expr: Union[hl.ArrayExpression, hl.NumericExpression],
) -> Union[hl.ArrayExpression, hl.NumericExpression]:
    """
    Get the expression for a weighted aggregate sum.

    The function will return the sum of the product of `expr` and `weight_expr`.

    The parameters `expr` and `weight_expr` can be either ArrayExpression or
    NumericExpression. If both are ArrayExpression, each element of the arrays will be
    multiplied together and summed. If one is ArrayNumericExpression and the other is
    NumericExpression, the NumericExpression will be multiplied by each element
    of the ArrayExpression and summed.

    :param expr: Expression to be weighted and summed.
    :param weight_expr: Expression to weight `expr` by.
    :return: Weighted aggregate sum expression.
    """
    expr_is_array = isinstance(expr, hl.expr.ArrayNumericExpression)
    weight_is_array = isinstance(weight_expr, hl.expr.ArrayNumericExpression)
    if expr_is_array and weight_is_array:
        return hl.agg.array_sum(hl.zip(expr, weight_expr).map(lambda x: x[0] * x[1]))
    elif not expr_is_array and not weight_is_array:
        return hl.agg.sum(expr * weight_expr)
    else:
        return hl.agg.array_sum(expr * weight_expr)


def get_counts_agg_expr(
    freq_expr: Optional[Union[hl.expr.ArrayExpression, hl.expr.StructExpression]],
    ht: Optional[hl.Table] = None,
    count_singletons: bool = False,
    max_af: Optional[float] = None,
    count_missing: bool = False,
) -> hl.expr.StructExpression:
    """
    Get aggregation expression for counting variants based on specified criteria.

    One of `ht` or `freq_expr` must be provided. If `freq_expr` is not provided, the
    function will use `ht.freq` if it exists.

    The aggregation expression will return a StructExpression with 'variant_count' and
    'singleton_count' annotations. 'variant_count' will be the count of variants that
    meet the specified criteria. 'singleton_count' will be the count of singleton
    variants if `count_singletons` is True.

    The criteria for counting variants are specified by `count_singletons` and `max_af`.

    The variant count will be 1 if:

        - `count_singletons` is True and the variant is a singleton
          (`freq_expr.AC == 1`).
        - `max_af` is not None and the variant has an allele frequency less than or
          equal to `max_af` and has an allele count greater than 0
          (`freq_expr.AF <= max_af` and `freq_expr.AC > 0`).
        - If `max_af` is None and `count_singletons` is False, and `freq_expr` is not
          None, the variant count will be 1 if the variant has an allele count greater
          than 0 (`freq_expr.AC > 0`).
        - If `max_af` is None and `count_singletons` is False, and `freq_expr` is None,
          the variant count will be 1.

    :param freq_expr: Optional ArrayExpression or StructExpression with 'AC' and 'AF'
        annotations. Default is None.
    :param ht: Optional Input Table. Default is None.
    :param count_singletons: Whether to count singletons. Default is False.
    :param max_af: Maximum variant allele frequency to keep. By default, no cutoff is
        applied.
    :param count_missing: Whether to count missing frequency values. Default is False.
    :return: StructExpression with 'variant_count' and 'singleton_count' annotations.
    """
    if ht is None and freq_expr is None:
        raise ValueError("Either `ht` or `freq_expr` must be provided.")

    if max_af or count_singletons:
        freq_expr = _get_expression_for_annotation(ht, "freq", freq_expr, "freq_expr")

    if max_af:
        logger.info(
            "The maximum variant allele frequency to be included in "
            "`variant_count` is %.3f.",
            max_af,
        )

    if count_singletons:
        logger.info(
            "Counting singleton variants and adding as 'singleton_count' annotation."
        )

    params = {"variant_count": {"singleton": False, "max_af": max_af}}
    if count_singletons:
        params["singleton_count"] = {"singleton": True}

    is_struct = isinstance(freq_expr, hl.expr.StructExpression)
    if is_struct:
        freq_expr = hl.array([freq_expr])

    count_expr = {
        k: freq_expr.map(
            lambda x: single_variant_count_expr(
                freq_expr=x, **p, count_missing=count_missing
            )
        )
        for k, p in params.items()
    }

    agg_expr = {k: hl.agg.array_sum(v) for k, v in count_expr.items()}
    agg_expr = {k: agg_expr[k][0] for k in agg_expr} if is_struct else agg_expr

    return hl.struct(**agg_expr)


def count_observed_and_possible_by_group(
    ht: hl.Table,
    possible_expr: hl.expr.Int32Expression,
    observed_expr: hl.expr.ArrayExpression,
    additional_grouping: Union[List[str]] = ("methylation_level",),
    partition_hint: int = 100,
    weight_exprs: Optional[
        Union[
            List[str],
            Dict[str, Union[hl.expr.ArrayExpression, hl.expr.NumericExpression]],
        ]
    ] = None,
    additional_agg_sum_exprs: Optional[
        Union[
            List[str],
            Dict[str, Union[hl.expr.ArrayExpression, hl.expr.NumericExpression]],
        ]
    ] = None,
) -> Union[hl.Table, Any]:
    if isinstance(weight_exprs, list):
        weight_exprs = {k: ht[k] for k in weight_exprs}
    if isinstance(additional_agg_sum_exprs, list):
        additional_agg_sum_exprs = {k: ht[k] for k in additional_agg_sum_exprs}

    weight_exprs = weight_exprs or {}
    additional_agg_sum_exprs = additional_agg_sum_exprs or {}

    # Build the grouping struct for the variant count aggregation.
    grouping = hl.struct(context=ht.context, ref=ht.ref, alt=ht.alt)
    grouping = grouping.annotate(
        **{g: ht[g] for g in additional_grouping if g not in grouping}
    )
    logger.info(
        "The following annotations will be used to group the input Table rows when"
        " counting variants: %s.",
        ", ".join(grouping.keys()),
    )

    agg_expr = {
        "observed_variants": hl.agg.array_sum(observed_expr),
        "possible_variants": hl.agg.sum(possible_expr),
    }

    # Update the possible variant count aggregation expression to include weighted sums
    # of possible variant counts.
    agg_expr.update(
        {
            k: weighted_agg_sum_expr(possible_expr, v)
            for k, v in weight_exprs.items()
        }
    )

    # Get sum aggregation expressions for requested fields.
    agg_expr.update(
        {
            k: (
                hl.agg.array_sum(v)
                if isinstance(v, hl.ArrayExpression)
                else hl.agg.sum(v)
            )
            for k, v in additional_agg_sum_exprs.items()
        }
    )

    # Apply each variant count aggregation in `agg_expr` to get counts for all
    # combinations of `grouping`.
    ht = ht.group_by(**grouping).partition_hint(partition_hint).aggregate(**agg_expr)

    return ht


# TODO: I think we should consider removing this or at least completely changing it
#  To remove pop and downsampling support, since that should just be handled as an
#  array, where the same thing is done.
def count_variants_by_group(
    ht: hl.Table,
    freq_expr: Optional[hl.expr.ArrayExpression] = None,
    freq_meta_expr: Optional[hl.expr.ArrayExpression] = None,
    count_singletons: bool = False,
    count_downsamplings: Tuple[str] = (),
    downsamplings: Optional[List[int]] = None,
    additional_grouping: Tuple[str] = (),
    partition_hint: int = 100,
    omit_methylation: bool = False,
    use_table_group_by: bool = False,
    singleton_expr: Optional[hl.expr.BooleanExpression] = None,
    max_af: Optional[float] = None,
) -> Union[hl.Table, Any]:
    """
    Count number of observed or possible variants by context, ref, alt, and optionally methylation_level.

    Performs variant count aggregations based on specified criteria
    (`count_singletons`, `count_downsamplings`, and `max_af`), and grouped by:
    'context', 'ref', 'alt', 'methylation_level' (optional), and all annotations
    provided in `additional_grouping`.

    If variant allele frequency information is required based on other parameter
    selections (described in detail below) and `freq_expr` is not supplied, `freq_expr`
    defaults to `ht.freq` if it exists.

    `freq_expr` should be an ArrayExpression of Structs with 'AC' and 'AF' annotations.
    This is the same format as the `freq` annotation that is created using
    `annotate_freq()`.

    Variant allele frequency information is needed when:
        - `max_af` is not None - `freq_expr[0].AF` is used to filter to only variants
          with a maximum allele frequency of `max_af` prior to counting variants. In
          the standard `freq` ArrayExpression annotated by `annotate_freq()`, this
          first element corresponds to the allele frequency information for high quality
          genotypes (adj).
        - `count_singletons` is True and `singleton_expr` is None - If singleton counts
          are requested and no expression is specified to determine whether a variant
          is a singleton, `singleton_expr` defaults to `freq_expr[0].AC == 1`. In the
          standard `freq` ArrayExpression annotated by `annotate_freq()`, this
          corresponds to allele count of only 1 in the callset after filtering to high
          quality genotypes.
        - `count_downsamplings` is not empty - When downsampling counts are requested,
          `freq_expr` needs to contain frequency information for downsamplings within
          each population requested. In addition to needing `freq_expr`, this also
          requires the use of `freq_meta_expr`. If `freq_meta_expr` is None,
          `freq_meta_expr` it defaults to `ht.freq_meta` if it exists. Similar to
          `freq_expr`, `freq_meta_expr` is expected to have the same format as
          the `freq_meta` global annotation that is created using `annotate_freq()`.
          `freq_meta_expr` is used to determine the index of allele frequency
          information within `freq_expr` for each population requested and it's
          downsamplings.

    This function will return a Table with annotations used for grouping ('context',
    'ref', 'alt', 'methylation_level' (optional), `additional_grouping`) and
    'variant_count' annotation.

    .. note::

        The following annotations should be present in `ht`:
            - ref - the reference allele
            - alt - the alternate base
            - context - trinucleotide genomic context
            - methylation_level - methylation level (optional if omit_methylation==True)
            - freq - allele frequency information (AC, AN, AF, homozygote count; not
              required if `freq_expr` is given)
            - freq_meta - an ordered list containing the frequency aggregation group
              for each element of the `freq` array row annotation (not required if
              `freq_meta_expr` is given)

    :param ht: Input Hail Table.
    :param freq_expr: ArrayExpression of Structs with 'AC' and 'AF' annotations. If
        `freq_expr` is None and any of `count_downsamplings`, `max_af`, and
        `count_singletons` is True, `freq_expr` would be `ht.freq`.
    :param freq_meta_expr: ArrayExpression of meta dictionaries corresponding to
        `freq_expr`. If `count_downsamplings` and `freq_meta_expr` is None,
        `freq_meta_expr` would be `ht.freq_meta`.
    :param count_singletons: Whether to count singletons (defined by `singleton_expr`).
        Default is False.
    :param count_downsamplings: Tuple of populations to use for downsampling counts.
        Default is ().
    :param downsamplings: Optional List of integers specifying what downsampling
        indices to obtain. Default is None, which will return all downsampling counts.
    :param additional_grouping: Additional features to group by. e.g. 'exome_coverage'.
        Default is ().
    :param partition_hint: Target number of partitions for aggregation. Default is 100.
    :param omit_methylation: Whether to omit 'methylation_level' from the grouping when
        counting variants. Default is False.
    :param use_table_group_by: Whether to group `ht` before aggregating the variant
        counts. If `use_table_group_by` is False, function will return a hl.
        StructExpression. Default is False.
    :param singleton_expr: Expression for defining a singleton. When `count_singletons`
        is True and `singleton_expr` is None, `singleton_expression` would be `freq_expr
        [0].AC == 1`. Default is None.
    :param max_af: Maximum variant allele frequency to keep. By default, no cutoff is
        applied.
    :return: Table including 'variant_count' annotation and if requested,
        `singleton_count` and downsampling counts.
    """
    if freq_expr is None and (
        count_downsamplings or max_af or (count_singletons and singleton_expr is None)
    ):
        logger.warning(
            "freq_expr was not provided, using 'freq' as the frequency annotation."
        )
        freq_expr = ht.freq
    if count_downsamplings and freq_meta_expr is None:
        logger.warning(
            "freq_meta_expr was not provided, using 'freq_meta' as the frequency"
            " metadata annotation."
        )
        freq_meta_expr = ht.freq_meta
    if count_singletons and singleton_expr is None:
        logger.warning(
            "count_singletons is True and singleton_expr was not provided, using"
            " freq_expr[0].AC == 1 as the singleton expression."
        )
        singleton_expr = freq_expr[0].AC == 1

    grouping = hl.struct(context=ht.context, ref=ht.ref, alt=ht.alt)
    if not omit_methylation:
        logger.info(
            "'methylation_level' annotation is included in the grouping when counting"
            " variants."
        )
        grouping = grouping.annotate(methylation_level=ht.methylation_level)
    for group in additional_grouping:
        grouping = grouping.annotate(**{group: ht[group]})
    logger.info(
        "The following annotations will be used to group the input Table rows when"
        " counting variants: %s.",
        ", ".join(grouping.keys()),
    )

    if max_af:
        logger.info(
            "The maximum variant allele frequency to be included in `variant_count` is"
            " %.3f.",
            max_af,
        )
        agg = {"variant_count": hl.agg.count_where(freq_expr[0].AF <= max_af)}
    else:
        agg = {"variant_count": hl.agg.count()}

    if count_singletons:
        logger.info(
            "Counting singleton variants and adding as 'singleton_count' annotation."
        )
        agg["singleton_count"] = hl.agg.count_where(singleton_expr)

    for pop in count_downsamplings:
        logger.info(
            "Counting variants in downsamplings for population '%s', and adding as"
            " 'downsampling_counts_%s' annotation.",
            pop,
            pop,
        )
        agg[f"downsampling_counts_{pop}"] = downsampling_counts_expr(
            freq_expr,
            freq_meta_expr,
            pop,
            max_af=max_af,
            downsamplings=downsamplings,
        )
        if count_singletons:
            logger.info(
                "Counting singleton variants in downsamplings for population '%s', and"
                " adding as 'singleton_downsampling_counts_%s' annotation.",
                pop,
                pop,
            )
            agg[f"singleton_downsampling_counts_{pop}"] = downsampling_counts_expr(
                freq_expr,
                freq_meta_expr,
                pop,
                max_af=max_af,
                downsamplings=downsamplings,
                singleton=True,
            )
    # Apply each variant count aggregation in `agg` to get counts for all
    # combinations of `grouping`.
    if use_table_group_by:
        return ht.group_by(**grouping).partition_hint(partition_hint).aggregate(**agg)
    else:
        return ht.aggregate(
            hl.struct(**{field: hl.agg.group_by(grouping, agg[field]) for field in agg})
        )


def get_downsampling_freq_indices(
    freq_meta_expr: hl.expr.ArrayExpression,
    pop: str = "global",
    variant_quality: str = "adj",
    genetic_ancestry_label: Optional[str] = None,
    subset: Optional[str] = None,
    downsamplings: Optional[List[int]] = None,
) -> hl.expr.ArrayExpression:
    """
    Get indices of dictionaries in meta dictionaries that only have the "downsampling" key with specified `genetic_ancestry_label` and "variant_quality" values.

    :param freq_meta_expr: ArrayExpression containing the set of groupings for each
        element of the `freq_expr` array (e.g., [{'group': 'adj'}, {'group': 'adj',
        'pop': 'nfe'}, {'downsampling': '5000', 'group': 'adj', 'pop': 'global'}]).
    :param pop: Population to use for filtering by the `genetic_ancestry_label` key in
        `freq_meta_expr`. Default is 'global'.
    :param variant_quality: Variant quality to use for filtering by the 'group' key in
        `freq_meta_expr`. Default is 'adj'.
    :param genetic_ancestry_label: Label defining the genetic ancestry groups. If None,
        "gen_anc" or "pop" is used (in that order of preference) if present. Default is
        None.
    :param subset: Subset to use for filtering by the 'subset' key in `freq_meta_expr`.
        Default is None, which will return all downsampling indices without a 'subset'
        key in `freq_meta_expr`.
    :param downsamplings: Optional List of integers specifying what downsampling
        indices to obtain. Default is None, which will return all downsampling indices.
    :return: ArrayExpression of indices of dictionaries in `freq_meta_expr` that only
        have the "downsampling" key with specified `genetic_ancestry_label` and
        "variant_quality" values.
    """
    if genetic_ancestry_label is None:
        gen_anc = ["gen_anc", "pop"]
    else:
        gen_anc = [genetic_ancestry_label]

    def _get_filter_expr(m: hl.expr.StructExpression) -> hl.expr.BooleanExpression:
        filter_expr = (
            (m.get("group") == variant_quality)
            & (hl.any([m.get(l, "") == pop for l in gen_anc]))
            & m.contains("downsampling")
        )
        if downsamplings is not None:
            filter_expr &= hl.literal(downsamplings).contains(
                hl.int(m.get("downsampling", "0"))
            )
        if subset is None:
            filter_expr &= ~m.contains("subset")
        else:
            filter_expr &= m.get("subset", "") == subset
        return filter_expr

    indices = hl.enumerate(freq_meta_expr).filter(lambda f: _get_filter_expr(f[1]))

    # Get an array of indices and meta dictionaries sorted by "downsampling" key.
    return hl.sorted(indices, key=lambda f: hl.int(f[1]["downsampling"]))


def downsampling_counts_expr(
    freq_expr: hl.expr.ArrayExpression,
    freq_meta_expr: hl.expr.ArrayExpression,
    pop: str = "global",
    variant_quality: str = "adj",
    singleton: bool = False,
    max_af: Optional[float] = None,
    genetic_ancestry_label: Optional[str] = None,
    subset: Optional[str] = None,
    downsamplings: Optional[List[int]] = None,
) -> hl.expr.ArrayExpression:
    """
    Return an aggregation expression to compute an array of counts of all downsamplings found in `freq_expr` where specified criteria is met.

    The frequency metadata (`freq_meta_expr`) should be in a similar format to the
    `freq_meta` annotation added by `annotate_freq()`. Each downsampling should have
    'group', `genetic_ancestry_label`, and 'downsampling' keys. Included downsamplings
    are those where 'group' == `variant_quality` and `genetic_ancestry_label` == `pop`.

    :param freq_expr: ArrayExpression of Structs with 'AC' and 'AF' annotations.
    :param freq_meta_expr: ArrayExpression containing the set of groupings for each
        element of the `freq_expr` array (e.g., [{'group': 'adj'}, {'group': 'adj',
        'pop': 'nfe'}, {'downsampling': '5000', 'group': 'adj', 'pop': 'global'}]).
    :param pop: Population to use for filtering by the `genetic_ancestry_label` key in
        `freq_meta_expr`. Default is 'global'.
    :param variant_quality: Variant quality to use for filtering by the 'group' key in
        `freq_meta_expr`. Default is 'adj'.
    :param singleton: Whether to filter to only singletons before counting (AC == 1).
        Default is False.
    :param max_af: Maximum variant allele frequency to keep. By default no allele
        frequency cutoff is applied.
    :param genetic_ancestry_label: Label defining the genetic ancestry groups. If None,
        "gen_anc" or "pop" is used (in that order of preference) if present. Default is
        None.
    :param subset: Subset to use for filtering by the 'subset' key in `freq_meta_expr`.
        Default is None, which will return all downsampling counts without a 'subset'
        key in `freq_meta_expr`. If specified, only downsamplings with the specified
        subset will be included.
    :param downsamplings: Optional List of integers specifying what downsampling
        indices to obtain. Default is None, which will return all downsampling counts.
    :return: Aggregation Expression for an array of the variant counts in downsamplings
        for specified population.
    """
    # Get an array of indices sorted by "downsampling" key.
    sorted_indices = get_downsampling_freq_indices(
        freq_meta_expr,
        pop,
        variant_quality,
        genetic_ancestry_label,
        subset,
        downsamplings,
    ).map(lambda x: x[0])

    def _get_criteria(i: hl.expr.Int32Expression) -> hl.expr.Int32Expression:
        """
        Return 1 when variant meets specified criteria (`singleton` or `max_af`), if requested, or with an AC > 0.

        :param i: The index of a downsampling.
        :return: Returns 1 if the variant in the downsampling with specified index met
            the criteria. Otherwise, returns 0.
        """
        if singleton:
            return hl.int(freq_expr[i].AC == 1)
        elif max_af:
            return hl.int((freq_expr[i].AC > 0) & (freq_expr[i].AF <= max_af))
        else:
            return hl.int(freq_expr[i].AC > 0)

    # Map `_get_criteria` function to each downsampling indexed by `sorted_indices` to
    # generate a list of 1's and 0's for each variant, where the length of the array is
    # the total number of downsamplings for the specified population and each element
    # in the array indicates if the variant in the downsampling indexed by
    # `sorted_indices` meets the specified criteria.
    # Return an array sum aggregation that aggregates arrays generated from mapping.
    return hl.agg.array_sum(hl.map(_get_criteria, sorted_indices))


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
                gen_anc=pop,
                downsampling=downsampling,
                **{
                    f"{metric}.{oe}": ht[metric][f"gen_anc_{oe}"][pop][i]
                    for oe in ["obs", "exp"]
                    for metric in metrics
                },
            )
            for pop, downsamplings in downsampling_meta.items()
            for i, downsampling in enumerate(downsamplings)
        ]
    )
    ht = ht.explode("_data")
    ht = ht.transmute(**ht._data)
    return ht


def annotate_mutation_type(
    t: Union[hl.MatrixTable, hl.Table],
    context_length: Optional[int] = None,
    num_scan_context_length: Optional[int] = 100,
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Annotate mutation types.

    The following annotations are added to the output Table:
        - cpg
        - transition
        - mutation_type - one of "CpG", "non-CpG transition", or "transversion"
        - mutation_type_model

    ..note:

        This function uses the term 'mutation_type' because 'variant_type' is already
        used in this repo to indicate a variant's multiallelic and SNP/indel status.

    :param t: Input Table or MatrixTable.
    :param context_length: Length of the 'context' annotation in 't'. If this is not
        specified, the value will be determined by examining the first
        `num_scan_context_length` values of the 'context' annotation. Default is None.
    :param num_scan_context_length: Number of values in the 'context' annotation to use
        for determining `context_length` if it is not specified. If set to None, all
        values in 'context' will be used. Default is 100.
    :return: Table with mutation type annotations added.
    """
    if context_length is None:
        # Determine the context length by collecting all the context lengths.
        if num_scan_context_length is None:
            context_lengths = t.aggregate(hl.agg.collect_as_set(hl.len(t.context)))
            msg = "all"
        else:
            context_lengths = hl.len(t.context).take(num_scan_context_length)
            msg = f"the first {num_scan_context_length}"
        context_lengths = list(filter(None, set(context_lengths)))
        if len(context_lengths) > 1:
            raise ValueError(
                f"More than one length was found among {msg} 'context' values. Length "
                "of 'context' should be consistent.",
            )
        else:
            context_length = context_lengths[0]
            logger.info(
                "Detected a length of %d for context length using %s 'context' values.",
                context_length,
                msg,
            )

    # Determine the middle index of the context annotation.
    if context_length == 3:
        mid_index = 1
    elif context_length == 7:
        mid_index = 3
    else:
        raise ValueError(
            "The length of context should be either 3 or 7, instead of"
            f" {context_length}."
        )

    transition_expr = hl.is_transition(t.ref, t.alt)
    cpg_expr = (
        (t.ref == "G") & (t.alt == "A") & (t.context[mid_index - 1 : mid_index] == "C")
    ) | (
        (t.ref == "C")
        & (t.alt == "T")
        & (t.context[mid_index + 1 : mid_index + 2] == "G")
    )
    if isinstance(t, hl.MatrixTable):
        t = t.annotate_rows(transition=transition_expr, cpg=cpg_expr)
    else:
        t = t.annotate(transition=transition_expr, cpg=cpg_expr)
    mutation_type_expr = (
        hl.switch(hl.len(t.context))
        .when(
            context_length,
            hl.case()
            .when(t.cpg, "CpG")
            .when(t.transition, "non-CpG transition")
            .default("transversion"),
        )
        .or_error("Found 'context' value with unexpected context length!")
    )
    mutation_type_model_expr = hl.if_else(t.cpg, t.context, "non-CpG")
    if isinstance(t, hl.MatrixTable):
        return t.annotate_rows(
            mutation_type=mutation_type_expr,
            mutation_type_model=mutation_type_model_expr,
        )
    else:
        return t.annotate(
            mutation_type=mutation_type_expr,
            mutation_type_model=mutation_type_model_expr,
        )


def trimer_from_heptamer(
    t: Union[hl.MatrixTable, hl.Table],
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Trim heptamer context to create trimer context.

    :param t: Input MatrixTable or Table with context annotation.
    :return: MatrixTable or Table with trimer context annotated.
    """
    trimer_expr = hl.if_else(hl.len(t.context) == 7, t.context[2:5], t.context)
    return (
        t.annotate_rows(context=trimer_expr)
        if isinstance(t, hl.MatrixTable)
        else t.annotate(context=trimer_expr)
    )


def collapse_strand(
    t: Union[hl.Table, hl.MatrixTable],
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Return the deduplicated context by collapsing DNA strands.

    Function returns the reverse complement for 'ref, 'alt', and 'context' if the
    reference allele is either 'G' or 'T'.

    The following annotations are added to the output Table:
        - was_flipped - whether the 'ref, 'alt', and 'context' were flipped (reverse
          complement taken)

    :param ht: Input Table.
    :return: Table with deduplicated context annotation (ref, alt, context,
        was_flipped).
    """
    ref_g_or_t_expr = (t.ref == "G") | (t.ref == "T")
    collapse_expr = {
        "ref": hl.if_else(ref_g_or_t_expr, hl.reverse_complement(t.ref), t.ref),
        "alt": hl.if_else(ref_g_or_t_expr, hl.reverse_complement(t.alt), t.alt),
        "context": hl.if_else(
            ref_g_or_t_expr,
            hl.reverse_complement(t.context),
            t.context,
        ),
        "was_flipped": ref_g_or_t_expr,
    }
    return (
        t.annotate(**collapse_expr)
        if isinstance(t, hl.Table)
        else t.annotate_rows(**collapse_expr)
    )


def transform_methylation_level(
    methylation_expr: Union[str, hl.expr.NumericExpression] = "methylation_level",
    methylation_cutoffs: Tuple[Union[int, float], Union[int, float]] = (0, 5),
    ht: Optional[hl.Table] = None,
) -> Union[hl.Table, hl.expr.NumericExpression]:
    """
    Transform methylation level into a 0 to 2 scale.

    The methylation level is transformed to a 0-2 scale based on the provided cutoffs.
    The methylation level is assigned a value of 2 if it is greater than the second
    cutoff, 1 if it is greater than the first cutoff (but less than or equal to the
    second), and 0 otherwise.

    :param methylation_expr: Methylation level expression or annotation name in `ht`. If
        `methylation_expr` is a string, `ht` must be provided.
    :param methylation_cutoffs: Tuple of two integers/floats representing the cutoffs
        for the methylation level transformation. Default is (0, 5).
    :param ht: Input Table. Default is None.
    :return: Table with methylation level annotation added or methylation level
        expression.
    """
    if isinstance(methylation_expr, str) and ht is None:
        raise ValueError("ht must be provided if methylation_expr is a string.")

    if isinstance(methylation_expr, str):
        methylation_expr = ht[methylation_expr]

    methylation_expr = (
        hl.case()
        .when(methylation_expr > methylation_cutoffs[1], 2)
        .when(methylation_expr > methylation_cutoffs[0], 1)
        .default(0)
    )

    if ht is not None:
        return ht.annotate(methylation_level=methylation_expr)

    return methylation_expr


def assemble_constraint_context_ht(
    ht: hl.Table,
    methylation_ht: hl.Table = None,
    gerp_ht: hl.Table = None,
    coverage_hts: Dict[str, hl.Table] = None,
    an_hts: Dict[str, hl.Table] = None,
    freq_hts: Dict[str, hl.Table] = None,
    filter_hts: Dict[str, hl.Table] = None,
    transformation_funcs: Optional[Dict[str, Callable]] = None,
) -> hl.Table:
    """
    Assemble context Table with necessary annotations for constraint calculations.

    .. note::

        Checks for 'was_split' annotation in Table. If not present, splits
        multiallelic sites.

    The following annotations are added to the output Table:

        - ref: Reference allele.
        - alt: Alternate allele.
        - context: Trimer context.
        - annotations added by `annotate_mutation_type()`, `collapse_strand()`, and
          `add_most_severe_csq_to_tc_within_vep_root()`.

    Depending on the annotations provided, the following annotations may also be added:

        - 'methylation_level': Methylation level annotation will be added if
          `methylation_ht` is provided.
        - 'gerp': GERP annotation will be added if `gerp_ht` is provided.
        - 'coverage': Coverage annotations will be added if `coverage_hts` is provided.
        - 'AN': Allele number annotations will be added if `an_hts` is provided.
        - 'freq': Frequency annotations will be added if `freq_hts` is provided.
        - 'filters': Filter annotations will be added if `filter_hts` is provided.

    The `transformation_funcs` parameter can be used to transform the HTs before adding
    annotations. The keys should be the annotation names ('coverage', 'gerp', etc.) and
    the values should be functions that take the annotation Table and the ht that is
    being annotated as input and return the transformed and keyed annotation. If not
    provided, the following default transformations are used:

        - 'methylation_level': Uses the 'methylation_level' annotation in the
          Methylation Table after transforming the methylation level to a 0-2 scale
          using `transform_grch37_methylation()` or `transform_grch38_methylation()`.
        - 'gerp': Uses the 'S' annotation in the GERP Table. If 'S' is missing, it
          defaults to 0.
        - 'coverage': If necessary, pulls out the first element of coverage statistics
          (which includes all samples). Relevant to v4, where coverage stats include
          additional elements to stratify by UK Biobank subset and platforms.
        - 'AN': Uses the 'AN' annotation in the allele number Table.
        - 'freq': Uses the 'freq' annotation in the frequency Table.
        - 'filters': Uses the 'filters' annotation in the filter Table.

    The following global annotations may also be added to the output Table:

        - 'an_globals': Global allele number annotations 'strata_sample_count' and
          'strata_meta' will be added if `an_hts` is provided.
        - 'freq_globals': Global frequency annotations 'freq_meta_sample_count' and
          'freq_meta' will be added if `freq_hts` is provided.

    :param ht: Input context Table with VEP annotation.
    :param methylation_ht: Optional Table with methylation annotation. Default is None.
    :param gerp_ht: Optional Table with GERP annotation. Default is None.
    :param coverage_hts: An optional Dictionary with key as one of 'exomes' or
        'genomes' and values as corresponding coverage Tables. Default is None.
    :param an_hts: A Dictionary with key as one of 'exomes' or 'genomes' and
        values as corresponding allele number Tables. Default is None.
    :param transformation_funcs: An optional Dictionary to transform the HTs before
        adding annotations. Default is None, which will use some default transformations
        as described above.
    :return: Table with sites split and necessary annotations.
    """
    ref = get_reference_genome(ht.locus)
    if ref == hl.get_reference("GRCh37"):
        from gnomad.resources.grch37 import transform_grch37_methylation as trans_methyl
    else:
        from gnomad.resources.grch38 import transform_grch38_methylation as trans_methyl

    # Check if context Table is split, and if not, split multiallelic sites.
    if "was_split" not in list(ht.row):
        ht = hl.split_multi_hts(ht)

    # Only need to keep the keys (locus, alleles) and the context and vep annotations.
    ht = ht.select("context", "vep")

    # Filter Table to only contigs 1-22, X, Y.
    ht = hl.filter_intervals(
        ht, [hl.parse_locus_interval(c, ref.name) for c in ref.contigs[:24]]
    )

    # Add annotations for 'ref' and 'alt'.
    ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1])

    # Trim heptamer context to create trimer context.
    ht = trimer_from_heptamer(ht)

    # Annotate mutation type (such as "CpG", "non-CpG transition", "transversion") and
    # collapse strands to deduplicate the context.
    ht = annotate_mutation_type(collapse_strand(ht))

    # Add most_severe_consequence annotation to 'transcript_consequences' within the
    # vep root annotation.
    csqs = ht.vep.transcript_consequences
    csqs = add_most_severe_consequence_to_consequence(csqs)
    vep_csq_fields = [
        "transcript_id",
        "gene_id",
        "gene_symbol",
        "biotype",
        "most_severe_consequence",
        "mane_select",
        "canonical",
        "lof",
        "lof_flags",
    ]
    vep_csq_fields = [x for x in vep_csq_fields if x in csqs.dtype.element_type]
    ht = ht.annotate(
        vep=ht.vep.select(
            "most_severe_consequence",
            transcript_consequences=csqs.map(lambda x: x.select(*vep_csq_fields)),
        )
    )

    # Add 'methylation_level' and 'gerp' annotations if specified.
    transformation_funcs = transformation_funcs or {}

    if "methylation_level" not in transformation_funcs:
        transformation_funcs["methylation_level"] = lambda x, t: hl.if_else(
            ht.cpg, trans_methyl(x)[t.locus].methylation_level, 0
        )

    if "gerp" not in transformation_funcs:
        transformation_funcs["gerp"] = lambda x, t: hl.if_else(
            hl.is_missing(x[t.locus].S), 0, x[t.locus].S
        )

    # If necessary, pull out first element of coverage statistics (which includes all
    # samples). Relevant to v4, where coverage stats include additional elements to
    # stratify by ukb subset and platforms.
    if "coverage" not in transformation_funcs:
        transformation_funcs["coverage"] = lambda x, t: (
            x[t.locus].coverage_stats[0] if "coverage_stats" in x.row else x[t.locus]
        )

    if "AN" not in transformation_funcs:
        transformation_funcs["AN"] = lambda x, t: x[t.locus].AN

    if "freq" not in transformation_funcs:
        transformation_funcs["freq"] = lambda x, t: x[t.key].freq

    if "filters" not in transformation_funcs:
        transformation_funcs["filters"] = lambda x, t: x[t.key].filters

    hts = {
        "methylation_level": methylation_ht,
        "gerp": gerp_ht,
        "coverage": coverage_hts,
        "AN": an_hts,
        "freq": freq_hts,
        "filters": filter_hts,
    }
    hts = {k: v for k, v in hts.items() if v is not None}
    exprs = {}
    for ann, ann_ht in hts.items():
        if isinstance(ann_ht, dict):
            ann_expr = hl.struct(
                **{k: transformation_funcs[ann](v, ht) for k, v in ann_ht.items()}
            )
        else:
            ann_expr = transformation_funcs[ann](ann_ht, ht)

        exprs[ann] = ann_expr

    ht = ht.annotate(**exprs)

    # Add global annotations for 'AN' and 'freq' if HTs are provided.
    global_anns = {
        "an_globals": (an_hts, ["strata_sample_count", "strata_meta"]),
        "freq_globals": (freq_hts, ["freq_meta_sample_count", "freq_meta"]),
    }
    global_anns = {k: v for k, v in global_anns.items() if v[0] is not None}
    ht = ht.annotate_globals(
        **{
            k: hl.struct(
                **{
                    data_type: ann_ht.select_globals(
                        *[x for x in g if x in ann_ht.globals]
                    ).index_globals()
                    for data_type, ann_ht in hts.items()
                }
            )
            for k, (hts, g) in global_anns.items()
        }
    )

    return ht


def calibration_model_group_expr(
    exomes_coverage_expr: hl.expr.Int32Expression,
    cpg_expr: hl.expr.BooleanExpression,
    low_cov_cutoff: Optional[int] = None,
    high_cov_cutoff: int = COVERAGE_CUTOFF,
    upper_cov_cutoff: Optional[int] = None,
    skip_coverage_model: bool = False,
    additional_grouping_exprs: Optional[Dict[str, hl.expr.StringExpression]] = None,
    cpg_in_high_only: bool = False,
) -> hl.expr.StructExpression:
    """
    Get the calibration model grouping annotation for a variant.

    The calibration model expression is a struct with the following fields:

        - genomic_region: The genomic region of the variant ("autosome_or_par",
          "chrx_nonpar", or "chry_nonpar").
        - high_or_low_coverage: Whether the variant belongs to the high or low coverage
          calibration model. The variant is assigned to the high coverage model if the
          exome coverage is greater than or equal to 'high_cov_cutoff' and less than or
          equal to 'upper_cov_cutoff' (if provided). The variant is assigned to the low
          coverage model if `skip_coverage_model` is False and the exome coverage is
          greater than 'low_cov_cutoff' (if provided) and less than 'high_cov_cutoff'.
        - cpg: Whether the variant is a CpG (`cpg_expr`).

    The global parameters for the calibration model are the values of the function
    parameters: `low_cov_cutoff`, `high_cov_cutoff`, `upper_cov_cutoff`, and
    `skip_coverage_model`.

    :param exomes_coverage_expr: Exome coverage expression.
    :param cpg_expr: CpG expression.
    :param low_cov_cutoff: Low coverage cutoff. Default is None.
    :param high_cov_cutoff: High coverage cutoff. Default is COVERAGE_CUTOFF.
    :param upper_cov_cutoff: Upper coverage cutoff. Default is None.
    :param skip_coverage_model: Whether to skip the coverage model. Default is False.
    :param additional_grouping_exprs: Optional Dictionary of additional expressions to
        group by. Default is None.
    :return: Tuple containing the calibration model expression and the globals.
    """
    high_cov_expr = exomes_coverage_expr >= high_cov_cutoff
    if upper_cov_cutoff is not None:
        high_cov_expr &= exomes_coverage_expr <= upper_cov_cutoff

    low_cov_expr = hl.bool(False) if skip_coverage_model else hl.bool(True)
    if low_cov_cutoff is not None:
        low_cov_expr &= exomes_coverage_expr > low_cov_cutoff

    # Define whether the variant should be included in the high or low coverage model.
    model_expr = (
        hl.case()
        .when(high_cov_expr, "high")
        .when(low_cov_expr,"low")
        .or_missing()
    )

    cpg_expr = hl.or_missing(model_expr == "high", cpg_expr) if cpg_in_high_only else cpg_expr
    return hl.or_missing(
        hl.is_defined(model_expr),
        hl.struct(
            high_or_low_coverage=model_expr,
            model_group=hl.struct(
                cpg=cpg_expr,
                **(additional_grouping_exprs or {}),
            ),
        )
    )


def _sum_agg_expr(
    fields_to_sum: Optional[List[str]] = None,
    exprs_to_sum: Optional[
        Union[hl.expr.StructExpression, Dict[str, hl.expr.NumericExpression]]
    ] = None,
    ht: Optional[hl.Table] = None,
) -> hl.expr.StructExpression:
    """
    Return an aggregation expression to sum fields or expressions in a Table.

    The aggregation expression is a struct with the sum or array_sum of the fields or
    expressions provided in `fields_to_sum` or `exprs_to_sum`.

    :param fields_to_sum: List of fields to sum in the Table. Default is None.
    :param exprs_to_sum: Dictionary of expressions to sum in the Table. Default is None.
    :param ht: Input Table. Default is None.
    :return: Aggregation expression to sum fields or expressions in the Table.
    """
    if fields_to_sum is None and exprs_to_sum is None:
        raise ValueError("Either 'fields_to_sum' or 'exprs_to_sum' must be provided.")
    if fields_to_sum is not None and ht is None:
        raise ValueError("ht must be provided if 'fields_to_sum' is provided.")

    exprs_to_sum = exprs_to_sum or {}
    exprs_to_sum = hl.struct(
        **exprs_to_sum,
        **{f: ht[f] for f in fields_to_sum or []}
    )

    return hl.struct(
        **{
            k: (
                hl.agg.array_sum(v)
                if isinstance(v, hl.ArrayExpression)
                else hl.agg.sum(v)
            )
            for k, v in exprs_to_sum.items()
        }
    )


# TODO: I have changed this so that it doesn't split up the pops anymore. I don't think
#  it is necessary to do this, and it makes the code more complicated. We should just
#  keep things the way we do for freq with a freq_meta. So we expect the
#  observed_variants and plateau_models_expr to be arrays of the same length.
def build_models(
    ht: hl.Table,
    coverage_expr: hl.expr.Int32Expression,
    weighted: bool = False,
    keys: Tuple[str] = (
        "context",
        "ref",
        "alt",
        "methylation_level",
    ),
    model_group_expr: Optional[hl.expr.StructExpression] = None,
    high_cov_definition: int = COVERAGE_CUTOFF,
    upper_cov_cutoff: Optional[int] = None,
    skip_coverage_model: bool = False,
    log10_coverage: bool = True,
    additional_grouping: Tuple[str] = (),
):
    """
    Build coverage and plateau models.

    This function builds models (plateau_models) using linear regression to calibrate
    mutation rate estimates against the proportion observed of each substitution,
    context, and methylation level in `ht`.

    Two plateau models are fit, one for CpG transitions, and one for the remainder of
    sites (transversions and non CpG transitions).

    The plateau models only consider high coverage sites, or sites above a median
    coverage of `high_cov_definition` and median coverage below `upper_cov_cutoff`.

    Plateau model: adjusts proportion of expected variation based on location in the
    genome and CpG status.

    The x and y of the plateau models:
        - x: `mu_snp` - mutation rate
        - y: proportion observed ('observed_variants' / 'possible_variants')

    This function also builds models (coverage models) to calibrate the proportion of
    expected variation at low coverage sites (sites below `high_cov_definition`).

    The coverage models are built by creating a scaling factor across all high coverage
    sites, applying this ratio to the low coverage sites, and running a linear
    regression.

    Coverage model: corrects proportion of expected variation at low coverage sites.
    Low coverage sites are defined as sites with median coverage < `high_cov_definition`.

    The x and y of the coverage model:

        - x: groupings of exome coverage at low coverage sites (log10 transformed if
          requested)
        - y: sum('observed_variants') / (`high_coverage_scale_factor` * sum('possible_variants' * 'mu_snp') at low coverage sites

    `high_coverage_scale_factor` = sum('observed_variants') / sum('possible_variants' * 'mu_snp') at high coverage sites

    .. note::

        This function expects that the input Table (`ht`) contains only high quality
        synonymous variants below 0.1% frequency.

    This function expects that the following fields are present in `coverage_ht`:

        - context: trinucleotide genomic context.
        - ref: the reference allele.
        - alt: the alternate allele.
        - methylation_level: methylation level.
        - mu_snp: mutation rate.
        - cpg: whether the variant is a CpG.
        - observed_variants: the number of observed variants in the dataset for each
          variant. Note that the term "variant" here refers to a specific substitution,
          context, methylation level, and coverage combination.
        - possible_variants: the number of possible variants in the dataset for each
          variant.

    :param ht: Input Table.
    :param coverage_expr: Expression that defines the coverage metric.
    :param weighted: Whether to weight the plateau models (a linear regression
        model) by 'possible_variants'. Default is False.
    :param keys: Annotations used to group observed and possible variant counts.
        Default is ("context", "ref", "alt", "methylation_level").
    :param model_group_expr: Optional Expression with `high_or_low_coverage` annotation
        to group variants into high or low coverage models. If not provided, the
        `calibration_model_group_expr()` function is used to define the grouping.
    :param high_cov_definition: Lower coverage cutoff. Sites with coverage above this
        cutoff are considered well covered. Default is `COVERAGE_CUTOFF`.
    :param upper_cov_cutoff: Upper coverage cutoff. Sites with coverage above this
        cutoff are excluded from the high coverage Table. Default is None.
    :param skip_coverage_model: Whether to skip generating the coverage model. If set
        to True, None is returned instead of the coverage model. Default is False.
    :param log10_coverage: Whether to convert coverage sites with log10 when building
        the coverage model. Default is True.
    :param additional_grouping: Optional Additional annotations to group by before
        counting the observed and possible variants. Default is ().
    :return: Coverage model and plateau models.
    """
    if model_group_expr is None:
        # Define whether the variant should be included in the high or low coverage
        # model.
        model_group_expr = calibration_model_group_expr(
            coverage_expr,
            ht.cpg,
            low_cov_cutoff=0,
            high_cov_cutoff=high_cov_definition,
            upper_cov_cutoff=upper_cov_cutoff,
            skip_coverage_model=skip_coverage_model,
            cpg_in_high_only=True,
        )

    grouping = keys + additional_grouping
    mu_type_fields = ("cpg", "transition", "mutation_type", "mutation_type_model")
    has_mu_type = all([x in ht.row for x in mu_type_fields])
    grouping += mu_type_fields if has_mu_type else ()

    grouping_exprs = {
        "build_model": model_group_expr.annotate(
            model_group=model_group_expr.model_group.select(
                *[k for k in model_group_expr.model_group if k != "sfs_bin"]
            )
        )
    }
    if not skip_coverage_model:
        grouping_exprs["exomes_coverage"] = hl.or_missing(
            model_group_expr.high_or_low_coverage == "low", coverage_expr
        )

    ht = ht.group_by(*grouping, **grouping_exprs).aggregate(
        mu_snp=hl.agg.take(ht.mu_snp, 1)[0],
        **_sum_agg_expr(
            fields_to_sum=["observed_variants", "possible_variants"], ht=ht
        ),
        sfs_bin_counts=hl.agg.group_by(
            model_group_expr.model_group,
            _sum_agg_expr(
                fields_to_sum=["observed_variants"], ht=ht
            ),
        ).items()
    ).key_by(*keys)

    ht = ht.explode("sfs_bin_counts")
    ht = ht.transmute(
        build_model=ht.build_model.annotate(model_group=ht.sfs_bin_counts[0]),
        observed_variants=hl.if_else(
            ht.sfs_bin_counts[0].sfs_bin == 0,
            ht.observed_variants,
            # ht.observed_variants.map(lambda x: ht.possible_variants-x),
            ht.sfs_bin_counts[1].observed_variants,
        ),
        possible_variants=ht.possible_variants,
    )

    if not has_mu_type:
        ht = annotate_mutation_type(ht)

    # Build plateau models.
    is_high_expr = ht.build_model.high_or_low_coverage == "high"
    agg_expr = {
        "plateau": hl.agg.filter(
            is_high_expr,
            build_plateau_models(
                ht.mu_snp,
                ht.observed_variants,
                ht.possible_variants,
                model_group_expr=ht.build_model.model_group,
                weighted=weighted,
            ),
        )
    }

    if not skip_coverage_model:
        # The coverage model is only built using the full dataset observed variants
        # so use the first element if the observed_variants is an array.
        obs_is_array = isinstance(ht.observed_variants, hl.expr.ArrayExpression)
        obs_expr = ht.observed_variants[0] if obs_is_array else ht.observed_variants

        # Create a metric that represents the relative mutability of the exome calculated
        # on high coverage sites and will be used as scaling factor when building the
        # coverage model.
        agg_expr["high_coverage_scale_factor"] = hl.agg.filter(
            is_high_expr,
            hl.agg.sum(obs_expr)
            / hl.agg.sum(ht.possible_variants * ht.mu_snp)
        )

        # Get the observed variant count and mu_snp for low coverage sites.
        agg_expr["coverage"] = hl.agg.filter(
            ht.build_model.high_or_low_coverage == "low",
            hl.agg.group_by(
                ht.exomes_coverage,
                hl.struct(
                    obs=hl.agg.sum(obs_expr),
                    mu_snp=hl.agg.sum(ht.possible_variants * ht.mu_snp)
                )
            )
        )

    models = ht.aggregate(hl.struct(**agg_expr), _localize=False)

    # Build coverage model.
    coverage_model = None
    if not skip_coverage_model:
        coverage_model = models.coverage.map_values(
            lambda x: x.annotate(
                low_coverage_oe=x.obs / (models.high_coverage_scale_factor * x.mu_snp)
            )
        )

        # TODO: consider weighting here as well.
        coverage_model = coverage_model.items().aggregate(
            lambda x: build_coverage_model(
                x[1].low_coverage_oe, x[0], log10_coverage=log10_coverage
            )
        ).beta

    return coverage_model, models.plateau


def build_plateau_models(
    mu_snp_expr: hl.expr.Float64Expression,
    observed_variants_expr: Union[hl.expr.ArrayExpression, hl.expr.Int64Expression],
    possible_variants_expr: Union[hl.expr.ArrayExpression, hl.expr.Int64Expression],
    weighted: bool = False,
    cpg_expr: Optional[hl.expr.BooleanExpression] = None,
    model_group_expr: Optional[hl.expr.StructExpression] = None,
) -> hl.expr.DictExpression:
    """
    Build plateau models to calibrate mutation rate to compute predicted proportion observed value.

    The x and y of the plateau models:

        - x: `mu_snp_expr`
        - y: `observed_variants_expr` / `possible_variants_expr`
          or `pops_observed_variants_array_expr`[index] / `possible_variants_expr`
          if `pops` is specified

    :param mu_snp_expr: Float64Expression of the mutation rate.
    :param observed_variants_expr: ArrayExpression or Int64Expression of the observed
        variant counts.
    :param possible_variants_expr: ArrayExpression or Int64Expression of the possible
        variant counts.
    :param weighted: Whether to generalize the model to weighted least squares using
        'possible_variants'. Default is False.
    :return: DictExpression of the intercepts and slopes of the plateau models for each
        group in `model_group_expr`.
    :param cpg_expr: BooleanExpression noting whether a site is a CPG site.
    :param model_group_expr: StructExpression to group by in the aggregation.

    """
    obs_is_array = isinstance(observed_variants_expr, hl.expr.ArrayExpression)
    pos_is_array = isinstance(possible_variants_expr, hl.expr.ArrayExpression)

    def _linreg(
        o: hl.expr.Int64Expression = observed_variants_expr,
        p: hl.expr.Int64Expression = possible_variants_expr,
    ) -> hl.expr.StructExpression:
        """
        Compute the linear regression of the observed variant counts over the possible variant counts.

        :param o: Int64Expression of the observed variant counts.
        :param p: Int64Expression of the possible variant counts.
        :return: StructExpression of the intercept and slope of the linear regression.
        """
        return hl.agg.linreg(
            o / p, [1, mu_snp_expr], weight=p if weighted else None
        ).beta

    if obs_is_array and pos_is_array:
        agg_expr = hl.agg.array_agg(
            lambda x: _linreg(*x),
            hl.zip(observed_variants_expr, possible_variants_expr)
        )
    elif obs_is_array:
        agg_expr = hl.agg.array_agg(lambda x: _linreg(o=x), observed_variants_expr)
    elif pos_is_array:
        agg_expr = hl.agg.array_agg(lambda x: _linreg(p=x), possible_variants_expr)
    else:
        agg_expr = _linreg()

    if model_group_expr is None and cpg_expr is None:
        return agg_expr

    model_group_expr = model_group_expr or hl.struct()
    if cpg_expr is not None:
        model_group_expr = model_group_expr.annotate(cpg=cpg_expr)

    return hl.agg.group_by(model_group_expr, agg_expr)


def build_coverage_model(
    low_coverage_oe_expr: hl.expr.Float64Expression,
    coverage_expr: hl.expr.Float64Expression,
    log10_coverage: bool = False,
) -> hl.expr.StructExpression:
    """
    Build coverage model.

    This function uses linear regression to build a model of coverage to correct
    proportion of expected variation at low coverage sites.

    The x and y of the coverage model:

        - x: `coverage_expr`
        - y: `low_coverage_oe_expr`

    :param low_coverage_oe_expr: The Float64Expression of observed:expected ratio
        for a given coverage level.
    :param coverage_expr: The Float64Expression of the coverage expression.
    :param log10_coverage: Whether to convert coverage sites by log10 when building the
        coverage model. Default is False.
    :return: StructExpression with intercept and slope of the model.
    """
    if log10_coverage:
        logger.info(
            "Converting coverage sites by log10 when building the coverage model."
        )
        coverage_expr = hl.log10(coverage_expr)

    return hl.agg.linreg(low_coverage_oe_expr, [1, coverage_expr])


def get_all_pop_lengths(
    ht: hl.Table,
    pops: Tuple[str],
    obs_expr: hl.expr.StructExpression,
) -> List[Tuple[str, str]]:
    """
    Get the minimum length of observed variant counts array for each population downsampling.

    The observed variant counts for each population in `pops` are specified by
    annotations on the `obs_expr` expression.

    The function also performs a check that arrays of variant counts within population
    downsamplings all have the same lengths.

    :param ht: Input Table containing `obs_expr`.
    :param pops: Populations used to categorize observed variant counts in downsamplings.
    :param obs_expr: Expression for the population observed variant counts. Should be a
        struct containing an array for each pop in `pops`.
    :return: A Dictionary with the minimum array length for each population.
    """
    # TODO: This function will be converted into doing just the length check if there
    #  is no usage of pop_lengths in the constraint pipeline.
    # Get minimum length of downsamplings for each population.
    pop_downsampling_lengths = ht.aggregate(
        [hl.agg.min(hl.len(obs_expr[pop])) for pop in pops]
    )

    # Zip population name with their downsampling length.
    pop_lengths = list(zip(pop_downsampling_lengths, pops))
    logger.info("Found: %s", "".join(map(str, pop_lengths)))

    assert ht.all(
        hl.all(
            lambda f: f,
            [hl.len(obs_expr[pop]) == length for length, pop in pop_lengths],
        )
    ), (
        "The arrays of variant counts within population downsamplings have different"
        " lengths!"
    )

    return pop_lengths


def get_constraint_grouping_expr(
    vep_annotation_expr: hl.StructExpression,
    coverage_expr: Optional[hl.Int32Expression] = None,
    include_transcript_group: bool = True,
    include_canonical_group: bool = True,
    include_mane_select_group: bool = False,
) -> Dict[str, Union[hl.StringExpression, hl.Int32Expression, hl.BooleanExpression]]:
    """
    Collect annotations used for constraint groupings.

    Function collects the following annotations:
        - annotation - 'most_severe_consequence' annotation in `vep_annotation_expr`
        - modifier - classic lof annotation from 'lof' annotation in
            `vep_annotation_expr`, LOFTEE annotation from 'lof' annotation in
            `vep_annotation_expr`, PolyPhen annotation from 'polyphen_prediction' in
            `vep_annotation_expr`, or "None" if neither is defined
        - gene - 'gene_symbol' annotation inside `vep_annotation_expr`
        - coverage - exome coverage if `coverage_expr` is specified
        - transcript - id from 'transcript_id' in `vep_annotation_expr` (added when
            `include_transcript_group` is True)
        - canonical from `vep_annotation_expr` (added when `include_canonical_group` is
            True)
        - mane_select from `vep_annotation_expr` (added when `include_mane_select_group` is
            True)

    .. note::
        This function expects that the following fields are present in
        `vep_annotation_expr`:
        - lof
        - polyphen_prediction
        - most_severe_consequence
        - gene_symbol
        - transcript_id (if `include_transcript_group` is True)
        - canonical (if `include_canonical_group` is True)
        - mane_select (if `include_mane_select_group` is True)

    :param vep_annotation_expr: StructExpression of VEP annotation.
    :param coverage_expr: Optional Int32Expression of exome coverage. Default is None.
    :param include_transcript_group: Whether to include the transcript annotation in the
        groupings. Default is True.
    :param include_canonical_group: Whether to include canonical annotation in the
        groupings. Default is True.
    :param include_mane_select_group: Whether to include mane_select annotation in the
        groupings. Default is False.
    :return: A dictionary with keys as annotation names and values as actual
        annotations.
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
        groupings["mane_select"] = hl.or_else(
            hl.is_defined(vep_annotation_expr.mane_select), False
        )

    return groupings


def annotate_exploded_vep_for_constraint_groupings(
    ht: hl.Table,
    coverage_expr: hl.expr.Int32Expression = None,
    vep_annotation: str = "transcript_consequences",
    include_canonical_group: bool = True,
    include_mane_select_group: bool = False,
) -> Tuple[Union[hl.Table, hl.MatrixTable], Tuple[str]]:
    """
    Annotate Table with annotations used for constraint groupings.

    Function explodes the specified VEP annotation (`vep_annotation`) and adds the following annotations:
        - annotation -'most_severe_consequence' annotation in `vep_annotation`
        - modifier - classic lof annotation from 'lof' annotation in
            `vep_annotation`, LOFTEE annotation from 'lof' annotation in
            `vep_annotation`, PolyPhen annotation from 'polyphen_prediction' in
            `vep_annotation`, or "None" if neither is defined
        - gene - 'gene_symbol' annotation inside `vep_annotation`
        - transcript - id from 'transcript_id' in `vep_annotation` (added when
            `include_transcript_group` is True)
        - canonical from `vep_annotation` (added when `include_canonical_group` is
            True)
        - mane_select from `vep_annotation` (added when `include_mane_select_group` is
            True)

    .. note::
        This function expects that the following annotations are present in `ht`:
        - vep
        - exome_coverage

    :param ht: Input Table or MatrixTable.
    :param coverage_expr: Expression that defines the coverage metric.
    :param vep_annotation: Name of annotation in 'vep' annotation (one of
        "transcript_consequences" and "worst_csq_by_gene") that will be used for
        obtaining constraint annotations. Default is "transcript_consequences".
    :param include_canonical_group: Whether to include 'canonical' annotation in the
        groupings. Default is True. Ignored unless `vep_annotation` is  "transcript_consequences".
    :param include_mane_select_group: Whether to include 'mane_select' annotation in the
        groupings. Default is False. Ignored unless `vep_annotation` is  "transcript_consequences".
    :return: A tuple of input Table or MatrixTable with grouping annotations added and
        the names of added annotations.
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


# TODO: Not totally sure this is needed anymore...
def compute_expected_variants(
    ht: hl.Table,
    plateau_models_expr: hl.StructExpression,
    mu_expr: hl.Float64Expression,
    cov_corr_expr: hl.Float64Expression,
    possible_variants_expr: hl.Int64Expression,
    cpg_expr: hl.BooleanExpression,
    pop: Optional[str] = None,
) -> Dict[str, Union[hl.Float64Expression, hl.Int64Expression]]:
    """
    Apply plateau models for all sites and for a population (if specified) to compute predicted proportion observed ratio and expected variant counts.

    :param ht: Input Table.
    :param plateau_models_expr: Linear models (output of `build_models()`, with the values
        of the dictionary formatted as a StructExpression of intercept and slope, that
        calibrates mutation rate to proportion observed for high coverage exomes. It
        includes models for CpG, non-CpG sites, and each population if specified.
    :param mu_expr: Float64Expression of mutation rate.
    :param possible_variants_expr: Int64Expression of possible variant counts.
    :param cov_corr_expr: Float64Expression of corrected coverage expression.
    :param cpg_expr: BooleanExpression noting whether a site is a CPG site.
    :param pop: Optional population to use when applying plateau model. Default is
        None.
    :return: A dictionary with predicted proportion observed ratio and expected variant
        counts.
    """
    if pop is None:
        pop = ""
        plateau_model = hl.literal(plateau_models_expr.total)[cpg_expr]
        slope = plateau_model[1]
        intercept = plateau_model[0]
        agg_func = hl.agg.sum
        ann_to_sum = ["observed_variants", "possible_variants"]
    else:
        plateau_model = hl.literal(plateau_models_expr[pop])
        slope = hl.map(lambda f: f[cpg_expr][1], plateau_model)
        intercept = hl.map(lambda f: f[cpg_expr][0], plateau_model)
        agg_func = hl.agg.array_sum
        pop = f"_{pop}"
        ann_to_sum = [f"downsampling_counts{pop}"]

    # Apply plateau models for specified population.
    ppo_expr = mu_expr * slope + intercept

    # Generate sum aggregators for 'predicted_proportion_observed' and
    # 'expected_variants', for specified population.
    agg_expr = {
        f"predicted_proportion_observed{pop}": agg_func(ppo_expr),
        f"expected_variants{pop}": agg_func(
            ppo_expr * cov_corr_expr * possible_variants_expr
        ),
    }

    # Generate sum aggregators for 'observed_variants' and 'possible_variants' on
    # the entire dataset if pop is None, and for 'downsampling_counts' for
    # specified population if pop is not None.
    agg_expr.update({ann: agg_func(ht[ann]) for ann in ann_to_sum})

    return agg_expr


# TODO: Can probably be modified some given my other changes.
def oe_aggregation_expr(
    ht: hl.Table,
    filter_expr: hl.expr.BooleanExpression,
    pops: Tuple[str] = (),
    exclude_mu_sum: bool = False,
) -> hl.expr.StructExpression:
    """
    Get aggregation expressions to compute the observed:expected ratio for rows defined by `filter_expr`.

    Return a Struct containing aggregation expressions to sum the number of observed
    variants, possible variants, expected variants, and mutation rate (if
    `exclude_mu_sum` is not True) for rows defined by `filter_expr`. The Struct also
    includes an aggregation expression for the observed:expected ratio.

    The following annotations are in the returned StructExpression:
        - obs - the sum of observed variants filtered to `filter_expr`.
        - mu - the sum of mutation rate of variants filtered to `filter_expr`.
        - possible - possible number of variants filtered to `filter_expr`.
        - exp - expected number of variants filtered to `filter_expr`.
        - oe - observed:expected ratio of variants filtered to `filter_expr`.

        If `pops` is specified:
            - gen_anc_exp - Struct with the expected number of variants per population (for
              all pop in `pops`) filtered to `filter_expr`.
            - gen_anc_obs - Struct with the observed number of variants per population (for
              all pop in `pops`) filtered to `filter_expr`.

    .. note::
        The following annotations should be present in `ht`:
            - observed_variants
            - mu
            - possible_variants
            - expected_variants
        If `pops` is specified, the following annotations should also be present:
            - expected_variants_{pop} for all pop in `pops`
            - downsampling_counts_{pop} for all pop in `pops`

    :param ht: Input Table to create observed:expected ratio aggregation expressions for.
    :param filter_expr: Boolean expression used to filter `ht` before aggregation.
    :param pops: List of populations to compute constraint metrics for. Default is ().
    :param exclude_mu_sum: Whether to exclude mu sum aggregation expression from
        returned struct. Default is False.
    :return: StructExpression with observed:expected ratio aggregation expressions.
    """
    # Create aggregators that sum the number of observed variants, possible variants,
    # and expected variants and compute observed:expected ratio.
    agg_expr = {
        "obs": hl.agg.sum(ht.observed_variants),
        "exp": hl.agg.sum(ht.expected_variants),
        "possible": hl.agg.sum(ht.possible_variants),
    }
    agg_expr["oe"] = divide_null(agg_expr["obs"], agg_expr["exp"])

    # Create an aggregator that sums the mutation rate.
    if not exclude_mu_sum:
        agg_expr["mu"] = hl.agg.sum(ht.mu)

    # Create aggregators that sum the number of observed variants
    # and expected variants for each population if pops is specified.
    if pops:
        agg_expr["gen_anc_exp"] = hl.struct(
            **{pop: hl.agg.array_sum(ht[f"expected_variants_{pop}"]) for pop in pops}
        )
        agg_expr["gen_anc_obs"] = hl.struct(
            **{pop: hl.agg.array_sum(ht[f"downsampling_counts_{pop}"]) for pop in pops}
        )

    agg_expr = hl.struct(**agg_expr)
    return hl.agg.group_by(filter_expr, agg_expr).get(True, hl.missing(agg_expr.dtype))


def apply_plateau_models(
    mu_expr: hl.Float64Expression,
    plateau_models_expr: hl.ArrayExpression,
) -> Union[hl.ArrayExpression, hl.Float64Expression]:
    """
    Compute the predicted probability observed.

    The predicted probability observed is computed as the mutation rate adjusted by the
    plateau model.

    :param mu_expr: Mutation rate expression.
    :param plateau_models_expr: This can be either a single plateau model, where the
        first element is the intercept and the second element is the slope, or an array
        of plateau models.
    :return: Predicted probability observed expression.
    """
    def _apply_model(plateau_model: hl.ArrayExpression) -> hl.Float64Expression:
        """
        Apply the plateau model to the mutation rate expression.

        :param plateau_model: ArrayExpression of the plateau model.
        :return: Predicted probability observed expression.
        """
        slope = plateau_model[1]
        intercept = plateau_model[0]
        ppo_expr = mu_expr * slope + intercept

        return ppo_expr

    if plateau_models_expr.dtype.element_type == hl.tarray(hl.tfloat64):
        return plateau_models_expr.map(lambda x: _apply_model(x))

    return _apply_model(plateau_models_expr)


def coverage_correction_expr(
    coverage_expr: hl.Float64Expression,
    coverage_model: Tuple[float, float],
    low_coverage_expr: Optional[hl.BooleanExpression] = None,
    coverage_cutoff: Optional[int] = None,
    log10_coverage: bool = False,
) -> hl.Float64Expression:
    """
    Compute the coverage correction expression.

    .. note::

        One and only one of `low_coverage_expr` or `coverage_cutoff` must be specified.

    The coverage correction expression is computed as follows:

        - If the coverage is 0, the coverage correction is 0.
        - If the low coverage expression (`low_coverage_expr`) is True, or the
          coverage (`coverage_expr`) is below the coverage cutoff (`coverage_cutoff`),
          the coverage model is applied to the coverage.
        - Otherwise, the coverage correction is 1.

    :param coverage_expr: Float64Expression of the coverage.
    :param coverage_model: Tuple of the intercept and slope of the coverage model.
    :param low_coverage_expr: Optional BooleanExpression indicating whether the site is
        a low coverage site, and the coverage model should be applied. Default is None.
    :param coverage_cutoff: Optional coverage cutoff. If specified, the coverage model
        is applied to sites with coverage below this cutoff. Default is None.
    :param log10_coverage: Whether to convert coverage sites by log10 when applying the
        coverage model. Default is False.
    :return: Float64Expression of the coverage correction.
    """
    if low_coverage_expr is None and coverage_cutoff is None:
        raise ValueError(
            "Either 'low_coverage_expr' or 'coverage_cutoff' must be specified!"
        )
    if low_coverage_expr is not None and coverage_cutoff is not None:
        raise ValueError(
            "Only one of 'low_coverage_expr' or 'coverage_cutoff' can be specified!"
        )

    if coverage_cutoff is not None:
        low_coverage_expr = coverage_expr < coverage_cutoff

    if log10_coverage:
        cov_corr_expr = hl.log10(coverage_expr)
    else:
        cov_corr_expr = coverage_expr

    return (
        hl.case()
        .when(coverage_expr == 0, 0)
        .when(
            low_coverage_expr,
            coverage_model[1] * cov_corr_expr + coverage_model[0],
        )
        .default(1)
    )


def apply_models(
    mu_expr: hl.expr.Float64Expression,
    plateau_models_expr: hl.ArrayExpression,
    possible_variants_expr: hl.expr.Int64Expression,
    coverage_model: Optional[Tuple[float, float]] = None,
    coverage_expr: Optional[hl.expr.Int32Expression] = None,
    cpg_expr: Optional[hl.expr.BooleanExpression] = None,
    model_group_expr: Optional[hl.expr.StructExpression] = None,
    high_cov_definition: int = COVERAGE_CUTOFF,
    log10_coverage: bool = True,
) -> hl.expr.StructExpression:
    """
    Apply calibration models to compute predicted proportion observed ratio and expected variant counts.

    :param mu_expr: Mutation rate expression.
    :param plateau_models_expr: This can be either a single plateau model, where the
        first element is the intercept and the second element is the slope, or an array
        of plateau models.
    :param possible_variants_expr: Int64Expression of possible variant counts to
        multiply the predicted probability observed by. Default is None.
    :param coverage_model: Tuple of the intercept and slope of the coverage model.
        Default is None.
    :param coverage_expr: Optional Int32Expression of the coverage. Default is None.
    :param cpg_expr: Optional BooleanExpression of whether the variant is a CpG site.
        Default is None.
    :param model_group_expr: Optional Expression with `high_or_low_coverage` annotation
        to group variants into high or low coverage models. If not provided, the
        `calibration_model_group_expr()` function is used to define the grouping.
        Default is None.
    :param high_cov_definition: Definition of high coverage. This is used to define
        whether the site/group of sites is high or low coverage. Default is
        `COVERAGE_CUTOFF`.
    :param log10_coverage: Whether to convert coverage sites by log10 when applying the
        coverage model. Default is True.
    :return: A dictionary with predicted proportion observed ratio and expected variant
        counts.
    """
    if coverage_model is not None and coverage_expr is None:
        raise ValueError(
            "If 'coverage_model' is specified, 'coverage_expr' must also be specified!"
        )
    if model_group_expr is None and (coverage_expr is None or cpg_expr is None):
        raise ValueError(
            "If 'model_group_expr' is not specified, 'coverage_expr' and 'cpg_expr' must"
            " be specified!"
        )

    if model_group_expr is None:
        # Get the annotations relevant for applying the calibration models.
        model_group_expr = calibration_model_group_expr(
            coverage_expr,
            cpg_expr,
            high_cov_cutoff=high_cov_definition,
            skip_coverage_model=coverage_model is None,
        )

    # Apply plateau models.
    ppo_expr = apply_plateau_models(mu_expr, plateau_models_expr)
    apply_expr = hl.struct(
        mu=mu_expr * possible_variants_expr,
        predicted_proportion_observed=ppo_expr,
        expected_variants=ppo_expr * possible_variants_expr,
    )

    # Get the coverage correction expression if a coverage model is provided.
    if coverage_model is not None:
        cov_corr_expr = coverage_correction_expr(
            coverage_expr,
            coverage_model,
            low_coverage_expr=model_group_expr.high_or_low_coverage == "low",
            log10_coverage=log10_coverage,
        )
        apply_expr = apply_expr.annotate(
            expected_variants=apply_expr.expected_variants * cov_corr_expr,
            coverage_correction=cov_corr_expr,
        )

    return apply_expr


def aggregate_expected_variants_expr(
    ht: hl.Table,
    additional_fields_to_sum: Optional[List[str]] = None,
    additional_exprs_to_sum: Optional[Dict[str, hl.expr.Expression]] = None,
) -> hl.expr.StructExpression:
    """
    Get an aggregation expression for the sum of expected variants and other fields.

    An aggregate sum or array sum is created for the following fields:

        - mu_snp
        - observed_variants
        - possible_variants
        - predicted_proportion_observed
        - coverage_correction
        - expected_variants

    :param ht: Input Table.
    :param additional_fields_to_sum: List of additional fields in `ht` to get an
        aggregate sum expression for. Default is None.
    :param additional_exprs_to_sum: Dictionary of additional expressions in `ht` to get
        an aggregate sum expression for. Field names are the keys and expressions are
        the values. Default is None.
    :return: StructExpression with the sum of expected variants and other fields.
    """
    return _sum_agg_expr(
        fields_to_sum=[
            "mu_snp",
            "observed_variants",
            "possible_variants",
            "predicted_proportion_observed",
            "coverage_correction",
            "expected_variants",
        ] + (additional_fields_to_sum or []),
        exprs_to_sum=additional_exprs_to_sum,
        ht=ht,
    )


def compute_pli(
    ht: hl.Table,
    obs_expr: hl.expr.Int64Expression,
    exp_expr: hl.expr.Float64Expression,
    expected_values: Optional[Dict[str, float]] = None,
    min_diff_convergence: float = 0.001,
) -> hl.StructExpression:
    """
    Compute the pLI score using the observed and expected variant counts.

    Full details on pLI can be found in the ExAC paper: Lek, M., Karczewski, K.,
    Minikel, E. et al. Analysis of protein-coding genetic variation in 60,706 humans.
    Nature 536, 285–291 (2016).

    pLI is the probability of being loss-of-function intolerant, and this function
    computes that probability using the expectation-maximization (EM) algorithm.

    We assume a 3 state model, where each gene fits into one of three categories
    with respect loss-of-function variation sensitivity:

        - Null: where protein truncating variation is completely tolerated by natural
          selection.
        - Recessive (Rec): where heterozygous pLoFs are tolerated but homozygous pLoFs
          are not.
        - Haploinsufficient (LI): where heterozygous pLoFs are not tolerated.

    The function requires the expected amount of loss-of-function depletion for each of
    these states. The default provided is based on the observed depletion of
    protein-truncating variation in the Blekhman autosomal recessive and ClinGen
    dosage sensitivity gene sets (Supplementary Information Table 12 of the above
    reference):

        - Null: 1.0, assume tolerant genes have the expected amount of truncating
          variation.
        - Rec: 0.463, derived from the empirical mean observed/expected rate of
          truncating variation for recessive disease genes (0.463).
        - LI: 0.089, derived from the empirical mean observed/expected rate of
          truncating variation for severe haploinsufficient genes.

    The output StructExpression will include the following annotations:

        - pLI: Probability of loss-of-function intolerance; probability that transcript
          falls into distribution of haploinsufficient genes.
        - pNull: Probability that transcript falls into distribution of unconstrained
          genes.
        - pRec: Probability that transcript falls into distribution of recessive genes.

    :param ht: Input Table containing `obs_expr` and `exp_expr`.
    :param obs_expr: Expression for the number of observed variants on each gene or
        transcript in `ht`.
    :param exp_expr: Expression for the number of expected variants on each gene or
        transcript in `ht`.
    :param expected_values: Dictionary containing the expected values for 'Null',
        'Rec', and 'LI' to use as starting values.
    :param min_diff_convergence: Minimum iteration change in LI to consider the EM
        model convergence criteria as met. Default is 0.001.
    :return: StructExpression for pLI scores.
    """
    if expected_values is None:
        expected_values = {"Null": 1.0, "Rec": 0.463, "LI": 0.089}

    # Set up initial values.
    last_pi = {k: 0 for k in expected_values.keys()}
    pi = {k: 1 / len(expected_values.keys()) for k in expected_values.keys()}

    dpois_expr = {
        k: hl.or_missing(
            exp_expr > 0, hl.dpois(obs_expr, exp_expr * expected_values[k])
        )
        for k in pi
    }
    _ht = ht.select(dpois=dpois_expr)
    # Checkpoint the temp HT because it will need to be aggregated several times.
    _ht = _ht.checkpoint(new_temp_file(prefix="compute_pli", extension="ht"))

    # Calculate pLI scores.
    while abs(pi["LI"] - last_pi["LI"]) > min_diff_convergence:
        last_pi = copy.deepcopy(pi)
        pi_expr = {k: v * _ht.dpois[k] for k, v in pi.items()}
        row_sum_expr = hl.sum([pi_expr[k] for k in pi])
        pi_expr = {k: pi_expr[k] / row_sum_expr for k, v in pi.items()}
        pi = _ht.aggregate({k: hl.agg.mean(pi_expr[k]) for k in pi.keys()})

    # Get expression for pLI scores.
    pli_expr = {k: v * dpois_expr[k] for k, v in pi.items()}
    row_sum_expr = hl.sum([pli_expr[k] for k in pi])

    return hl.struct(**{f"p{k}": pli_expr[k] / row_sum_expr for k in pi.keys()})


def oe_confidence_interval(
    obs_expr: hl.expr.Int64Expression,
    exp_expr: hl.expr.Float64Expression,
    alpha: float = 0.05,
) -> hl.expr.StructExpression:
    """
    Determine the confidence interval around the observed:expected ratio.

    For a given pair of observed (`obs_expr`) and expected (`exp_expr`) values, the
    function computes the density of the Poisson distribution (performed using Hail's
    `dpois` module) with fixed k (`x` in `dpois` is set to the observed number of
    variants) over a range of lambda (`lamb` in `dpois`) values, which are given by the
    expected number of variants times a varying parameter ranging between 0 and 2 (the
    observed:expected ratio is typically between 0 and 1, so we want to extend the
    upper bound of the confidence interval to capture this). The cumulative density
    function of the Poisson distribution density is computed and the value of the
    varying parameter is extracted at points corresponding to `alpha` (defaults to 5%)
    and 1-`alpha` (defaults to 95%) to indicate the lower and upper bounds of the
    confidence interval.

    The following annotations are in the output StructExpression:
        - lower - the lower bound of confidence interval
        - upper - the upper bound of confidence interval

    :param obs_expr: Expression for the observed variant counts of pLoF, missense, or
        synonymous variants in `ht`.
    :param exp_expr: Expression for the expected variant counts of pLoF, missense, or
        synonymous variants in `ht`.
    :param alpha: The significance level used to compute the confidence interval.
        Default is 0.05.
    :return: StructExpression for the confidence interval lower and upper bounds.
    """
    # Set up range between 0 and 2.
    range_expr = hl.range(0, 2000).map(lambda x: hl.float64(x) / 1000)
    range_dpois_expr = range_expr.map(lambda x: hl.dpois(obs_expr, exp_expr * x))

    # Compute cumulative density function of the Poisson distribution density.
    cumulative_dpois_expr = hl.cumulative_sum(range_dpois_expr)
    max_cumulative_dpois_expr = cumulative_dpois_expr[-1]
    norm_dpois_expr = cumulative_dpois_expr.map(lambda x: x / max_cumulative_dpois_expr)

    # Extract the value of the varying parameter within specified range.
    lower_idx_expr = hl.argmax(
        norm_dpois_expr.map(lambda x: hl.or_missing(x < alpha, x))
    )
    upper_idx_expr = hl.argmin(
        norm_dpois_expr.map(lambda x: hl.or_missing(x > 1 - alpha, x))
    )
    return hl.struct(
        lower=hl.if_else(obs_expr > 0, range_expr[lower_idx_expr], 0),
        upper=range_expr[upper_idx_expr],
    )


def calculate_raw_z_score(
    obs_expr: hl.expr.Int64Expression,
    exp_expr: hl.expr.Float64Expression,
) -> hl.expr.Float64Expression:
    """
    Compute the signed raw z-score using observed and expected variant counts.

    The raw z-scores are positive when the transcript had fewer variants than expected,
    and are negative when transcripts had more variants than expected.

    :param obs_expr: Observed variant count expression.
    :param exp_expr: Expected variant count expression.
    :return: Raw z-score expression.
    """
    chisq_expr = divide_null((obs_expr - exp_expr) ** 2, exp_expr)
    return hl.sqrt(chisq_expr) * hl.if_else(obs_expr > exp_expr, -1, 1)


def get_constraint_flags(
    exp_expr: hl.expr.Float64Expression,
    raw_z_expr: hl.expr.Float64Expression,
    raw_z_lower_threshold: Optional[Union[float, hl.expr.Float64Expression]] = -5.0,
    raw_z_upper_threshold: Optional[Union[float, hl.expr.Float64Expression]] = 5.0,
    flag_postfix: str = "",
) -> Dict[str, hl.expr.Expression]:
    """
    Determine the constraint flags that define why constraint will not be calculated.

    Flags which are added:
        - "no_exp_{flag_postfix}" - for genes that have missing or zero expected variants.
        - "outlier_{flag_postfix}" - for genes that are raw z-score outliers:
          (`raw_z_expr` < `raw_z_lower_threshold`) or (`raw_z_expr` >
          `raw_z_upper_threshold`).

    :param exp_expr: Expression for the expected variant counts of pLoF, missense, or
        synonymous variants.
    :param raw_z_expr: Expression for the signed raw z-score of pLoF, missense, or
        synonymous variants.
    :param raw_z_lower_threshold: Lower threshold for the raw z-score. When `raw_z_expr`
        is less than this threshold it is considered an 'outlier'. Default is -5.0.
    :param raw_z_upper_threshold: Upper threshold for the raw z-score. When `raw_z_expr`
        is greater than this threshold it is considered an 'outlier'. Default is 5.0.
    :param flag_postfix: Postfix to add to the end of the constraint flag names.
    :return: Dictionary containing expressions for constraint flags.
    """
    outlier_expr = False
    if raw_z_lower_threshold is not None:
        outlier_expr |= hl.or_else(raw_z_expr < raw_z_lower_threshold, False)
    if raw_z_upper_threshold is not None:
        outlier_expr |= hl.or_else(raw_z_expr > raw_z_upper_threshold, False)

    if flag_postfix:
        flag_postfix = f"_{flag_postfix}"

    constraint_flags = {
        f"no_exp{flag_postfix}": hl.or_else(exp_expr <= 0, True),
        f"outlier{flag_postfix}": hl.or_else(outlier_expr, False),
    }

    return constraint_flags


def calculate_raw_z_score_sd(
    raw_z_expr: hl.expr.Float64Expression,
    flag_expr: hl.expr.StringExpression,
    mirror_neg_raw_z: bool = True,
) -> hl.expr.Expression:
    """
    Calculate the standard deviation of the raw z-score.

    When using `mirror_neg_raw_z` is True, all the negative raw z-scores (defined by
    `raw_z_expr`) are combined with those same z-scores multiplied by -1 (to create a
    mirrored distribution).

    :param raw_z_expr: Expression for the raw z-score.
    :param flag_expr: Expression for the constraint flags. z-score will not be
        calculated if flags are present.
    :param mirror_neg_raw_z: Whether the standard deviation should be computed using a
        mirrored distribution of negative `raw_z_expr`.
    :return: StructExpression containing standard deviation of the raw z-score and
        the z-score.
    """
    filter_expr = (hl.len(flag_expr) == 0) & hl.is_defined(raw_z_expr)

    if mirror_neg_raw_z:
        filter_expr &= raw_z_expr < 0
        sd_expr = hl.agg.explode(
            lambda x: hl.agg.stats(x), [raw_z_expr, -raw_z_expr]
        ).stdev
    else:
        sd_expr = hl.agg.stats(raw_z_expr).stdev

    return hl.agg.filter(filter_expr, sd_expr)


def add_gencode_transcript_annotations(
    ht: hl.Table,
    gencode_ht: hl.Table,
    annotations: List[str] = ["level", "transcript_type"],
) -> hl.Table:
    """
    Add GENCODE annotations to Table based on transcript id.

    .. note::
        Added annotations by default are:
        - level
        - transcript_type

        Computed annotations are:
        - chromosome
        - cds_length
        - num_coding_exons

    :param ht: Input Table.
    :param gencode_ht: Table with GENCODE annotations.
    :param annotations: List of GENCODE annotations to add. Default is ["level", "transcript_type"].
        Added annotations also become keys for the group by when computing "cds_length" and "num_coding_exons".
    :return: Table with transcript annotations from GENCODE added.
    """
    gencode_ht = gencode_ht.annotate(
        length=gencode_ht.interval.end.position
        - gencode_ht.interval.start.position
        + 1,
        chromosome=gencode_ht.interval.start.contig,
    )

    # Obtain CDS annotations from GENCODE file and calculate CDS length and
    # number of exons.
    annotations_to_add = set(annotations + ["chromosome", "transcript_id", "length"])

    gencode_cds = (
        gencode_ht.filter(gencode_ht.feature == "CDS")
        .select(*annotations_to_add)
        .key_by("transcript_id")
        .drop("interval")
    )

    annotations_to_add.remove("length")

    gencode_cds = (
        gencode_cds.group_by(*annotations_to_add)
        .aggregate(
            cds_length=hl.agg.sum(gencode_cds.length), num_coding_exons=hl.agg.count()
        )
        .key_by("transcript_id")
    )

    gencode_cds = gencode_cds.checkpoint(
        new_temp_file(prefix="gencode_cds", extension="ht")
    )

    # Add GENCODE annotations to input Table.
    ht = ht.annotate(**gencode_cds[ht.transcript])

    return ht
