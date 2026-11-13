"""
Functions to annotate the context Table for constraint.

Covers mutation rates, mutation types, trinucleotide context, methylation levels, and
GERP cutoffs.
"""

import logging
from typing import Callable, Dict, Optional, Tuple, Union

import hail as hl

from gnomad.utils.constraint.counts import _resolve_row_annotation_expr
from gnomad.utils.reference_genome import get_primary_contigs, get_reference_genome
from gnomad.utils.vep import add_most_severe_consequence_to_consequence

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_utils")
logger.setLevel(logging.INFO)


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
        ht,
        [hl.parse_locus_interval(c, ref.name) for c in get_primary_contigs(ref.name)],
    )

    # Add annotations for 'ref' and 'alt'.
    ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1])

    # Trim heptamer context to create trimer context.
    ht = trimer_from_heptamer(ht)
    ht = ht.filter(ht.context.matches(f"[ATCG]{{{3}}}"))

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
        "sift_score",
        "polyphen_score",
        "domains",
        "uniprot_isoform",
        "amino_acids",
        "codons",
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
                    data_type: (
                        ann_ht.select_globals(
                            *[x for x in g if x in ann_ht.globals]
                        ).index_globals()
                    )
                    for data_type, ann_ht in hts.items()
                }
            )
            for k, (hts, g) in global_anns.items()
        }
    )

    return ht


def calculate_gerp_cutoffs(
    ht: hl.Table,
    gerp_expr: Optional[hl.expr.Float64Expression] = None,
    lower_percentile: float = 0.05,
    upper_percentile: float = 0.95,
) -> Tuple[float, float]:
    """
    Find GERP cutoffs at the given percentile thresholds.

    .. note::

        Uses ``hl.agg.approx_quantiles``, so results are approximate.

    :param ht: Input Table.
    :param gerp_expr: GERP score expression. Default is ``ht.gerp``.
    :param lower_percentile: Lower percentile threshold (0-1). Default is 0.05.
    :param upper_percentile: Upper percentile threshold (0-1). Default is 0.95.
    :return: Tuple of (lower cutoff, upper cutoff) GERP scores.
    """
    gerp_expr = _resolve_row_annotation_expr(ht, "gerp", gerp_expr, "gerp_expr")

    cutoffs = ht.aggregate(
        hl.agg.approx_quantiles(gerp_expr, [lower_percentile, upper_percentile])
    )
    return cutoffs[0], cutoffs[1]
