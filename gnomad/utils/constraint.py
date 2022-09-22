# noqa: D100


import logging
from typing import Any, Optional, Tuple, Union

import hail as hl

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_utils")
logger.setLevel(logging.INFO)


def annotate_with_mu(
    ht: hl.Table,
    mutation_ht: hl.Table,
    mu_annotation: str = "mu_snp",
) -> hl.Table:
    """
    Annotate SNP mutation rate for the input Table.

    .. note::
        `ht` is expected to include annotations that `mutation_ht` is keyed by, but these annotations don't need to be the keys of `ht`.

    :param ht: Input Table to annotate.
    :param mutation_ht: Mutation rate Table.
    :param mu_annotation: The name of mutation rate annotation in `mutation_ht`. Defaults to 'mu_snp'.
    :return: Table with mutational rate annotation added (default name for annotation is 'mu_snp').
    """
    mu = mutation_ht.index(*[ht[k] for k in mutation_ht.key])[mu_annotation]
    return ht.annotate(
        **{mu_annotation: hl.case().when(hl.is_defined(mu), mu).or_error("Missing mu")}
    )


def count_variants_by_group(
    ht: hl.Table,
    freq_expr: hl.expr.ArrayExpression = None,
    freq_meta_expr: hl.expr.ArrayExpression = None,
    count_singletons: bool = False,
    count_downsamplings: Optional[Tuple[str]] = (),
    additional_grouping: Optional[Tuple[str]] = (),
    partition_hint: int = 100,
    omit_methylation: bool = False,
    use_table_group_by: bool = False,
    singleton_expr: hl.expr.BooleanExpression = None,
    max_af: Optional[float] = None,
) -> Union[hl.Table, Any]:
    """
    Count number of observed or possible variants by context, ref, alt, and optionally methylation_level.

    Function counts variants in `freq_expr` based on specified criteria (singleton or in downsamplings or None), and aggregates variant counts based on
    'context', 'ref', 'alt', 'methylation_level' (optional), and `additional_grouping` provided.

    If `freq_expr` is not given and variants in downsamplings or singleton variants need to be counted, `freq_expr` defaults to `ht.freq`. If variants in downsamplings
    needs to be counted and `freq_meta_expr` is None, `freq_meta_expr` defaults to `ht.freq_meta`. If variants in singleton needs to be counted and singleton_expr is None,
    `singleton_expr` defaults to `freq_expr[0].AC == 1`.

    Function will return a Table with annotations used for grouping ('context', 'ref', 'alt', 'methylation_level' (optional), `additional_grouping`) and 'variant_count' annotation.

    .. note::
        The following annotations should be present in `ht`:
        - ref - the middle base of `context`
        - alt - the alternate base
        - context - trinucleotide genomic context
        - methylation_level - methylation level
        - freq - Allele frequency information (AC, AN, AF, homozygote count) in gnomAD release.
        - freq_meta - gnomAD frequency metadata. An ordered list containing the frequency aggregation group for each element of the gnomad_freq array row annotation.

    :param ht: Input Hail Table.
    :param freq_expr: ArrayExpression of Structs with with 'AC' and 'AF' annotations. If `freq_expr` is None and any of `count_downsamplings`, `max_af`, and `count_singletons` is True, `freq_expr` would be `ht.freq`.
    :param freq_meta_expr: ArrayExpression of meta dictionaries corresponding to freq_expr. If `count_downsamplings` and `freq_meta_expr` is None, `freq_meta_expr` would be `ht.freq_meta`.
    :param count_singletons: Whether to count singletons. Defaults to False.
    :param count_downsamplings: List of populations to use for downsampling counts. Defaults to ().
    :param additional_grouping: Additional features to group by. e.g. 'exome_coverage'. Defaults to ().
    :param partition_hint: Target number of partitions for aggregation. Defaults to 100.
    :param omit_methylation: Whether to omit 'methylation_level' from the grouping when counting variants. Defaults to False.
    :param use_table_group_by: Whether to group the variant counts before aggregate the variant counts. If `use_table_group_by` is False, function will return a hl.StructExpression. Defaults to False.
    :param singleton_expression: Expression for defining a singleton. When `count_singletons` is True and `singleton_expression` is None, `singleton_expression` would be `freq_expr[0].AC == 1`. Defaults to None.
    :param max_af: Maximum variant AF to keep. By default, no AF cutoff is applied.
    :return: Table including 'variant_count' and downsampling counts if requested.
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
            "freq_meta_expr was not provided, using 'freq_meta' as the frequency metadata annotation."
        )
        freq_meta_expr = ht.freq_meta
    if count_singletons and singleton_expr is None:
        logger.warning(
            "count_singletons is True and singleton_expr was not provided, using freq_expr[0].AC == 1 as the singleton expression."
        )
        # Slower, but more flexible (allows for downsampling agg's)
        singleton_expr = freq_expr[0].AC == 1

    grouping = hl.struct(context=ht.context, ref=ht.ref, alt=ht.alt)
    if not omit_methylation:
        grouping = grouping.annotate(methylation_level=ht.methylation_level)
    for group in additional_grouping:
        grouping = grouping.annotate(**{group: ht[group]})

    agg = {
        "variant_count": hl.agg.count_where(freq_expr[0].AF <= max_af)
        if max_af
        else hl.agg.count()
    }
    if count_singletons:
        agg["singleton_count"] = hl.agg.count_where(singleton_expr)

    for pop in count_downsamplings:
        agg[f"downsampling_counts_{pop}"] = downsampling_counts_expr(
            freq_expr, freq_meta_expr, pop, max_af=max_af
        )
        if count_singletons:
            agg[f"singleton_downsampling_counts_{pop}"] = downsampling_counts_expr(
                freq_expr, freq_meta_expr, pop, singleton=True
            )

    if use_table_group_by:
        return ht.group_by(**grouping).partition_hint(partition_hint).aggregate(**agg)
    else:
        return ht.aggregate(
            hl.struct(**{field: hl.agg.group_by(grouping, agg[field]) for field in agg})
        )


def downsampling_counts_expr(
    freq_expr: hl.expr.ArrayExpression,
    freq_meta_expr: hl.expr.ArrayExpression,
    pop: str = "global",
    variant_quality: str = "adj",
    singleton: bool = False,
    max_af: Optional[float] = None,
) -> hl.expr.ArrayExpression:
    """
    Count variants in `freq_expr` of downsamplings with specified criteria in their meta data.

    The meta data (`freq_meta_expr`) of each downsampling should have 'group',
    'pop', and 'downsampling' keys, where 'group' key is 'variant_quality' and 'pop' key is 'pop'.

    :param freq_expr: ArrayExpression of Structs with with 'AC' and 'AF' annotations.
    :param freq_meta_expr: ArrayExpression of meta dictionaries corresponding to freq_expr.
    :param pop: Population. Defaults to 'global'.
    :param variant_quality: Variant quality for "group" key in `freq_meta_expr`. Defaults to 'adj'.
    :param singleton: Whether to sum only alleles that are singletons. Defaults to False.
    :param max_af: Maximum variant allele frequency to keep. By default no allele frequency cutoff is applied.
    :return: Variant counts in downsamplings for specified population.
    """
    # get indices of dictionaries in meta dictionaries that only have the "downsampling" key with required "group" and "pop" keys
    indices = hl.enumerate(freq_meta_expr).filter(
        lambda f: (f[1].size() == 3)
        & (f[1].get("group") == variant_quality)
        & (f[1].get("pop") == pop)
        & f[1].contains("downsampling")
    )
    # get an array of indices sorted by "downsampling" key
    sorted_indices = hl.sorted(indices, key=lambda f: hl.int(f[1]["downsampling"])).map(
        lambda x: x[0]
    )

    def _get_criteria(i: hl.expr.Int32Expression) -> hl.expr.ArrayNumericExpression:
        if singleton:
            return hl.int(freq_expr[i].AC == 1)
        elif max_af:
            return hl.int((freq_expr[i].AC > 0) & (freq_expr[i].AF <= max_af))
        else:
            return hl.int(freq_expr[i].AC > 0)

    return hl.agg.array_sum(hl.map(_get_criteria, sorted_indices))


def annotate_mutation_type(
    t: Union[hl.MatrixTable, hl.Table]
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Annotate mutation types.

    The following annotations are added to the output Table:
        - cpg
        - transition
        - mutation_type - one of "CpG", "non-CpG transition", or "transversion"
        - mutation_type_model

    :param t: Input Table or MatrixTable.
    :return: Table with mutation type annotations added.
    """
    # Determine the middle index of context by sampling the first 100 values of 'context'
    context_lengths = list(filter(None, set(hl.len(t.context).take(100))))
    if len(context_lengths) > 1:
        raise ValueError(
            "More than one length was found among the first 100 'context' values. Length of 'context' should be consistent."
        )
    else:
        context_length = context_lengths[0]
        logger.info("Detected a length of %d for context length", context_length)

    if context_length == 3:
        mid_index = 1
    elif context_length == 7:
        mid_index = 3
    else:
        raise ValueError(
            f"The length of context should be either 3 or 7, instead of {context_length}."
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
        hl.case()
        .when(t.cpg, "CpG")
        .when(t.transition, "non-CpG transition")
        .default("transversion")
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
    t: Union[hl.MatrixTable, hl.Table]
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
    t: Union[hl.Table, hl.MatrixTable]
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Return the deduplicated context by collapsing DNA strands.

    Function returns the reverse complement for 'ref, 'alt', and 'context' if the reference allele is either 'G' or 'T'.

    The following annotations are added to the output Table:
        - was_flipped - whether the 'ref, 'alt', and 'context' were flipped (reverse complement taken)

    :param ht: Input Table.
    :return: Table with deduplicated context annotation (ref, alt, context, was_flipped).
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
