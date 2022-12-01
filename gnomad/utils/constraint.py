"""Script containing generic constraint functions that may be used in the constraint pipeline."""

import copy
import logging
from typing import Any, Dict, List, Optional, Tuple, Union

import hail as hl

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_utils")
logger.setLevel(logging.INFO)

COVERAGE_CUTOFF = 40
"""
Minimum median exome coverage differentiating high coverage sites from low coverage sites.

Low coverage sites require an extra calibration when computing the proportion of expected variation.
"""


def annotate_with_mu(
    ht: hl.Table,
    mutation_ht: hl.Table,
    mu_annotation: str = "mu_snp",
) -> hl.Table:
    """
    Annotate SNP mutation rate for the input Table.

    .. note::

        Function expects that`ht` includes`mutation_ht`'s key fields. Note that these
        annotations don't need to be the keys of `ht`.

    :param ht: Input Table to annotate.
    :param mutation_ht: Mutation rate Table.
    :param mu_annotation: The name of mutation rate annotation in `mutation_ht`.
        Default is 'mu_snp'.
    :return: Table with mutational rate annotation added.
    """
    mu = mutation_ht.index(*[ht[k] for k in mutation_ht.key])[mu_annotation]
    return ht.annotate(
        **{mu_annotation: hl.case().when(hl.is_defined(mu), mu).or_error("Missing mu")}
    )


def count_variants_by_group(
    ht: hl.Table,
    freq_expr: Optional[hl.expr.ArrayExpression] = None,
    freq_meta_expr: Optional[hl.expr.ArrayExpression] = None,
    count_singletons: bool = False,
    count_downsamplings: Tuple[str] = (),
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
            freq_expr, freq_meta_expr, pop, max_af=max_af
        )
        if count_singletons:
            logger.info(
                "Counting singleton variants in downsamplings for population '%s', and"
                " adding as 'singleton_downsampling_counts_%s' annotation.",
                pop,
                pop,
            )
            agg[f"singleton_downsampling_counts_{pop}"] = downsampling_counts_expr(
                freq_expr, freq_meta_expr, pop, singleton=True
            )
    # Apply each variant count aggregation in `agg` to get counts for all
    # combinations of `grouping`.
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
    Return an aggregation expression to compute an array of counts of all downsamplings found in `freq_expr` where specified criteria is met.

    The frequency metadata (`freq_meta_expr`) should be in a similar format to the
    `freq_meta` annotation added by `annotate_freq()`. Each downsampling should have
    'group', 'pop', and 'downsampling' keys. Included downsamplings are those where
    'group' == `variant_quality` and 'pop' == `pop`.

    :param freq_expr: ArrayExpression of Structs with 'AC' and 'AF' annotations.
    :param freq_meta_expr: ArrayExpression containing the set of groupings for each
        element of the `freq_expr` array (e.g., [{'group': 'adj'}, {'group': 'adj',
        'pop': 'nfe'}, {'downsampling': '5000', 'group': 'adj', 'pop': 'global'}]).
    :param pop: Population to use for filtering by the 'pop' key in `freq_meta_expr`.
        Default is 'global'.
    :param variant_quality: Variant quality to use for filtering by the 'group' key in
        `freq_meta_expr`. Default is 'adj'.
    :param singleton: Whether to filter to only singletons before counting (AC == 1).
        Default is False.
    :param max_af: Maximum variant allele frequency to keep. By default no allele
        frequency cutoff is applied.
    :return: Aggregation Expression for an array of the variant counts in downsamplings
        for specified population.
    """
    # Get indices of dictionaries in meta dictionaries that only have the
    # "downsampling" key with specified "group" and "pop" values.
    indices = hl.enumerate(freq_meta_expr).filter(
        lambda f: (f[1].size() == 3)
        & (f[1].get("group") == variant_quality)
        & (f[1].get("pop") == pop)
        & f[1].contains("downsampling")
    )
    # Get an array of indices sorted by "downsampling" key.
    sorted_indices = hl.sorted(indices, key=lambda f: hl.int(f[1]["downsampling"])).map(
        lambda x: x[0]
    )

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

    ..note:

        This function uses the term 'mutation_type' because 'variant_type' is already
        used in this repo to indicate a variant's multiallelic and SNP/indel status.

    :param t: Input Table or MatrixTable.
    :return: Table with mutation type annotations added.
    """
    # Determine the context length by collecting all the context lengths.
    context_lengths = list(filter(None, set(hl.len(t.context).collect())))
    if len(context_lengths) > 1:
        raise ValueError(
            "More than one length was found among the first 100 'context' values."
            " Length of 'context' should be consistent."
        )
    else:
        context_length = context_lengths[0]
        logger.info("Detected a length of %d for context length", context_length)
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


def build_models(
    coverage_ht: hl.Table,
    weighted: bool = False,
    pops: Tuple[str] = (),
    keys: Tuple[str] = (
        "context",
        "ref",
        "alt",
        "methylation_level",
        "mu_snp",
    ),
    cov_cutoff: int = COVERAGE_CUTOFF,
) -> Tuple[Tuple[float, float], hl.expr.StructExpression]:
    """
    Build coverage and plateau models.

    This function builds models (plateau_models) using linear regression to calibrate
    mutation rate estimates against the proportion observed of each substitution,
    context, and methylation level in `coverage_ht`.

    Two plateau models are fit, one for CpG transitions, and one for the remainder of
    sites (transversions and non CpG transitions).

    The plateau models only consider high coverage sites, or sites above a median
    coverage of `cov_cutoff`.

    Plateau model: adjusts proportion of expected variation based on location in the
    genome and CpG status.
    The x and y of the plateau models:
    - x: `mu_snp` - mutation rate
    - y: proportion observed ('observed_variants' or 'observed_{pop}' / 'possible_variants')

    This function also builds models (coverage models) to calibrate the proportion of
    expected variation at low coverage sites (sites below `cov_cutoff`).

    The coverage models are built by creating a scaling factor across all high coverage
    sites, applying this ratio to the low coverage sites, and running a linear
    regression.

    Coverage model: corrects proportion of expected variation at low coverage sites.
    Low coverage sites are defined as sites with median coverage < `cov_cutoff`.

    The x and y of the coverage model:
    - x: log10 groupings of exome coverage at low coverage sites
    - y: sum('observed_variants')/ (`high_coverage_scale_factor` * sum('possible_variants' * 'mu_snp') at low coverage sites

    `high_coverage_scale_factor` = sum('observed_variants') /
                        sum('possible_variants' * 'mu_snp') at high coverage sites

    .. note::

        This function expects that the input Table(`coverage_ht`) was created using
        `get_proportion_observed_by_coverage`, which means that `coverage_ht` should
        contain only high quality synonymous variants below 0.1% frequency.

        This function also expects that the following fields are present in
        `coverage_ht`:
        - context - trinucleotide genomic context
        - ref - the reference allele
        - alt - the alternate allele
        - methylation_level - methylation level
        - cpg - whether the site is CpG site
        - exome_coverage - median exome coverage at integer values between 1-100
        - observed_variants - the number of observed variants in the dataset for each
        variant. Note that the term "variant" here refers to a specific substitution,
        context, methylation level, and coverage combination
        - downsampling_counts_{pop} (optional) - array of observed variant counts per
        population after downsampling. Used only when `pops` is specified.
        - mu_snp - mutation rate
        - possible_variants - the number of possible variants in the dataset for each
        variant

    :param coverage_ht: Input coverage Table.
    :param weighted: Whether to weight the plateau models (a linear regression
        model) by 'possible_variants'. Default is False.
    :param pops: List of populations used to build plateau models.
        Default is ().
    :param keys: Annotations used to group observed and possible variant counts.
        Default is ("context", "ref", "alt", "methylation_level", "mu_snp").
    :param cov_cutoff: Median coverage cutoff. Sites with coverage above this cutoff
        are considered well covered. Default is `COVERAGE_CUTOFF`.
    :return: Coverage model and plateau models.
    """
    # Filter to sites with coverage above `cov_cutoff`.
    high_cov_ht = coverage_ht.filter(coverage_ht.exome_coverage >= cov_cutoff)
    agg_expr = {
        "observed_variants": hl.agg.sum(high_cov_ht.observed_variants),
        "possible_variants": hl.agg.sum(high_cov_ht.possible_variants),
    }
    for pop in pops:
        agg_expr[f"observed_{pop}"] = hl.agg.array_sum(
            high_cov_ht[f"downsampling_counts_{pop}"]
        )

    # Generate a Table with all necessary annotations (x and y listed above)
    # for the plateau models.
    high_cov_group_ht = high_cov_ht.group_by(*keys).aggregate(**agg_expr)
    high_cov_group_ht = annotate_mutation_type(high_cov_group_ht)

    # Build plateau models.
    plateau_models_agg_expr = build_plateau_models(
        cpg_expr=high_cov_group_ht.cpg,
        mu_snp_expr=high_cov_group_ht.mu_snp,
        observed_variants_expr=high_cov_group_ht.observed_variants,
        possible_variants_expr=high_cov_group_ht.possible_variants,
        pops_observed_variants_array_expr=[
            high_cov_group_ht[f"observed_{pop}"] for pop in pops
        ],
        weighted=weighted,
    )
    if pops:
        # Map the models to their corresponding populations if pops is specified.
        _plateau_models = dict(
            high_cov_group_ht.aggregate(hl.struct(**plateau_models_agg_expr))
        )
        pop_models = _plateau_models["pop"]
        plateau_models = {
            pop: hl.literal(pop_models[idx]) for idx, pop in enumerate(pops)
        }
        plateau_models["total"] = _plateau_models["total"]
        plateau_models = hl.struct(**plateau_models)
    else:
        plateau_models = high_cov_group_ht.aggregate(
            hl.struct(**plateau_models_agg_expr)
        )

    # Filter to sites with coverage below `cov_cutoff` and larger than 0.
    low_cov_ht = coverage_ht.filter(
        (coverage_ht.exome_coverage < cov_cutoff) & (coverage_ht.exome_coverage > 0)
    )

    # Create a metric that represents the relative mutability of the exome calculated
    # on high coverage sites and will be used as scaling factor when building the
    # coverage model.
    high_coverage_scale_factor = high_cov_ht.aggregate(
        hl.agg.sum(high_cov_ht.observed_variants)
        / hl.agg.sum(high_cov_ht.possible_variants * high_cov_ht.mu_snp)
    )

    # Generate a Table with all necessary annotations (x and y listed above)
    # for the coverage model.
    low_cov_group_ht = low_cov_ht.group_by(
        log_coverage=hl.log10(low_cov_ht.exome_coverage)
    ).aggregate(
        low_coverage_oe=hl.agg.sum(low_cov_ht.observed_variants)
        / (
            high_coverage_scale_factor
            * hl.agg.sum(low_cov_ht.possible_variants * low_cov_ht.mu_snp)
        )
    )

    # Build the coverage model.
    # TODO: consider weighting here as well
    coverage_model_expr = build_coverage_model(
        low_coverage_oe_expr=low_cov_group_ht.low_coverage_oe,
        log_coverage_expr=low_cov_group_ht.log_coverage,
    )
    coverage_model = tuple(low_cov_group_ht.aggregate(coverage_model_expr).beta)
    return coverage_model, plateau_models


def build_plateau_models(
    cpg_expr: hl.expr.BooleanExpression,
    mu_snp_expr: hl.expr.Float64Expression,
    observed_variants_expr: hl.expr.Int64Expression,
    possible_variants_expr: hl.expr.Int64Expression,
    pops_observed_variants_array_expr: List[hl.expr.ArrayExpression] = [],
    weighted: bool = False,
) -> Dict[str, Union[Dict[bool, hl.expr.ArrayExpression], hl.ArrayExpression]]:
    """
    Build plateau models to calibrate mutation rate to compute predicted proportion observed value.

    The x and y of the plateau models:
    - x: `mu_snp_expr`
    - y: `observed_variants_expr` / `possible_variants_expr`
    or `pops_observed_variants_array_expr`[index] / `possible_variants_expr`
    if `pops` is specified

    :param cpg_expr: BooleanExpression noting whether a site is a CPG site.
    :param mu_snp_expr: Float64Expression of the mutation rate.
    :param observed_variants_expr: Int64Expression of the observed variant counts.
    :param possible_variants_expr: Int64Expression of the possible variant counts.
    :param pops_observed_variants_array_expr: Nested ArrayExpression with all observed
        variant counts ArrayNumericExpressions for specified populations. e.g., `[[1,1,
        1],[1,1,1]]`. Default is None.
    :param weighted: Whether to generalize the model to weighted least squares using
        'possible_variants'. Default is False.
    :return: A dictionary of intercepts and slopes of plateau models. The keys are
        'total' (for all sites) and 'pop' (optional; for populations). The values for
        'total' is a dictionary (e.g., <DictExpression of type dict<bool,
        array<float64>>>), and the value for 'pop' is a nested list of dictionaries (e.
        g., <ArrayExpression of type array<array<dict<bool, array<float64>>>>>). The
        key of the dictionary in the nested list is CpG status (BooleanExpression), and
        the value is an ArrayExpression containing intercept and slope values.
    """
    # Build plateau models for all sites
    plateau_models_agg_expr = {
        "total": hl.agg.group_by(
            cpg_expr,
            hl.agg.linreg(
                observed_variants_expr / possible_variants_expr,
                [1, mu_snp_expr],
                weight=possible_variants_expr if weighted else None,
            ).beta,
        )
    }
    if pops_observed_variants_array_expr:
        # Build plateau models using sites in population downsamplings if
        # population is specified.
        plateau_models_agg_expr["pop"] = hl.agg.array_agg(
            lambda pop_obs_var_array_expr: hl.agg.array_agg(
                lambda pop_observed_variants: hl.agg.group_by(
                    cpg_expr,
                    hl.agg.linreg(
                        pop_observed_variants / possible_variants_expr,
                        [1, mu_snp_expr],
                        weight=possible_variants_expr,
                    ).beta,
                ),
                pop_obs_var_array_expr,
            ),
            pops_observed_variants_array_expr,
        )
    return plateau_models_agg_expr


def build_coverage_model(
    low_coverage_oe_expr: hl.expr.Float64Expression,
    log_coverage_expr: hl.expr.Float64Expression,
) -> hl.expr.StructExpression:
    """
    Build coverage model.

    This function uses linear regression to build a model of log10(coverage) to correct
    proportion of expected variation at low coverage sites.

    The x and y of the coverage model:
    - x: `log_coverage_expr`
    - y: `low_coverage_oe_expr`

    :param low_coverage_oe_expr: The Float64Expression of observed:expected ratio
        for a given coverage level.
    :param log_coverage_expr: The Float64Expression of log10 coverage.
    :return: StructExpression with intercept and slope of the model.
    """
    return hl.agg.linreg(low_coverage_oe_expr, [1, log_coverage_expr])


def get_all_pop_lengths(
    ht: hl.Table,
    pops: Tuple[str],
    prefix: str = "observed_",
) -> List[Tuple[str, str]]:
    """
    Get the minimum length of observed variant counts array for each population downsamping.

    The annotations are specified by the combination of `prefix` and each population in
    `pops`.

    :param ht: Input Table used to build population plateau models.
    :param pops: Populations used to categorize observed variant counts in downsampings.
    :param prefix: Prefix of population observed variant counts. Default is `observed_`.
    :return: A Dictionary with the minimum array length for each population.
    """
    # TODO: This function will be converted into doing just the length check if there
    # is no usage of pop_lengths in the constraint pipeline.
    # Get minimum length of downsamplings for each population.
    pop_downsampling_lengths = ht.aggregate(
        [hl.agg.min(hl.len(ht[f"{prefix}{pop}"])) for pop in pops]
    )

    # Zip population name with their downsampling length.
    pop_lengths = list(zip(pop_downsampling_lengths, pops))
    logger.info("Found: %s", "".join(map(str, pop_lengths)))

    assert ht.all(
        hl.all(
            lambda f: f,
            [hl.len(ht[f"{prefix}{pop}"]) == length for length, pop in pop_lengths],
        )
    ), (
        "The arrays of variant counts within population downsamplings have different"
        " lengths!"
    )

    return pop_lengths


def compute_oe_per_transcript(
    ht: hl.Table,
    annotation_name: str,
    keys: Tuple[str] = ("gene", "transcript", "canonical"),
    pops: Tuple[str] = (),
) -> hl.Table:
    """
    Compute the observed:expected ratio for pLoF variants, missense variants, or synonymous variants.

    Function sums the number of observed variants, possible variants, and expected
    variants across all the combinations of `keys`, and uses the expected variant
    counts and observed variant counts to compute the observed:expected ratio.

    The following annotations are added to the output Table:
        - obs_{annotation_name} - the sum of observed variants grouped by `keys`
        - mu_{annotation_name} - the sum of mutation rate grouped by `keys`
        - possible_{annotation_name} - possible number of variants grouped by `keys`
        - exp_{annotation_name} - expected number of variants grouped by `keys`
        - exp_{annotation_name}_{pop} (pop defaults to `POPS`) - expected number of
            variants per population grouped by `keys`
        - obs_{annotation_name}_{pop} (pop defaults to `POPS`) - observed number of
            variants per population grouped by `keys`
        - oe_{annotation_name} - observed:expected ratio

    .. note::
        The following annotations should be present in `ht`:
            - observed_variants
            - mu
            - possible_variants
            - expected_variants

    :param ht: Input Table with pLoF variants, missense variants, or synonymous
        variants.
    :param annotation_name: Annotation name used for constraint metrics to distinguish
        mutation types.
    :param keys: The keys of the output Table. Default is ('gene', 'transcript',
        'canonical').
    :param pops: List of populations used to compute constraint metrics. Default is ().
    :return: Table with observed:expected ratio.
    """
    # Create aggregators that sum the number of observed variants, possible variants,
    # and expected variants and compute observed:expected ratio.
    agg_expr = {
        f"obs_{annotation_name}": hl.agg.sum(ht.observed_variants),
        f"exp_{annotation_name}": hl.agg.sum(ht.expected_variants),
        f"oe_{annotation_name}": hl.agg.sum(ht.observed_variants)
        / hl.agg.sum(ht.expected_variants),
        f"possible_{annotation_name}": hl.agg.sum(ht.possible_variants),
    }

    # Create an aggregator that sums the mutation rate.
    if annotation_name != "mis_pphen":
        agg_expr[f"mu_{annotation_name}"] = hl.agg.sum(ht.mu)

    # Create aggregators that sum the number of observed variants
    # and expected variants for each population if pops` is specified.
    for pop in pops:
        agg_expr[f"exp_{annotation_name}_{pop}"] = hl.agg.array_sum(
            ht[f"expected_variants_{pop}"]
        )
        agg_expr[f"obs_{annotation_name}_{pop}"] = hl.agg.array_sum(
            ht[f"downsampling_counts_{pop}"]
        )
    return ht.group_by(*keys).aggregate(**agg_expr)


def compute_all_pLI_scores(
    ht: hl.Table,
    annotation_name: str,
    pops: Tuple[str] = (),
    keys: Tuple[str] = ("gene", "transcript", "canonical"),
    calculate_pop_pLI=False,
    n_ds_to_skip: int = 8,
) -> hl.Table:
    """
    Compute the pLI scores for pLoF variants.

    :param ht: Input Table with pLoF variants.
    :param annotation_name: Annotation name used for constraint metrics to distinguish
        mutation types.
    :param pops: Populations that need to compute pLI scores.
    :param keys: The keys of the output Table. Default is ('gene', 'transcript',
        'canonical').
    :param calculate_pop_pLI: Whether to compute the pLI scores for each population.
        Default is False.
    :param n_ds_to_skip: The number of small downsamplings to skip when calculating pLI scores for populations. Default is 8.
    :return: Table with pLI, pNull, and pRec scores.
    """
    # Filter to only variants with expected variant counts larger than 0.
    ht = ht.filter(ht[f"exp_{annotation_name}"] > 0)

    # Calculate pLI scores for each population if specified.
    if calculate_pop_pLI:
        pop_lengths = get_all_pop_lengths(ht, pops, f"obs_{annotation_name}_{pop}")
        for pop_length, pop in pop_lengths:
            logger.info(f"Calculating pLI for {pop}...")
            plis = []
            for i in range(n_ds_to_skip, pop_length):
                logger.info(i)
                ht = ht.filter(ht[f"exp_{annotation_name}_{pop}"][i] > 0)
                pli_ht = pLI(
                    ht,
                    ht[f"obs_{annotation_name}_{pop}"][i],
                    ht[f"exp_{annotation_name}_{pop}"][i],
                )
                plis.append(pli_ht[ht.key])
            ht = ht.annotate(
                **{
                    f"pLI_{pop}": [pli.pLI for pli in plis],
                    f"pRec_{pop}": [pli.pRec for pli in plis],
                    f"pNull_{pop}": [pli.pNull for pli in plis],
                }
            )
    # Calculate pLI scores for all variants.
    pli_result = pLI(ht, ht[f"obs_{annotation_name}"], ht[f"exp_{annotation_name}"])[
        ht.key
    ]
    pli_result = pli_result.rename(
        {x: f"{x}_{annotation_name}" for x in list(pli_result.keys())}
    )
    return ht.annotate(**pli_result).key_by(*keys)


def pLI(
    ht: hl.Table, obs: hl.expr.Int32Expression, exp: hl.expr.Float32Expression
) -> hl.Table:
    """
    Compute the pLI score using the observed and expected variant counts.

    The output Table will include the following annotations:
        - pLI - Probability of loss-of-function intolerance; probability that transcript falls into
            distribution of haploinsufficient genes
        - pNull - Probability that transcript falls into distribution of unconstrained genes
        - pRec - Probability that transcript falls into distribution of recessive genes

    :param ht: Input Table.
    :param obs: Expression for the number of observed variants on each gene or transcript in `ht`.
    :param exp: Expression for the number of expected variants on each gene or transcript in `ht`.
    :return: Table with pLI scores.
    """
    # Set up initial values.
    last_pi = {"Null": 0, "Rec": 0, "LI": 0}
    pi = {"Null": 1 / 3, "Rec": 1 / 3, "LI": 1 / 3}
    expected_values = {"Null": 1, "Rec": 0.463, "LI": 0.089}
    ht = ht.annotate(_obs=obs, _exp=exp)

    # Calculate pLI scores.
    while abs(pi["LI"] - last_pi["LI"]) > 0.001:
        last_pi = copy.deepcopy(pi)
        ht = ht.annotate(
            **{
                k: v * hl.dpois(ht._obs, ht._exp * expected_values[k])
                for k, v in pi.items()
            }
        )
        ht = ht.annotate(row_sum=hl.sum([ht[k] for k in pi]))
        ht = ht.annotate(**{k: ht[k] / ht.row_sum for k, v in pi.items()})
        pi = ht.aggregate({k: hl.agg.mean(ht[k]) for k in pi.keys()})

    # Annotate pLI scores.
    ht = ht.annotate(
        **{
            k: v * hl.dpois(ht._obs, ht._exp * expected_values[k])
            for k, v in pi.items()
        }
    )
    ht = ht.annotate(row_sum=hl.sum([ht[k] for k in pi]))
    return ht.select(**{f"p{k}": ht[k] / ht.row_sum for k in pi.keys()})


def oe_confidence_interval(
    ht: hl.Table,
    obs: hl.expr.Int32Expression,
    exp: hl.expr.Float32Expression,
    prefix: str = "oe",
    alpha: float = 0.05,
    select_only_ci_metrics: bool = True,
) -> hl.Table:
    """
    Determine the confidence interval around the observed:expected ratio.

    For a given pair of observed (`obs`) and expected (`exp`) values, the function 
    computes the density of the Poisson distribution (performed using Hail's `dpois` 
    module) with fixed k (`x` in `dpois` is set to the observed number of variants) 
    over a range of lambda (`lamb` in `dpois`) values, which are given by the expected 
    number of variants times a varying parameter ranging between 0 and 2. The 
    cumulative density function of the Poisson distribution density is computed and the 
    value of the varying parameter is extracted at points corresponding to `alpha` 
    (defaults to 5%) and 1-`alpha`(defaults to 95%) to indicate the lower and upper 
    bounds of the confidence interval.

    Function will have following annotations in the output Table in addition to keys:
    - {prefix}_lower - the lower bound of confidence interval
    - {prefix}_upper - the upper bound of confidence interval

    :param ht: Input Table with the observed and expected variant counts for either
        pLoF, missense, or synonymous variants.
    :param obs: Expression for the observed variant counts of pLoF, missense, or
        synonymous variants in `ht`.
    :param exp: Expression for the expected variant counts of pLoF, missense, or
        synonymous variants in `ht`.
    :param prefix: Prefix of upper and lower bounds, defaults to 'oe'.
    :param alpha: The significance level used to compute the confidence interval.
        Default is 0.05.
    :param select_only_ci_metrics: Whether to return only upper and lower bounds
        instead of keeping all the annotations except `_exp`, defaults to True.
    :return: Table with the confidence interval lower and upper bounds.
    """
    ht = ht.annotate(_obs=obs, _exp=exp)
    oe_ht = ht.annotate(_range=hl.range(0, 2000).map(lambda x: hl.float64(x) / 1000))
    oe_ht = oe_ht.annotate(
        _range_dpois=oe_ht._range.map(lambda x: hl.dpois(oe_ht._obs, oe_ht._exp * x))
    )

    oe_ht = oe_ht.transmute(_cumulative_dpois=hl.cumulative_sum(oe_ht._range_dpois))
    max_cumulative_dpois = oe_ht._cumulative_dpois[-1]
    oe_ht = oe_ht.transmute(
        _norm_dpois=oe_ht._cumulative_dpois.map(lambda x: x / max_cumulative_dpois)
    )
    oe_ht = oe_ht.transmute(
        _lower_idx=hl.argmax(
            oe_ht._norm_dpois.map(lambda x: hl.or_missing(x < alpha, x))
        ),
        _upper_idx=hl.argmin(
            oe_ht._norm_dpois.map(lambda x: hl.or_missing(x > 1 - alpha, x))
        ),
    )
    oe_ht = oe_ht.transmute(
        **{
            f"{prefix}_lower": hl.if_else(
                oe_ht._obs > 0, oe_ht._range[oe_ht._lower_idx], 0
            ),
            f"{prefix}_upper": oe_ht._range[oe_ht._upper_idx],
        }
    )

    return oe_ht.select(f"{prefix}_lower", f"{prefix}_upper")


def calculate_z(
    input_ht: hl.Table,
    obs: hl.expr.Int32Expression,
    exp: hl.expr.Float32Expression,
    z_score_output_annotation: str = "z_raw",
) -> hl.Table:
    """
    Compute the signed raw z score using observed and expected variant counts.

    The raw z scores are positive when the transcript had fewer variants than expected, and are negative when transcripts had more variants than expected.

    The following annotation is included in the output Table in addition to the `input_ht` keys:
        - `output` - the raw z score

    :param input_ht: Input Table.
    :param obs: Observed variant count expression.
    :param exp: Expected variant count expression.
    :param z_score_output_annotation: The annotation label to use for the raw z score output, defaults to 'z_raw'.
    :return: Table with raw z scores.
    """
    ht = input_ht.select(_obs=obs, _exp=exp)
    ht = ht.annotate(_chisq=(ht._obs - ht._exp) ** 2 / ht._exp)
    return ht.select(
        **{
            z_score_output_annotation: hl.sqrt(ht._chisq)
            * hl.if_else(ht._obs > ht._exp, -1, 1)
        }
    )


def calculate_z_scores(ht: hl.Table, mutation_type: str) -> hl.Table:
    """
    Calculate z scores for synomynous variants, missense variants, or pLoF variants.

    z score = {variant_annotation}_z_raw / {variant_annotation}_sd (variant_annotation could be syn, mis, or lof)

    Function will add the following annotations to output Table:
        - {mutation_type}_sd (global) - standard deviation of raw z score
        - constraint_flag - Reason gene does not have constraint metrics. One of:
            no variants: Zero observed synonymous, missense, pLoF variants
            no_exp_syn: Zero expected synonymous variants
            no_exp_mis: Zero expected missense variants
            no_exp_lof: Zero expected pLoF variants
            syn_outlier: Too many or too few synonymous variants; synonymous z score < -5 or synonymous z score > 5
            mis_too_many: Too many missense variants; missense z score < -5
            lof_too_many: Too many pLoF variants; pLoF z score < -5
        - syn_z - z score of synonymous variants
        - mis_z - z score of missense variants
        - lof_z - z score of pLoF variants

    :param ht: Input Table with observed and expected variant counts for one of
        synomynous variants, missense variants, and pLoF variants.
    :param mutation_type: The mutation type of variants in the input Table, one of
        'lof', 'mis', and 'syn'.
    :return: Table with z scores.
    """
    # Calculate z scores.
    ht = ht.annotate(
        **calculate_z(
            ht,
            ht[f"obs_{mutation_type}"],
            ht[f"exp_{mutation_type}"],
            f"{mutation_type}_z_raw",
        )[ht.key]
    )

    # Create constraint flag to annotate if the variant is defined, has expected variant count, and is outlier.
    reasons = hl.empty_set(hl.tstr)

    reasons = hl.if_else(
        hl.or_else(ht[f"obs_{mutation_type}"], 0) == 0,
        reasons.add("no_variants"),
        reasons,
    )

    reasons = hl.if_else(
        ht[f"exp_{mutation_type}"] > 0,
        reasons,
        reasons.add(f"no_exp_{mutation_type}"),
        missing_false=True,
    )

    # Compute the standard deviation after filtering to variant that is defined, is not
    # outlier, has expected variants, has raw z score, and its raw z score is less than
    # 0 for LoF and missense variants.
    if mutation_type != "syn":
        reasons = hl.if_else(
            ht[f"{mutation_type}_z_raw"] < -5,
            reasons.add(f"{mutation_type}_outlier"),
            reasons,
            missing_false=True,
        )
        ht = ht.annotate(constraint_flag=reasons)
        sds = ht.aggregate(
            hl.struct(
                **{
                    f"{mutation_type}_sd": hl.agg.filter(
                        ~hl.len(ht.constraint_flag)
                        == 0
                        & hl.is_defined(ht[f"{mutation_type}_z_raw"])
                        & (ht[f"{mutation_type}_z_raw"] < 0),
                        hl.agg.explode(
                            lambda x: hl.agg.stats(x),
                            [
                                ht[f"{mutation_type}_z_raw"],
                                -ht[f"{mutation_type}_z_raw"],
                            ],
                        ),
                    ).stdev
                }
            )
        )
    else:
        reasons = hl.if_else(
            hl.abs(ht.syn_z_raw) > 5,
            reasons.add("syn_outlier"),
            reasons,
            missing_false=True,
        )
        ht = ht.annotate(constraint_flag=reasons)
        sds = ht.aggregate(
            hl.struct(
                syn_sd=hl.agg.filter(
                    ~hl.len(ht.constraint_flag) == 0 & hl.is_defined(ht.syn_z_raw),
                    hl.agg.stats(ht.syn_z_raw),
                ).stdev,
            )
        )
    logger.info(sds)
    ht = ht.annotate_globals(**sds)
    return ht.transmute(
        **{
            f"{mutation_type}_z": ht[f"{mutation_type}_z_raw"]
            / sds[f"{mutation_type}_sd"]
        }
    )
