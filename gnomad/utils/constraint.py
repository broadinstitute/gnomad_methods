# noqa: D100

import logging
from typing import Dict, List, Tuple, Union

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

        This function uses the term 'mutation_type' because 'variant_type' is already used in this repo to indicate a variant's multiallelic and SNP/indel status.

    :param t: Input Table or MatrixTable.
    :return: Table with mutation type annotations added.
    """
    # Determine the middle index of context by collecting all the context lengths
    context_lengths = list(filter(None, set(hl.len(t.context).collect())))
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

    This function builds models (plateau_models) using linear regression to calibrate mutation rate
    estimates against the proportion observed of each substitution, context, and methylation level in
    `coverage_ht`.

    Two plateau models are fit, one for CpG transitions, and one for the remainder of sites
    (transversions and non CpG transitions).

    The plateau models only consider high coverage sites, or sites above a median coverage of
    `cov_cutoff`.

    Plateau model: adjusts proportion of expected variation based on location in the genome and CpG status.
    The x and y of the plateau models:
        - x: `mu_snp` - mutation rate
        - y: proportion observed ('observed_variants' or 'observed_{pop}' / 'possible_variants')

    This function also builds models (coverage models) to calibrate the proportion of expected variation
    at low coverage sites (sites below `cov_cutoff`).

    The coverage models are built by creating a scaling factor across all high coverage sites,
    applying this ratio to the low coverage sites, and running a linear regression.

    Coverage model: corrects proportion of expected variation at low coverage sites.
    Low coverage sites are defined as sites with median coverage < `cov_cutoff`.

    The x and y of the coverage model:
        - x: log10 groupings of exome coverage at low coverage sites
        - y: sum('observed_variants')/ (`high_coverage_scale_factor` * sum('possible_variants' * 'mu_snp') at low coverage sites

        `high_coverage_scale_factor` = sum('observed_variants') / sum('possible_variants' * 'mu_snp') at high coverage sites


    .. note::

        This function expects that the input Table(`coverage_ht`) was created using
        `get_proportion_observed_by_coverage`, which means that `coverage_ht` should contain only high quality
        synonymous variants below 0.1% frequency.

        This function also expects that the following fields are present in `coverage_ht`:
        - context - trinucleotide genomic context
        - ref - the reference allele
        - alt - the alternate base
        - methylation_level - methylation level
        - cpg - whether the site is CpG site
        - exome_coverage - median exome coverage at integer values between 1-100
        - observed_variants - the number of observed variants in the dataset for each variant. Note that the term
        "variant" here refers to a specific substitution, context, methylation level, and coverage combination.
        - downsampling_counts_{pop} (optional) - array of observed variant counts per population after downsampling when `pops` is specified.
        - mu_snp - mutation rate
        - possible_variants - the number of possible variants in the dataset for each variant.

    :param coverage_ht: Input coverage Table.
    :param weighted: Whether to weight the high coverage model (a linear regression model) by 'possible_variants'. Defaults to False.
    :param pops: List of populations used to build coverage and plateau models. Defaults to ().
    :param keys: Keys to group observed and possible variant counts.
        Default is ("context","ref", "alt", "methylation_level", "mu_snp").
    :param cov_cutoff: Median coverage cutoff. Sites with coverage above this cutoff are considered well covered
    and will be used to build plateau models. Sites below this cutoff have low coverage and will be used to build
    coverage models. Defaults to `COVERAGE_CUTOFF`.
    :return: Coverage model and plateau models.
    """
    # Filter to sites with coverage above `cov_cutoff`
    all_high_coverage_ht = coverage_ht.filter(coverage_ht.exome_coverage >= cov_cutoff)
    agg_expr = {
        "observed_variants": hl.agg.sum(all_high_coverage_ht.observed_variants),
        "possible_variants": hl.agg.sum(all_high_coverage_ht.possible_variants),
    }
    for pop in pops:
        agg_expr[f"observed_{pop}"] = hl.agg.array_sum(
            all_high_coverage_ht[f"downsampling_counts_{pop}"]
        )

    # Generate a Table with all necessary annotations (x and y listed above) for the plateau models
    high_coverage_ht = all_high_coverage_ht.group_by(*keys).aggregate(**agg_expr)
    high_coverage_ht = annotate_mutation_type(high_coverage_ht)

    # Build plateau models
    plateau_models_agg_expr = build_plateau_models(
        cpg_expr=high_coverage_ht.cpg,
        mu_snp_expr=high_coverage_ht.mu_snp,
        observed_variants_expr=high_coverage_ht.observed_variants,
        possible_variants_expr=high_coverage_ht.possible_variants,
        pop_observed_variants_exprs={
            pop: high_coverage_ht[f"observed_{pop}"] for pop in pops
        },
        weighted=weighted,
    )
    plateau_models = high_coverage_ht.aggregate(hl.struct(**plateau_models_agg_expr))

    # Filter to sites with coverage below `cov_cutoff` and larger than 0
    all_low_coverage_ht = coverage_ht.filter(
        (coverage_ht.exome_coverage < cov_cutoff) & (coverage_ht.exome_coverage > 0)
    )

    # Metric that represents the relative mutability of the exome calculated on high coverage sites
    # This metric is used as scaling factor when building the coverage model
    high_coverage_scale_factor = all_high_coverage_ht.aggregate(
        hl.agg.sum(all_high_coverage_ht.observed_variants)
        / hl.agg.sum(
            all_high_coverage_ht.possible_variants * all_high_coverage_ht.mu_snp
        )
    )

    # Generate a Table with all necessary annotations (x and y listed above) for the coverage model
    low_coverage_ht = all_low_coverage_ht.group_by(
        log_coverage=hl.log10(all_low_coverage_ht.exome_coverage)
    ).aggregate(
        low_coverage_obs_exp=hl.agg.sum(all_low_coverage_ht.observed_variants)
        / (
            high_coverage_scale_factor
            * hl.agg.sum(
                all_low_coverage_ht.possible_variants * all_low_coverage_ht.mu_snp
            )
        )
    )

    # Build the coverage model
    # TODO: consider weighting here as well
    coverage_model_expr = build_coverage_model(
        low_coverage_obs_exp=low_coverage_ht.low_coverage_obs_exp,
        log_coverage=low_coverage_ht.log_coverage,
    )
    coverage_model = tuple(low_coverage_ht.aggregate(coverage_model_expr).beta)
    return coverage_model, plateau_models


def build_plateau_models(
    cpg_expr: hl.expr.BooleanExpression,
    mu_snp_expr: hl.expr.Float64Expression,
    observed_variants_expr: hl.expr.Int64Expression,
    possible_variants_expr: hl.expr.Int64Expression,
    pop_observed_variants_exprs: List[hl.ArrayNumericExpression] = [],
    weighted: bool = False,
) -> Dict[str, Dict[bool, hl.expr.ArrayExpression]]:
    """
    Build plateau models for all the sites in `ht` and sites in each `pop` population downsampling to calibrate mutation rate to compute predicted proportion observed value.

    The x and y of the plateau models:
        - x: `mu_snp` - mutation rate
        - y: proportion observed ('observed_variants' or 'observed_{pop}' / 'possible_variants')

    :param cpg_expr: The BooleanExpression noting whether a site is a CPG site.
    :param mu_snp_expr: The Float64Expression of the mutation rate.
    :param observed_variants_expr: The Int64Expression of the observed variant counts for each combination of keys in `ht`.
    :param possible_variants_expr: The Int64Expression of the possible variant counts for each combination of keys in `ht`.
    :param pop_observed_variants_exprs: List of ArrayNumericExpression of observed variant counts for specified populations. Default is [].
    :param weighted: Whether to generalize the model to weighted least squares using 'possible_variants'. Default is False.
    :return: A Dictionary of intercepts and slopes for plateau models of each population.
    """
    # Build a plateau model using all the sites in the Table
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
    # Build plateau models using sites in population downsamplings if population is specified
    for pop, pop_observed_variants_expr in pop_observed_variants_exprs.items():
        plateau_models_agg_expr[pop] = hl.agg.array_agg(
            lambda pop_observed_variants: hl.agg.group_by(
                cpg_expr,
                hl.agg.linreg(
                    pop_observed_variants / possible_variants_expr,
                    [1, mu_snp_expr],
                    weight=possible_variants_expr if weighted else None,
                ).beta,
            ),
            pop_observed_variants_expr,
        )
    return plateau_models_agg_expr


def build_coverage_model(
    low_coverage_obs_exp_expr: hl.expr.Float64Expression,
    log_coverage_expr: hl.expr.Float64Expression,
) -> hl.expr.StructExpression:
    """
    Build coverage model.

    This function uses linear regression to build a model of log10(coverage) to correct proportion of expected variation at low coverage sites.

    The x and y of the coverage model:
        - x: log10('exome_coverage') at low coverage site
        - y: sum('observed_variants')/ (`high_coverage_scale_factor` * sum('possible_variants' * 'mu_snp') at low coverage site
        where `high_coverage_scale_factor` = sum('observed_variants') / sum('possible_variants' * 'mu_snp') at high coverage site

    :param low_coverage_obs_exp_expr: The Float64Expression of observed:expected ratio for a given coverage level.
    :param log_coverage_expr: The Float64Expression of log10 coverage.
    :return: Tuple with intercept and slope of the model.
    """
    return hl.agg.linreg(low_coverage_obs_exp_expr, [1, log_coverage_expr])


def get_all_pop_lengths(
    ht,
    pops: Tuple[str],
    prefix: str = "observed_",
) -> List[Tuple[str, str]]:
    """
    Get the minimum length of observed variant counts array for each population downsamping.

    The minimum array length will be used to determine the number of plateau models built for each population.

    The annotations are specified by the combination of `prefix` and each population in `pops`.

    :param ht: Input Table used to build population plateau models.
    :param pops: Populations used to categorize observed variant counts in downsampings.
    :param prefix: Prefix of population observed variant counts. Defaults to 'observed_'.
    :return: A Dictionary with the minimum array length for each population.
    """
    # Get minimum length of downsamplings for each population
    pop_downsampling_lengths = ht.aggregate(
        [hl.agg.min(hl.len(ht[f"{prefix}{pop}"])) for pop in pops]
    )

    # Zip population name with their downsampling length
    pop_lengths = list(zip(pop_downsampling_lengths, pops))
    logger.info("Found: %s", "".join(map(str, pop_lengths)))

    assert ht.all(
        hl.all(
            lambda f: f,
            [hl.len(ht[f"{prefix}{pop}"]) == length for length, pop in pop_lengths],
        )
    ), "The arrays of variant counts within population downsamplings have different lengths!"

    return pop_lengths
