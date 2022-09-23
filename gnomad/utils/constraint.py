# noqa: D100

import logging
from typing import List, Optional, Tuple, Union

import hail as hl

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_utils")
logger.setLevel(logging.INFO)

COVERAGE_CUTOFF = 40


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


def build_models(
    coverage_ht: hl.Table,
    weighted: bool = False,
    pops: Optional[Tuple[str]] = (),
    keys: Tuple[str] = (
        "context",
        "ref",
        "alt",
        "methylation_level",
        "mu_snp",
    ),
    cov_cutoff: int = COVERAGE_CUTOFF,
) -> Tuple[Tuple[float, float], hl.Struct]:
    """
    Build coverage and plateau models.

    Coverage model: corrects proportion of expected variation at low coverage sites. Low coverage sites are defined as sites with median coverage < `HIGH_COVERAGE_CUTOFF`.
    The x and y of the plateau models:
    - x: `mu_snp` - mutation rate
    - y: proportion observed ('observed_variants' or 'observed_{pop}' / 'possible_variants')

    This function builds plateau models to calibrate mutation rate estimates against the proportion observed
    of each substitution, context, and methylation level in `coverage_ht` considering only high coverage sites,
    or sites above a median coverage of `HIGH_COVERAGE_CUTOFF`. Since the output of `get_proportion_observed_by_coverage`
    is used as `coverage_ht`, the proportion observed will be high-quality variants below 0.1% frequency at synonymous
    sites. Two plateau models are fit, one for CpG transitions and one for the remainder of sites (transversions and non-CpG transitions).

    Plateau model: adjusts proportion of expected variation based on location in the genome and CpG status.
    The x and y of the coverage model:
    - x: log10('exome_coverage') at low coverage site
    - y: sum('observed_variants')/ (`high_coverage_scale_factor` * sum('possible_variants' * 'mu_snp') at low coverage site
    where `high_coverage_scale_factor` = sum('observed_variants') / sum('possible_variants' * 'mu_snp') at high coverage site

    For low coverage sites, or sites below `HIGH_COVERAGE_CUTOFF`, this function performs a base-level resolution rather than exon-level to compute a coverage correction factor
    to reduce the inaccuracy of expected variant counts caused by low coverage on each base. The coverage models are built
    by first defining a metric that is derived by dividing the number of observed variants with the total number of
    possible variants times the mutation rate summed across all substitutions, contexts, and methylation level. Since
    the output of `get_proportion_observed_by_coverage` is used as `coverage_ht`, the number of observed variants and possible
    variants will be at synonymous sites. The function computes this metric for high coverage sites as a global scaling
    factor, and divides this metric at low coverage sites by this scaling factor to create an observed:expected ratio. Then the
    coverage model is built of log10(coverage) to this scaled ratio as a correction factor for low coverage sites.

    .. note::
        This function expects that the input `coverage_ht` is the output of `get_proportion_observed_by_coverage`, and
        therefore the following fields should be present in `coverage_ht`:
        - context - trinucleotide genomic context
        - ref - the reference allele
        - alt - the alternate base
        - methylation_level - methylation level
        - exome_coverage - median exome coverage at integer values between 1-100
        - observed_variants - the number of observed variants in the dataset for each variant. Note that the term "variant" here refers to a specific substitution, context, methylation level, and coverage combination.
        - downsampling_counts_{pop} (pop defaults to ()) - array of observed variant counts per population after downsampling
        - mu_snp - mutation rate
        - possible_variants - the number of possible variants in the dataset for each variant. Note that the term "variant" here refers to a specific substitution, context, methylation level, and coverage combination.

    :param coverage_ht: Input coverage Table.
    :param weighted: Whether to weight the high coverage model (a linear regression model) by 'possible_variants'. Defaults to False.
    :param pops: List of populations used to build coverage and plateau models. Defaults to ().
    :param keys: Keys to group observed and possible variant counts.
    :param cov_cutoff: Median coverage cutoff. Sites with coverage above this cutoff are considered well covered and will be used to build plateau models. Sites below this cutoff have low coverage and will be used to build coverage models. Defaults to `HIGH_COVERAGE_CUTOFF`.
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
    high_coverage_ht = all_high_coverage_ht.group_by(*keys).aggregate(**agg_expr)
    high_coverage_ht = annotate_mutation_type(high_coverage_ht)

    # Build plateau models
    plateau_models = build_plateau_models_total(
        high_coverage_ht,
        cpg="cpg",
        mu_snp="mu_snp",
        observed_variants="observed_variants",
        possible_variants="possible_variants",
        weighted=weighted,
    )
    # Build plateau models for populations if `pops` is given
    if pops:
        pop_lengths = get_all_pop_lengths(high_coverage_ht, pops)
        for length, pop in pop_lengths:
            pop_plateau_model = build_plateau_models_pop(
                high_coverage_ht,
                cpg="cpg",
                mu_snp="mu_snp",
                observed_variants=f"observed_{pop}",
                possible_variants="possible_variants",
                pop=pop,
                length=length,
                weighted=weighted,
            )
            plateau_models = plateau_models.annotate(**pop_plateau_model[pop])

    # Build coverage model
    all_low_coverage_ht = coverage_ht.filter(
        (coverage_ht.exome_coverage < cov_cutoff) & (coverage_ht.exome_coverage > 0)
    )
    # a metric that represents the relative mutability of the exome for high coverage sites as a global scaling factor
    high_coverage_scale_factor = all_high_coverage_ht.aggregate(
        hl.agg.sum(all_high_coverage_ht.observed_variants)
        / hl.agg.sum(
            all_high_coverage_ht.possible_variants * all_high_coverage_ht.mu_snp
        )
    )
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
    coverage_model = build_coverage_model(
        low_coverage_ht,
        low_coverage_obs_exp="low_coverage_obs_exp",
        log_coverage="log_coverage",
    )
    # TODO: consider weighting here as well

    return coverage_model, plateau_models


def build_plateau_models_pop(
    ht: hl.Table,
    cpg: str,
    mu_snp: str,
    observed_variants: str,
    possible_variants: str,
    pop: str,
    length: int,
    weighted: bool = False,
) -> hl.Struct:
    """
    Build plateau model for each `pop` population to calibrate mutation rate to compute predicted proportion observed value.

    The x and y of the plateau models:
    - x: `mu_snp` - mutation rate
    - y: proportion observed (observed_{pop}' / 'possible_variants')

    :param ht: High coverage Table.
    :param cpg: The annotation name of booleans that determine whether the site is CPG site or not.
    :param mu_snp: The annotation name of the mutation rate.
    :param observed_variants: The annotation name of the observed variant counts for each combination of keys in `ht`.
    :param possible_variants: The annotation name of the possible variant counts for each combination of keys in `ht`.
    :param pop: List of populations. Defaults to None.
    :param length: The minimum array length of variant counts array in downsampings for `pop`.
    :param weighted: Whether to generalize the model to weighted least squares using 'possible_variants', defaults to False.
    :return: A Dictionary of intercepts and slopes for plateau models of each population.
    """
    agg_expr = {
        pop: [
            hl.agg.group_by(
                ht[cpg],
                hl.agg.linreg(
                    ht[observed_variants][i] / ht[possible_variants],
                    [1, ht[mu_snp]],
                    weight=ht[possible_variants] if weighted else None,
                ).beta,
            )
            for i in range(length)
        ]
    }
    return ht.aggregate(hl.struct(**agg_expr))


def build_plateau_models_total(
    ht: hl.Table,
    cpg: str,
    mu_snp: str,
    observed_variants: str,
    possible_variants: str,
    weighted: bool = False,
) -> hl.Struct:
    """
    Build plateau model for all release samples without downsamplings to calibrate mutation rate to compute predicted proportion observed value.

    The x and y of the plateau models:
    - x: `mu_snp` - mutation rate
    - y: proportion observed ('observed_variants'/ 'possible_variants')

    :param ht: High coverage Table.
    :param cpg: The annotation name of booleans that determine whether the site is CPG site or not.
    :param mu_snp: The annotation name of the mutation rate.
    :param observed_variants: The annotation name of the observed variant counts for each combination of keys in `ht`.
    :param possible_variants: The annotation name of the possible variant counts for each combination of keys in `ht`.
    :param weighted: Whether to generalize the model to weighted least squares using 'possible_variants', defaults to False.
    :return: A Dictionary of intercepts and slopes for a plateau model.
    """
    agg_expr = {
        "total": hl.agg.group_by(
            ht[cpg],
            hl.agg.linreg(
                ht[observed_variants] / ht[possible_variants],
                [1, ht[mu_snp]],
                weight=ht[possible_variants] if weighted else None,
            ).beta,
        )
    }
    return ht.aggregate(hl.struct(**agg_expr))


def build_coverage_model(
    ht: hl.Table,
    low_coverage_obs_exp: str,
    log_coverage: str,
) -> Tuple[float, float]:
    """
    Build coverage model.

    This function uses linear regression to build a model of log10(coverage) to this scaled ratio as a correction
    factor for predicted proportion observed/expected variant counts at low coverage sites.

    The x and y of the coverage model:
    - x: log10('exome_coverage') at low coverage site
    - y: sum('observed_variants')/ (`high_coverage_scale_factor` * sum('possible_variants' * 'mu_snp') at low coverage site
    where `high_coverage_scale_factor` = sum('observed_variants') / sum('possible_variants' * 'mu_snp') at high coverage site

    :param ht: Low coverage Table.
    :param low_cov_oe_expr: The annotation name of observed:expected ratio for a given coverage level.
    :param log_cov_expr: The annotation name of log10 coverage.
    :return: Tuple with intercept and slope of the model.
    """
    return tuple(
        ht.aggregate(
            hl.agg.linreg(ht[low_coverage_obs_exp], [1, ht[log_coverage]])
        ).beta
    )


def get_all_pop_lengths(
    ht, pops: Tuple[str], prefix: str = "observed_"
) -> List[Tuple[str, str]]:
    """
    Get the minimum array length of variant counts array in downsampings for specific per population annotations in `ht`.

    The minimum array length will be used to determine the number of plateau models built for each population.

    The annotations are specified by the combination of `prefix` and each population in `pops`.

    :param ht: Input Table used to build population plateau models.
    :param pops: List of populations.
    :param prefix: Prefix of population variant count. Defaults to 'observed_'.
    :return: A Dictionary with the minimum array length for each population.
    """
    pop_downsampling_lengths = ht.aggregate(
        [hl.agg.min(hl.len(ht[f"{prefix}{pop}"])) for pop in pops]
    )

    pop_lengths = list(zip(pop_downsampling_lengths, pops))
    logger.info("Found: %s", "".join(map(str, pop_lengths)))

    # check all the arrays of variant counts within a population have the same length
    assert ht.all(
        hl.all(
            lambda f: f,
            [hl.len(ht[f"{prefix}{pop}"]) == length for length, pop in pop_lengths],
        )
    )
    return pop_lengths
