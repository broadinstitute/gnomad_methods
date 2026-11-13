"""
Functions to compute constraint metrics.

Covers pLI, the o/e confidence interval, z-scores, flags, and GENCODE annotations.
"""

import copy
from typing import Dict, List, Optional, Tuple, Union

import hail as hl
from hail.utils.misc import divide_null, new_temp_file

from gnomad.utils.constraint.constants import DEFAULT_GENCODE_ANNOTATIONS


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


def _oe_ci_discretized_poisson(
    obs_expr: hl.expr.Int64Expression,
    exp_expr: hl.expr.Float64Expression,
    alpha: float = 0.05,
) -> hl.expr.StructExpression:
    """
    Compute OE confidence interval via discretized Poisson CDF.

    Sweeps the OE ratio parameter over [0, 2) in steps of 0.001, evaluates
    ``dpois(obs, exp * x)`` at each point, normalises the cumulative sum, and
    reads off the bounds at ``alpha`` and ``1 - alpha``.

    :param obs_expr: Observed variant count expression.
    :param exp_expr: Expected variant count expression.
    :param alpha: Significance level. Default is 0.05.
    :return: Struct with ``lower`` and ``upper`` bounds.
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


def _oe_ci_gamma(
    obs_expr: hl.expr.Int64Expression,
    exp_expr: hl.expr.Float64Expression,
    alpha: float = 0.05,
) -> hl.expr.StructExpression:
    """
    Compute OE confidence interval using the Gamma distribution.

    Uses Hail's ``hl.qgamma`` quantile function.

    :param obs_expr: Observed variant count expression.
    :param exp_expr: Expected variant count expression.
    :param alpha: Significance level. Default is 0.05.
    :return: Struct with ``lower`` and ``upper`` bounds.
    """
    shape = obs_expr + hl.literal(1.0)
    scale = divide_null(hl.literal(1.0), exp_expr)
    return hl.struct(
        lower=hl.qgamma(hl.literal(alpha), shape, scale),
        upper=hl.qgamma(hl.literal(1.0 - alpha), shape, scale),
    )


def oe_confidence_interval(
    obs_expr: hl.expr.Int64Expression,
    exp_expr: hl.expr.Float64Expression,
    alpha: float = 0.05,
    method: str = "gamma",
) -> hl.expr.StructExpression:
    """
    Compute a confidence interval around the observed/expected ratio.

    Two methods are available:

    - ``"gamma"`` (default): uses ``hl.qgamma`` to compute exact quantiles of
      the Gamma posterior. Fast and precise.
    - ``"poisson"``: sweeps a discretized Poisson likelihood over the OE
      parameter space [0, 2). Retained for backwards compatibility.

    .. note::

        The default changed from the discretized Poisson computation to
        ``"gamma"``, so calls that do not pass ``method`` return different bounds
        than before. ``"gamma"`` requires Hail >= 0.2.137, which added
        ``hl.qgamma``; on older Hail versions pass ``method="poisson"``.

    :param obs_expr: Observed variant count expression.
    :param exp_expr: Expected variant count expression.
    :param alpha: Significance level for the confidence interval. Default is
        0.05 (90% CI).
    :param method: CI method — ``"gamma"`` or ``"poisson"``. Default is
        ``"gamma"``.
    :return: Struct with ``lower`` and ``upper`` bounds.
    """
    if method == "gamma":
        if not hasattr(hl, "qgamma"):
            raise RuntimeError(
                "method='gamma' requires hl.qgamma, available in Hail >= 0.2.137. "
                "Use method='poisson' or upgrade Hail."
            )
        return _oe_ci_gamma(obs_expr, exp_expr, alpha)
    elif method == "poisson":
        return _oe_ci_discretized_poisson(obs_expr, exp_expr, alpha)
    else:
        raise ValueError(f"Unknown CI method: {method!r}. Use 'gamma' or 'poisson'.")


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
    annotations: Union[List[str], Tuple[str, ...]] = DEFAULT_GENCODE_ANNOTATIONS,
    remove_y_par: bool = True,
) -> hl.Table:
    """
    Add GENCODE annotations to Table based on transcript id.

    In addition to the annotations specified by ``annotations``, the following
    computed annotations are always added:

        - chromosome
        - cds_length
        - num_coding_exons

    :param ht: Input Table.
    :param gencode_ht: Table with GENCODE annotations.
    :param annotations: GENCODE annotations to add. Default is
        ``DEFAULT_GENCODE_ANNOTATIONS``.
    :param remove_y_par: Whether to remove features for the Y chromosome PAR regions.
        Default is True because the Y chromosome PAR regions are typically not included
        in the constraint calculations and both chrX and chrY will have the same
        'transcript_id' field for these regions. This parameter can only be True if
        ``gencode_ht`` includes a 'transcript_id_version' field because Y_PAR (or
        PAR_Y, depending on the GENCODE release) is included in the version of the
        transcript, which has been stripped from the 'transcript_id' field.
    :return: Table with transcript annotations from GENCODE added.
    """
    annotations = list(annotations)

    if remove_y_par and "transcript_id_version" not in gencode_ht.row:
        raise ValueError(
            "remove_y_par is True but 'transcript_id_version' is not in gencode_ht"
        )

    if remove_y_par:
        gencode_ht = gencode_ht.filter(
            ~gencode_ht.transcript_id_version.endswith("Y_PAR")
            & ~gencode_ht.transcript_id_version.endswith("PAR_Y")
        )

    gencode_ht = gencode_ht.annotate(
        length=gencode_ht.interval.end.position
        - gencode_ht.interval.start.position
        + 1,
        chromosome=gencode_ht.interval.start.contig,
        start_position=gencode_ht.interval.start.position,
        end_position=gencode_ht.interval.end.position,
    )

    # Add transcript annotations to input Table.
    annotations_to_add = set(annotations + ["chromosome", "transcript_id"])
    gencode_transcript_ht = (
        gencode_ht.filter(gencode_ht.feature == "transcript")
        .select(*annotations_to_add)
        .key_by("transcript_id")
        .drop("interval")
    )

    # Obtain CDS annotations from GENCODE file and calculate CDS length and
    # number of exons.
    annotations_to_add = set(annotations + ["chromosome", "transcript_id", "length"])

    gencode_cds_ht = (
        gencode_ht.filter(gencode_ht.feature == "CDS")
        .select(*annotations_to_add)
        .key_by("transcript_id")
        .drop("interval")
    )

    annotations_to_add.remove("length")

    gencode_cds_ht = gencode_cds_ht.group_by("transcript_id").aggregate(
        cds_length=hl.agg.sum(gencode_cds_ht.length),
        num_coding_exons=hl.agg.count(),
    )

    # Join the transcript and CDS annotations.
    gencode_transcript_ht = gencode_transcript_ht.join(gencode_cds_ht).checkpoint(
        new_temp_file(prefix="gencode", extension="ht")
    )

    # Add GENCODE annotations to input Table.
    ht = ht.annotate(**gencode_transcript_ht[ht.transcript])

    return ht
