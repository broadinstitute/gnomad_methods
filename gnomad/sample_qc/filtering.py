import logging
from typing import Dict, Iterable, List, Optional, Tuple

import hail as hl
from gnomad.utils.gen_stats import get_median_and_mad_expr, merge_stats_counters_expr
from hail.utils.misc import divide_null

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def compute_qc_metrics_residuals(
    ht: hl.Table,
    pc_scores: hl.expr.ArrayNumericExpression,
    qc_metrics: Dict[str, hl.expr.NumericExpression],
    use_pc_square: bool = True,
    n_pcs: Optional[int] = None,
    regression_sample_inclusion_expr: hl.expr.BooleanExpression = hl.bool(True),
) -> hl.Table:
    """
    Computes QC metrics residuals after regressing out PCs (and optionally PC^2)

    .. note::

        The `regression_sample_inclusion_expr` can be used to select a subset of the samples to include in the regression calculation.
        Residuals are always computed for all samples.

    :param ht: Input sample QC metrics HT
    :param pc_scores: The expression in the input HT that stores the PC scores
    :param qc_metrics: A dictionary with the name of each QC metric to compute residuals for and their corresponding expression in the input HT.
    :param use_pc_square: Whether to  use PC^2 in the regression or not
    :param n_pcs: Numer of PCs to use. If not set, then all PCs in `pc_scores` are used.
    :param regression_sample_inclusion_expr: An optional expression to select samples to include in the regression calculation.
    :return: Table with QC metrics residuals
    """

    # Annotate QC HT with fields necessary for computation
    _sample_qc_ht = ht.select(
        **qc_metrics, scores=pc_scores, _keep=regression_sample_inclusion_expr
    )

    # If n_pcs wasn't provided, use all PCs
    if n_pcs is None:
        n_pcs = _sample_qc_ht.aggregate(hl.agg.min(hl.len(_sample_qc_ht.scores)))

    logger.info(
        "Computing regressed QC metrics filters using {} PCs for metrics: {}".format(
            n_pcs, ", ".join(qc_metrics)
        )
    )

    # Prepare regression variables, adding 1.0 first for the intercept
    # Adds square of variables if use_pc_square is true
    x_expr = [1.0] + [_sample_qc_ht.scores[i] for i in range(0, n_pcs)]
    if use_pc_square:
        x_expr.extend(
            [_sample_qc_ht.scores[i] * _sample_qc_ht.scores[i] for i in range(0, n_pcs)]
        )

    # Compute linear regressions
    lms = _sample_qc_ht.aggregate(
        hl.struct(
            **{
                metric: hl.agg.filter(
                    _sample_qc_ht._keep,
                    hl.agg.linreg(y=_sample_qc_ht[metric], x=x_expr),
                )
                for metric in qc_metrics
            }
        )
    )

    _sample_qc_ht = _sample_qc_ht.annotate_globals(lms=lms).persist()

    # Compute residuals
    def get_lm_prediction_expr(metric: str):
        lm_pred_expr = _sample_qc_ht.lms[metric].beta[0] + hl.sum(
            hl.range(n_pcs).map(
                lambda i: _sample_qc_ht.lms[metric].beta[i + 1]
                * _sample_qc_ht.scores[i]
            )
        )
        if use_pc_square:
            lm_pred_expr = lm_pred_expr + hl.sum(
                hl.range(n_pcs).map(
                    lambda i: _sample_qc_ht.lms[metric].beta[i + n_pcs + 1]
                    * _sample_qc_ht.scores[i]
                    * _sample_qc_ht.scores[i]
                )
            )
        return lm_pred_expr

    residuals_ht = _sample_qc_ht.select(
        **{
            f"{metric}_residual": _sample_qc_ht[metric] - get_lm_prediction_expr(metric)
            for metric in _sample_qc_ht.lms
        }
    )

    return residuals_ht.persist()


def compute_stratified_metrics_filter(
    ht: hl.Table,
    qc_metrics: Dict[str, hl.expr.NumericExpression],
    strata: Optional[Dict[str, hl.expr.Expression]] = None,
    lower_threshold: float = 4.0,
    upper_threshold: float = 4.0,
    metric_threshold: Optional[Dict[str, Tuple[float, float]]] = None,
    filter_name: str = "qc_metrics_filters",
) -> hl.Table:
    """
    Compute median, MAD, and upper and lower thresholds for each metric used in outlier filtering

    :param ht: HT containing relevant sample QC metric annotations
    :param qc_metrics: list of metrics (name and expr) for which to compute the critical values for filtering outliers
    :param strata: List of annotations used for stratification. These metrics should be discrete types!
    :param lower_threshold: Lower MAD threshold
    :param upper_threshold: Upper MAD threshold
    :param metric_threshold: Can be used to specify different (lower, upper) thresholds for one or more metrics
    :param filter_name: Name of resulting filters annotation
    :return: Table grouped by strata, with upper and lower threshold values computed for each sample QC metric
    """

    _metric_threshold = {
        metric: (lower_threshold, upper_threshold) for metric in qc_metrics
    }
    if metric_threshold is not None:
        _metric_threshold.update(metric_threshold)

    def make_filters_expr(
        ht: hl.Table, qc_metrics: Iterable[str]
    ) -> hl.expr.SetExpression:
        return hl.set(
            hl.filter(
                lambda x: hl.is_defined(x),
                [hl.or_missing(ht[f"fail_{metric}"], metric) for metric in qc_metrics],
            )
        )

    if strata is None:
        strata = {}

    ht = ht.select(**qc_metrics, **strata).key_by("s").persist()

    agg_expr = hl.struct(
        **{
            metric: hl.bind(
                lambda x: x.annotate(
                    lower=x.median - _metric_threshold[metric][0] * x.mad,
                    upper=x.median + _metric_threshold[metric][1] * x.mad,
                ),
                get_median_and_mad_expr(ht[metric]),
            )
            for metric in qc_metrics
        }
    )

    if strata:
        ht = ht.annotate_globals(
            qc_metrics_stats=ht.aggregate(
                hl.agg.group_by(hl.tuple([ht[x] for x in strata]), agg_expr),
                _localize=False,
            )
        )
        metrics_stats_expr = ht.qc_metrics_stats[hl.tuple([ht[x] for x in strata])]
    else:
        ht = ht.annotate_globals(
            qc_metrics_stats=ht.aggregate(agg_expr, _localize=False)
        )
        metrics_stats_expr = ht.qc_metrics_stats

    fail_exprs = {
        f"fail_{metric}": (ht[metric] <= metrics_stats_expr[metric].lower)
        | (ht[metric] >= metrics_stats_expr[metric].upper)
        for metric in qc_metrics
    }
    ht = ht.transmute(**fail_exprs)
    stratified_filters = make_filters_expr(ht, qc_metrics)
    return ht.annotate(**{filter_name: stratified_filters})


def compute_stratified_sample_qc(
    mt: hl.MatrixTable,
    strata: Dict[str, hl.expr.BooleanExpression],
    tmp_ht_prefix: Optional[str],
    gt_col: Optional[str] = None,
) -> hl.Table:
    """
    Runs hl.sample_qc on different strata and then also merge the results into a single expression.
    Note that strata should be non-overlapping, e.g. SNV vs indels or bi-allelic vs multi-allelic

    :param mt: Input MT
    :param strata: Strata names and filtering expressions
    :param tmp_ht_prefix: Optional path prefix to write the intermediate strata results to (recommended for larger datasets)
    :param gt_col: Name of entry field storing the genotype. Default: 'GT'
    :return: Sample QC table, including strat-specific numbers
    """
    mt = mt.select_rows(**strata)

    if gt_col is not None:
        mt = mt.select_entries(GT=mt[gt_col])
    else:
        mt = mt.select_entries("GT")

    strat_hts = {}
    for strat in strata:
        strat_sample_qc_ht = hl.sample_qc(mt.filter_rows(mt[strat])).cols()
        if tmp_ht_prefix is not None:
            strat_sample_qc_ht = strat_sample_qc_ht.checkpoint(
                tmp_ht_prefix + f"_{strat}.ht", overwrite=True
            )
        else:
            strat_sample_qc_ht = strat_sample_qc_ht.persist()
        strat_hts[strat] = strat_sample_qc_ht

    sample_qc_ht = strat_hts.pop(list(strata)[0])
    sample_qc_ht = sample_qc_ht.select(
        **{f"{list(strata)[0]}_sample_qc": sample_qc_ht.sample_qc},
        **{
            f"{strat}_sample_qc": strat_hts[strat][sample_qc_ht.key].sample_qc
            for strat in list(strata)[1:]
        },
    )
    sample_qc_ht = sample_qc_ht.annotate(
        sample_qc=merge_sample_qc_expr(list(sample_qc_ht.row_value.values()))
    )

    return sample_qc_ht


def merge_sample_qc_expr(
    sample_qc_exprs: List[hl.expr.StructExpression],
) -> hl.expr.StructExpression:
    """
    Creates an expression that merges results from non-overlapping strata of hail.sample_qc

    E.g.:

    - Compute autosomes and sex chromosomes metrics separately, then merge results
    - Compute bi-allelic and multi-allelic metrics separately, then merge results

    Note regarding the merging of ``dp_stats`` and ``gq_stats``:
    Because ``n`` is needed to aggregate ``stdev``, ``n_called`` is used for this purpose.
    This should work very well on a standard GATK VCF and it essentially assumes that:

    - samples that are called have `DP` and `GQ` fields
    - samples that are not called do not have `DP` and `GQ` fields

    Even if these assumptions are broken for some genotypes, it shouldn't matter too much.

    :param sample_qc_exprs: List of sample QC struct expressions for each stratification
    :return: Combined sample QC results
    """

    # List of metrics that can be aggregated by summing
    additive_metrics = [
        "n_called",
        "n_not_called",
        "n_filtered",
        "n_hom_ref",
        "n_het",
        "n_hom_var",
        "n_snp",
        "n_insertion",
        "n_deletion",
        "n_singleton",
        "n_transition",
        "n_transversion",
        "n_star",
    ]

    # List of metrics that are ratio of summed metrics (name, nominator, denominator)
    ratio_metrics = [
        ("call_rate", "n_called", "n_not_called"),
        ("r_ti_tv", "n_transition", "n_transversion"),
        ("r_het_hom_var", "n_het", "n_hom_var"),
        ("r_insertion_deletion", "n_insertion", "n_deletion"),
    ]

    # List of metrics that are struct generated by a stats counter
    stats_metrics = ["gq_stats", "dp_stats"]

    # Gather metrics present in sample qc fields
    sample_qc_fields = set(sample_qc_exprs[0])
    for sample_qc_expr in sample_qc_exprs[1:]:
        sample_qc_fields = sample_qc_fields.union(set(sample_qc_expr))

    # Merge additive metrics in sample qc fields
    merged_exprs = {
        metric: hl.sum([sample_qc_expr[metric] for sample_qc_expr in sample_qc_exprs])
        for metric in additive_metrics
        if metric in sample_qc_fields
    }

    # Merge ratio metrics in sample qc fields
    merged_exprs.update(
        {
            metric: hl.float64(divide_null(merged_exprs[nom], merged_exprs[denom]))
            for metric, nom, denom in ratio_metrics
            if nom in sample_qc_fields and denom in sample_qc_fields
        }
    )

    # Merge stats counter metrics in sample qc fields
    # Use n_called as n for DP and GQ stats
    if "n_called" in sample_qc_fields:
        merged_exprs.update(
            {
                metric: merge_stats_counters_expr(
                    [
                        sample_qc_expr[metric].annotate(n=sample_qc_expr.n_called)
                        for sample_qc_expr in sample_qc_exprs
                    ]
                ).drop("n")
                for metric in stats_metrics
                if metric in sample_qc_fields
            }
        )

    return hl.struct(**merged_exprs)
