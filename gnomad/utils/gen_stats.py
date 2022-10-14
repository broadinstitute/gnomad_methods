# noqa: D100

import logging

import hail as hl

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def to_phred(linear_expr: hl.expr.NumericExpression) -> hl.expr.Float64Expression:
    """
    Compute the phred-scaled value of the linear-scale input.

    :param linear_expr: input
    :return: Phred-scaled value
    """
    return -10 * hl.log10(linear_expr)


def from_phred(
    phred_score_expr: hl.expr.NumericExpression,
) -> hl.expr.Float64Expression:
    """
    Compute the linear-scale value of the phred-scaled input.

    :param phred_score_expr: phred-scaled value
    :return: linear-scale value of the phred-scaled input.
    """
    return 10 ** -(phred_score_expr / 10)


def get_median_and_mad_expr(
    metric_expr: hl.expr.ArrayNumericExpression, k: float = 1.4826
) -> hl.expr.StructExpression:
    """
    Compute the median and median absolute deviation (MAD) for the given expression.

    ..note::

        The default value of k assumes normally distributed data.

    :param metric_expr: Expression to compute median and MAD for
    :param k: The scaling factor for MAD calculation. Default assumes normally distributed data.
    :return: Struct with median and MAD
    """
    return hl.bind(
        lambda x: hl.struct(median=x[1], mad=k * hl.median(hl.abs(x[0] - x[1]))),
        hl.bind(lambda x: hl.tuple([x, hl.median(x)]), hl.agg.collect(metric_expr)),
    )


def merge_stats_counters_expr(
    stats: hl.expr.ArrayExpression,
) -> hl.expr.StructExpression:
    """
    Merge multiple stats counters, assuming that they were computed on non-overlapping data.

    Examples:
    - Merge stats computed on indel and snv separately
    - Merge stats computed on bi-allelic and multi-allelic variants separately
    - Merge stats computed on autosomes and sex chromosomes separately

    :param stats: An array of stats counters to merge
    :return: Merged stats Struct
    """

    def add_stats(
        i: hl.expr.StructExpression, j: hl.expr.StructExpression
    ) -> hl.expr.StructExpression:
        """
        Merge two stats counters together.

        Assumes that all stats counter fields are present in the struct.

        :param i: accumulator: struct with mean, n and variance
        :param j: new element: stats_struct -- needs to contain mean, n and variance
        :return: Accumulation over all elements: struct with mean, n and variance
        """
        delta = j.mean - i.mean
        n_tot = i.n + j.n
        return hl.struct(
            min=hl.min(i.min, j.min),
            max=hl.max(i.max, j.max),
            mean=(i.mean * i.n + j.mean * j.n) / n_tot,
            variance=i.variance + j.variance + (delta * delta * i.n * j.n) / n_tot,
            n=n_tot,
            sum=i.sum + j.sum,
        )

    # Gather all metrics present in all stats counters
    metrics = set(stats[0])
    dropped_metrics = set()
    for stat_expr in stats[1:]:
        stat_expr_metrics = set(stat_expr)
        dropped_metrics = dropped_metrics.union(stat_expr_metrics.difference(metrics))
        metrics = metrics.intersection(stat_expr_metrics)
    if dropped_metrics:
        logger.warning(
            "The following metrics will be dropped during stats counter merging as they"
            " do not appear in all counters: %s",
            ", ".join(dropped_metrics),
        )

    # Because merging standard deviation requires having the mean and n,
    # check that they are also present if `stdev` is. Otherwise remove stdev
    if "stdev" in metrics:
        missing_fields = [x for x in ["n", "mean"] if x not in metrics]
        if missing_fields:
            logger.warning(
                "Cannot merge `stdev` from given stats counters since they are missing"
                " the following fields: %s",
                ",".join(missing_fields),
            )
            metrics.remove("stdev")

    # Create a struct with all possible stats for merging.
    # This step helps when folding because we can rely on the struct schema
    # Note that for intermediate merging, we compute the variance rather than the stdev
    all_stats = hl.array(stats).map(
        lambda x: hl.struct(
            min=x.min if "min" in metrics else hl.null(hl.tfloat64),
            max=x.max if "max" in metrics else hl.null(hl.tfloat64),
            mean=x.mean if "mean" in metrics else hl.null(hl.tfloat64),
            variance=x.stdev * x.stdev if "stdev" in metrics else hl.null(hl.tfloat64),
            n=x.n if "n" in metrics else hl.null(hl.tfloat64),
            sum=x.sum if "sum" in metrics else hl.null(hl.tfloat64),
        )
    )

    # Merge the stats
    agg_stats = all_stats[1:].fold(add_stats, all_stats[0])

    # Return only the metrics that were present in all independent stats counters
    # If `stdev` is present, then compute it from the variance
    return agg_stats.select(
        **{
            metric: agg_stats[metric]
            if metric != "stdev"
            else hl.sqrt(agg_stats.variance)
            for metric in metrics
        }
    )
