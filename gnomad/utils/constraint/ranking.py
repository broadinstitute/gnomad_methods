"""Functions to rank transcripts by a metric and assign percentile bins."""

import logging
from typing import Callable, Dict, List, Optional, Tuple, Union

import hail as hl

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_utils")
logger.setLevel(logging.INFO)


def rank_and_assign_bins(
    value_expr: hl.expr.Float64Expression,
    bin_granularities: Optional[Dict[str, int]] = None,
    prefix: str = "",
) -> hl.StructExpression:
    """
    Rank rows by a numeric expression and assign bin labels.

    **Rank-based binning**: every row receives a unique position in the sorted
    order, and bins are derived from that position. This differs from
    threshold-based binning (see :func:`annotate_bins_by_threshold`), where
    bins are assigned by comparing values against pre-computed boundary values.

    Rows are ordered ascending by ``value_expr``. Each row is assigned a
    0-based ``{prefix}rank`` and a ``{prefix}bin_{name}`` field for every
    entry in ``bin_granularities``, computed as
    ``hl.int(rank * multiplier / n_rows)``.

    Rows where ``value_expr`` is missing are excluded from the ranking and
    receive missing rank and bin annotations. ``n_rows`` is therefore the
    number of rows with a defined ``value_expr``, so bins are evenly sized
    and the full range of bins is populated.

    Used by :func:`rank_array_element_metrics` to rank metrics within array
    elements.

    .. note::

        This function evaluates its input when it is called. It needs the number
        of rows with a defined ``value_expr`` for the bin denominator, so it counts
        them, and it caches the filtered and the ranked rows so the filter and the
        sort are not recomputed. In Hail, ``Table.cache()`` ignores the storage
        level and writes a checkpoint to a ``persist_Table*`` path under the
        temporary directory set by ``hl.init(tmp_dir=...)``, so with a GCS
        ``tmp_dir`` these Tables are written to GCS. They hold only ``value_expr``
        and the key, and are not removed when the function returns.

    :param value_expr: Numeric expression to rank by (ascending).
    :param bin_granularities: Mapping of bin name to multiplier. Each entry
        produces a ``{prefix}bin_{name}`` field. Default is
        ``{"percentile": 100, "decile": 10, "sextile": 6}``.
    :param prefix: String prepended to ``rank`` and ``bin_{name}`` field
        names. Default is ``""`` (no prefix).
    :return: Struct with ``{prefix}rank`` and ``{prefix}bin_{name}`` fields
        for each entry in ``bin_granularities``.
    """
    if bin_granularities is None:
        bin_granularities = {"percentile": 100, "decile": 10, "sextile": 6}

    ht = value_expr._indices.source
    source_key = list(ht.key)
    ranked_ht = ht.select(_=value_expr)

    # Rank only rows with a defined value. `order_by` sorts missing values after
    # non-missing ones, so retaining them would assign them the highest ranks and
    # inflate the denominator, shifting every row into a lower bin and leaving the
    # top bins empty. Rows dropped here get missing rank/bin fields from the join.
    # Cache before counting so the filter is not recomputed when ranking below.
    ranked_ht = ranked_ht.filter(hl.is_defined(ranked_ht._)).cache()
    n_rows = ranked_ht.count()

    ranked_ht = ranked_ht.order_by("_").add_index("rank")
    ranked_ht = ranked_ht.select(
        *source_key,
        **{f"{prefix}rank": ranked_ht.rank},
        **{
            f"{prefix}bin_{name}": hl.int(ranked_ht.rank * multiplier / n_rows)
            for name, multiplier in bin_granularities.items()
        },
    ).cache()

    return ranked_ht.key_by(*source_key)[ht.key]


def compute_percentile_thresholds(
    ht: hl.Table,
    metric_expr: hl.expr.Float64Expression,
    outlier_expr: Optional[hl.expr.BooleanExpression] = None,
    transcript_filter_expr: Optional[hl.expr.BooleanExpression] = None,
    percentiles: Tuple[float, ...] = (1, 5, 10, 15, 25, 50, 75),
    quantile_k: int = 1000,
) -> Dict[float, float]:
    """
    Compute approximate percentile thresholds for a metric expression.

    **Threshold-based binning, step 1**: computes the boundary values that
    define bin edges. Its values, keyed by ``(granularity, metric)`` instead of
    percentile, form the ``thresholds`` argument of
    :func:`annotate_bins_by_threshold` (step 2), which assigns each row to a bin;
    see that function for an example.

    This two-step approach differs from rank-based binning (see
    :func:`rank_and_assign_bins`), where every row receives a unique position
    in the sorted order. Threshold-based binning allows thresholds to be
    computed on a filtered subset (e.g., representative transcripts) and then
    applied to all rows, so multiple rows can share the same bin.

    Optionally filters to a subset of rows and excludes outliers, then
    computes approximate quantile thresholds at the requested percentiles in
    a single aggregation pass.

    .. note::

        Uses ``hl.agg.approx_quantiles``, so results are approximate. Increase
        ``quantile_k`` for higher accuracy.

    :param ht: Input Table.
    :param metric_expr: Float expression to compute thresholds for. Must be
        defined on ``ht``.
    :param outlier_expr: Optional boolean expression that is ``True`` for rows
        to exclude. When None (default), no outlier filtering is applied.
    :param transcript_filter_expr: Optional boolean expression that is ``True``
        for rows to include. When None (default), all rows are included.
    :param percentiles: Percentile values (0-100) at which to compute
        thresholds. Default is (1, 5, 10, 15, 25, 50, 75).
    :param quantile_k: Accuracy parameter for
        :func:`hail.expr.aggregators.approx_quantiles`. Default is 1000.
    :return: Dict mapping each percentile to its threshold value.
    """
    qs = [p / 100.0 for p in percentiles]
    filter_expr = hl.is_defined(metric_expr)
    if outlier_expr is not None:
        filter_expr = filter_expr & ~outlier_expr
    if transcript_filter_expr is not None:
        filter_expr = filter_expr & transcript_filter_expr

    result = ht.aggregate(
        hl.struct(
            thresholds=hl.agg.filter(
                filter_expr, hl.agg.approx_quantiles(metric_expr, qs, k=quantile_k)
            ),
            n=hl.agg.count_where(filter_expr),
        )
    )
    logger.info(
        "Computed percentile thresholds on %d transcripts.",
        result.n,
    )

    return dict(zip(percentiles, result.thresholds))


def annotate_bins_by_threshold(
    ht: hl.Table,
    metric_exprs: Dict[str, hl.expr.Float64Expression],
    thresholds: Dict[Tuple[str, str], List[float]],
    granularities: Union[List[str], Tuple[str, ...]],
    field_name: str = "constraint_bins",
) -> hl.Table:
    """
    Annotate rows with bin assignments using pre-computed thresholds.

    **Threshold-based binning, step 2**: assigns each row to a bin by
    comparing its metric value against boundary values produced by
    :func:`compute_percentile_thresholds` (step 1). For each
    ``(granularity, metric)`` pair, the bin equals the number of threshold
    boundaries the value exceeds. Bin 0 is the most constrained (below all
    thresholds); bin N means the value exceeds all N boundaries.

    This differs from rank-based binning (see :func:`rank_and_assign_bins`),
    where every row gets a unique position. Here, multiple rows can share a
    bin, and the thresholds may have been derived from a different subset of
    rows than those being annotated.

    :func:`compute_percentile_thresholds` returns a dict of percentile to value,
    so ``thresholds`` is built from it by selecting the percentiles for each
    granularity. Example use, with decile bins for two metrics:

    .. code-block:: python

        deciles = [10, 20, 30, 40, 50, 60, 70, 80, 90]
        metric_exprs = {"lof": ht.lof_oe_upper, "mis": ht.mis_oe_upper}
        thresholds = {}
        for metric, expr in metric_exprs.items():
            pct_to_value = compute_percentile_thresholds(ht, expr, percentiles=deciles)
            thresholds[("decile", metric)] = [pct_to_value[p] for p in deciles]

        ht = annotate_bins_by_threshold(ht, metric_exprs, thresholds, ["decile"])

    :param ht: Input Table.
    :param metric_exprs: Mapping of metric name to the Float64Expression to
        bin (e.g. ``{"lof": ht.lof_oe_upper, "mis": ht.mis_oe_upper}``).
    :param thresholds: Mapping of ``(granularity, metric)`` to a list of threshold
        values in ascending order, built from :func:`compute_percentile_thresholds`
        output (see the example above).
    :param granularities: Granularity names to iterate over (e.g.
        ``["decile", "ventile"]``). Each must appear as the first element of
        at least one key in ``thresholds``.
    :param field_name: Name of the struct field to annotate on ``ht``.
        Default is ``"constraint_bins"``.
    :return: Table with an added struct field containing per-granularity,
        per-metric bin assignments.
    """
    miss = hl.missing(hl.tint32)

    def _bin_expr(
        value_expr: hl.expr.Float64Expression,
        threshold_list: List[float],
    ) -> hl.expr.Int32Expression:
        """
        Count how many thresholds a value meets or exceeds.

        :param value_expr: Metric value to bin.
        :param threshold_list: Boundary values in ascending order.
        :return: Number of boundaries met or exceeded (0 = below all thresholds).
        """
        arr = hl.literal(threshold_list)
        return hl.sum(arr.map(lambda t: hl.int(value_expr >= t)))

    granularities_expr = {
        gran: hl.struct(
            **{
                metric: hl.if_else(
                    hl.is_defined(metric_exprs[metric]),
                    _bin_expr(metric_exprs[metric], thresholds[(gran, metric)]),
                    miss,
                )
                for metric in metric_exprs
            }
        )
        for gran in granularities
    }

    return ht.annotate(**{field_name: hl.struct(**granularities_expr)})


def rank_array_element_metrics(
    ht: hl.Table,
    array_field: str,
    element_value_fn: Callable[
        [hl.expr.StructExpression], Dict[str, hl.expr.Float64Expression]
    ],
    filter_fn: Optional[Callable[[hl.Table], hl.expr.BooleanExpression]] = None,
    bin_granularities: Optional[Dict[str, int]] = None,
    rank_field_prefix: str = "",
) -> hl.Table:
    """
    Rank metrics within array elements and annotate rank structs back.

    **Rank-based binning for array fields**: applies
    :func:`rank_and_assign_bins` independently to each element of an array
    field. For each element, ``element_value_fn`` extracts named metric
    values, which are ranked across rows (optionally on a filtered subset
    via ``filter_fn``). Each array element is then annotated with
    ``{metric_name}_rank`` structs containing rank and bin fields.

    Rows not matching ``filter_fn`` get missing rank annotations.

    This differs from threshold-based binning (see
    :func:`compute_percentile_thresholds` and
    :func:`annotate_bins_by_threshold`), where bins are assigned by
    comparing values against pre-computed boundaries rather than sorted
    position.

    .. note::

        This function evaluates its input when it is called. It caches the full
        input Table (with an added row index) and the extracted metric values,
        reads one row to get the array length, and each
        :func:`rank_and_assign_bins` call it makes caches its own ranked rows.
        Caching avoids recomputing ``ht`` and ``filter_fn`` for every metric, but
        in Hail ``Table.cache()`` writes a checkpoint under the temporary directory
        (see :func:`rank_and_assign_bins`), so this writes a copy of ``ht`` there.
        The copies are not removed when the function returns.

    :param ht: Input Table.
    :param array_field: Name of the array field on ``ht``.
    :param element_value_fn: Function that takes an array element
        (StructExpression) and returns a dict mapping metric names to
        Float64Expressions to rank. Applied identically to every element.
    :param filter_fn: Optional function that takes a Table and returns a
        BooleanExpression to filter rows before ranking. When None, all
        rows are ranked.
    :param bin_granularities: Bin granularities passed to
        :func:`rank_and_assign_bins`.
    :param rank_field_prefix: Prefix for the ``rank`` and ``bin_{name}``
        sub-fields within each rank struct. Passed through to
        :func:`rank_and_assign_bins`. Default is ``""`` (no prefix).
    :return: Table with ``{metric_name}_rank`` structs added to each array
        element. The table is returned with its original key restored.
    """
    logger.warning(
        "rank_array_element_metrics evaluates its input and checkpoints it, and the"
        " ranked values, to the Hail temporary directory (via .cache()); see the"
        " function's docstring."
    )

    original_key = list(ht.key)
    ht = ht.add_index("_rank_idx").key_by("_rank_idx").cache()

    subset_ht = ht.filter(filter_fn(ht)) if filter_fn is not None else ht

    # Extract values to rank using element_value_fn applied via .map().
    subset_ht = subset_ht.select(
        _rank_values=subset_ht[array_field].map(
            lambda elem: hl.struct(**element_value_fn(elem))
        )
    ).cache()

    # Metric names come from the type, but the array length is not part of the type,
    # so the element count is taken from a sample row.
    metric_names = list(subset_ht._rank_values.dtype.element_type.fields)
    sample = subset_ht.take(1)
    if not sample:
        raise ValueError("No rows left to rank after applying 'filter_fn'.")
    n_elements = len(sample[0]._rank_values)

    # Rank each metric within each array element.
    subset_ht = subset_ht.annotate(
        _rank_values=[
            hl.struct(
                **{
                    name: rank_and_assign_bins(
                        subset_ht._rank_values[i][name],
                        bin_granularities,
                        prefix=rank_field_prefix,
                    )
                    for name in metric_names
                }
            )
            for i in range(n_elements)
        ]
    ).cache()

    # Join ranks back to the original table. Use or_missing so that
    # unranked rows get correctly-typed missing rank annotations without
    # needing a manually-constructed missing struct for if_else.
    ht = ht.annotate(_ranks=subset_ht[ht._rank_idx]._rank_values)
    ht = ht.annotate(
        **{
            array_field: [
                ht[array_field][i].annotate(
                    **{
                        f"{name}_rank": hl.or_missing(
                            hl.is_defined(ht._ranks),
                            ht._ranks[i][name],
                        )
                        for name in metric_names
                    }
                )
                for i in range(n_elements)
            ]
        }
    ).drop("_ranks")

    return ht.key_by(*original_key).drop("_rank_idx")
