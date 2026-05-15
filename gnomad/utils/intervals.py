# noqa: D100

import logging
from typing import List, Optional, Union

import hail as hl

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def sort_intervals(intervals: List[hl.Interval]):
    """
    Sort an array of intervals by start contig, then start position, then end contig, then end position.

    :param intervals: Intervals to sort
    :return: Sorted interval list
    """
    return sorted(
        intervals,
        key=lambda interval: (
            interval.start.reference_genome.contigs.index(interval.start.contig),
            interval.start.position,
            interval.end.reference_genome.contigs.index(interval.end.contig),
            interval.end.position,
        ),
    )


def union_intervals(intervals: List[hl.Interval], is_sorted: bool = False):
    """
    Generate a list with the union of all intervals in the input list by merging overlapping intervals.

    :param intervals: Intervals to merge
    :param is_sorted: If set, assumes intervals are already sorted, otherwise will sort.
    :return: List of merged intervals
    """
    sorted_intervals = intervals if is_sorted else sort_intervals(intervals)
    merged_intervals = sorted_intervals[:1]
    for interval in sorted_intervals[1:]:
        if merged_intervals[-1].start.contig == interval.start.contig:
            if merged_intervals[-1].end.position < interval.end.position:
                if interval.start.position <= merged_intervals[-1].end.position:
                    merged_intervals[-1] = hl.Interval(
                        merged_intervals[-1].start, interval.end
                    )
                else:
                    merged_intervals.append(interval)
        else:
            merged_intervals.append(interval)

    return merged_intervals


def interval_length(interval: hl.Interval) -> int:
    """
    Return the total number of bases in an Interval.

    :param interval: Input interval
    :return: Total length of the interval
    """
    if interval.start.contig != interval.end.contig:
        ref = interval.start.reference_genome
        return (
            ref.contig_length(interval.start.contig)
            - interval.start.position
            + sum(
                ref.contig_length(contig)
                for contig in ref.contigs[
                    ref.contigs.index(interval.start.contig)
                    + 1 : ref.contigs.index(interval.end.contig)
                ]
            )
            + interval.end.position
        )
    else:
        return interval.end.position - interval.start.position


def pad_intervals(
    intervals: Union[
        hl.expr.IntervalExpression,
        hl.Interval,
        List[hl.Interval],
    ],
    padding_bp: int,
) -> Union[hl.expr.IntervalExpression, List[hl.expr.IntervalExpression]]:
    """
    Add padding to interval(s).

    :param intervals: Interval(s) to add padding to. This can be a single Interval,
        a list of Intervals, or an IntervalExpression.
    :param padding_bp: Number of base pairs to add to each side of the interval.
    :return: Interval(s) with padding added.
    """

    def _add_padding(
        interval: Union[hl.expr.IntervalExpression, hl.Interval],
    ) -> hl.expr.IntervalExpression:
        """
        Add padding to an interval.

        :param interval: Interval to add padding to.
        :return: Interval with padding added.
        """
        return hl.locus_interval(
            interval.start.contig,
            interval.start.position - padding_bp,
            interval.end.position + padding_bp,
            includes_start=interval.includes_start,
            includes_end=interval.includes_end,
            reference_genome=interval.start.dtype.reference_genome,
        )

    if isinstance(intervals, List):
        return [_add_padding(i) for i in intervals]
    else:
        return _add_padding(intervals)


def interval_to_pos_range(
    interval_expr: hl.expr.IntervalExpression,
) -> hl.expr.ArrayExpression:
    """
    Convert an IntervalExpression to an ArrayExpression of integer positions.

    Handles both inclusive and exclusive interval endpoints.

    :param interval_expr: IntervalExpression to convert.
    :return: ArrayExpression of integer positions within the interval.
    """
    start = hl.if_else(
        interval_expr.includes_start,
        interval_expr.start.position,
        interval_expr.start.position + 1,
    )
    end = hl.if_else(
        interval_expr.includes_end,
        interval_expr.end.position + 1,
        interval_expr.end.position,
    )
    return hl.range(start, end)


def explode_intervals_to_loci(
    intervals: Union[hl.Table, hl.expr.IntervalExpression, hl.expr.ArrayExpression],
    interval_field: Optional[str] = None,
    keep_intervals: Optional[bool] = False,
    deduplicate: bool = True,
    flatten: bool = True,
) -> Union[hl.Table, hl.expr.ArrayExpression]:
    """
    Expand interval(s) to loci.

    If input is a Table, function will expand intervals to loci and key Table by loci.

    If input is an IntervalExpression or an ArrayExpression of IntervalExpressions,
    function will return an ArrayExpression containing all loci within the input
    interval(s).

    .. warning::
        - Overlapping intervals will produce duplicate loci. Use ``deduplicate=True``
          (default) to remove them.
        - When ``keep_intervals=True`` on a Table input, deduplication is not possible
          because duplicate rows with different interval annotations may exist; a warning
          is displayed instead.
        - Caution when using this function on very large intervals (e.g., whole
          chromosomes), as it will create extremely large arrays, which may cause
          performance issues.

    Note that intervals that cross chromosomes are currently not supported.

    :param intervals: Table, IntervalExpression, or ArrayExpression of
        IntervalExpressions.
    :param interval_field: Name of the interval field. Only required if input is a Hail
        Table. Default is None.
    :param keep_intervals: If True, keep the original intervals as a column in output.
        Only applies if input is a Hail Table. Default is False.
    :param deduplicate: If True, remove duplicate loci produced by overlapping
        intervals. For Table input with ``keep_intervals=True``, deduplication is
        skipped with a warning. For an ArrayExpression of IntervalExpressions, the
        returned ArrayExpression will have duplicate loci removed. Has no effect when
        ``flatten=False``. Default is True.
    :param flatten: If True, flatten the per-interval loci arrays into a single
        ArrayExpression. Only applies when input is an ArrayExpression of
        IntervalExpressions. Default is True.
    :return: If input is a Hail Table, returns exploded Table keyed by locus. If input
        is an IntervalExpression or ArrayExpression of IntervalExpressions, returns
        ArrayExpression containing loci within input interval(s). If input is an
        ArrayExpression and ``flatten=False``, returns a nested ArrayExpression of
        per-interval loci arrays.
    """
    assert (
        isinstance(intervals, hl.Table)
        or isinstance(intervals, hl.expr.IntervalExpression)
        or isinstance(intervals, hl.expr.ArrayExpression)
    ), (
        "Input must be a Table, IntervalExpression, or ArrayExpression of"
        " IntervalExpressions!"
    )

    if isinstance(intervals, hl.Table) and (
        not interval_field or keep_intervals is None
    ):
        raise ValueError(
            "`interval_field` and `keep_intervals` must be defined if input is a Table!"
        )
    if isinstance(intervals, hl.Table):
        if interval_field not in intervals.row:
            raise ValueError(
                "`interval_field` must be an annotation present on input Table!"
            )
        if not isinstance(intervals[interval_field], hl.expr.IntervalExpression):
            raise ValueError(
                f"`interval_field` '{interval_field}' must be an IntervalExpression"
                " in the input Table!"
            )

    def _make_loci_array(interval_expr):
        return interval_to_pos_range(interval_expr).map(
            lambda pos: hl.locus(
                interval_expr.start.contig,
                pos,
                reference_genome=interval_expr.start.dtype.reference_genome,
            )
        )

    if isinstance(intervals, hl.expr.ArrayExpression):
        logger.info(
            "Input is an ArrayExpression of IntervalExpressions, so function will"
            " return an ArrayExpression of loci within all input intervals."
        )
        if not deduplicate and flatten:
            logger.warning(
                "Overlapping intervals in the input array may produce duplicate loci"
                " in the returned ArrayExpression. Set `deduplicate=True` to remove"
                " them."
            )

        loci_arrays = intervals.map(lambda i: _make_loci_array(i))
        if flatten:
            result = hl.flatten(loci_arrays)
            if deduplicate:
                result = hl.array(hl.set(result))
            result = hl.sorted(result)
            return result
        else:
            return loci_arrays

    intervals_expr = (
        intervals
        if isinstance(intervals, hl.expr.IntervalExpression)
        else intervals[interval_field]
    )
    loci_array = _make_loci_array(intervals_expr)

    if isinstance(intervals, hl.Table):
        intervals = intervals.annotate(_loci=loci_array).explode("_loci")
        intervals = intervals.key_by(locus=intervals._loci)

        fields_to_drop = ["_loci"]
        if not keep_intervals:
            fields_to_drop.append(interval_field)

        intervals = intervals.drop(*fields_to_drop)

        if deduplicate:
            if keep_intervals:
                logger.warning(
                    "`deduplicate=True` has no effect when `keep_intervals=True`"
                    " because rows with different interval annotations cannot be safely"
                    " collapsed. Duplicate loci may be present in the output if"
                    " intervals overlap."
                )
            else:
                intervals = intervals.distinct()
                logger.warning(
                    "The `distinct()` call will arbitrarily select one row for each"
                    " duplicated locus. If the input table has annotations such as gene"
                    " or transcript ID, the values retained for duplicated loci are not"
                    " deterministic."
                )

        return intervals

    logger.info(
        "Input is an IntervalExpression, so function will return an ArrayExpression of"
        " loci within the input interval."
    )
    return loci_array
