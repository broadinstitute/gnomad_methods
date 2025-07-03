# noqa: D100

from typing import List, Union

import hail as hl


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


def explode_intervals_to_loci(
    ht: hl.Table,
    keep_intervals: bool = False,
) -> hl.Table:
    """
    Expand intervals to loci.

    :param ht: Hail Table with an interval field.
    :param keep_intervals: If True, keep the original interval as a column in output.
        Default is False.
    :return: Hail Table keyed by loci and intervals as optional field.
    """
    ht = ht.annotate(
        pos=hl.range(ht.interval.start.position, ht.interval.end.position + 1),
    ).explode("pos")
    ht = ht.annotate(
        locus=hl.locus(
            ht.interval.start.contig,
            ht.pos,
            reference_genome=ht.interval.start.dtype.reference_genome,
        ),
    ).key_by("locus")
    return ht.drop("interval", "pos") if not keep_intervals else ht.drop("pos")
