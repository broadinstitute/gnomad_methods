"""Utils module containing generic functions that are useful for adding transcript expression-aware annotations."""
from typing import Callable, Optional, Tuple

import hail as hl


def summarize_rsem_mt(
    rsem_mt: hl.MatrixTable,
    rsem_expr: hl.expr.NumericExpression,
    tissue_expr: hl.expr.StringExpression,
    summary_agg_func: Optional[Callable] = None,
    tissue_as_row: bool = False,
) -> Tuple[hl.Table, hl.Table]:
    """
    Summarize an RSEM table with ENSTs and ENSGs as rows and samples as columns by tissue.

    The `summary_agg_func` argument allows the user to specify a Hail aggregation
    function to use to summarize the expression by tissue. By default, the median is
    used.

    .. note::

        The outputs can be returned in one of the following formats:

        - A Table with a field containing an array of summarized expression
          values by tissue, where the order of tissues in the array is indicated by
          the "tissues" global annotation (`tissue_as_row` set to False).
        - A Table with a row annotation for each tissue containing the summarized
          tissue expression value (`tissue_as_row` set to True).

    :param rsem_mt: MatrixTable of RSEM quantifications.
    :param tissue_expr: Column expression indicating tissue type.
    :param rsem_expr: Entry expression indicating RSEM quantification.
    :param summary_agg_func: Optional aggregation function to use to summarize the RSEM
        values by tissue. Default is None, which will use a median aggregation.
    :param tissue_as_row: If True, return a Table with a row annotation for each tissue
        instead of an array of RSEM values. Default is False.
    :return: A Table of summarized transcript expression and a Table of summarized
        gene expression.
    """
    if summary_agg_func is None:
        summary_agg_func = lambda x: hl.median(hl.agg.collect(x))

    rsem_mt = rsem_mt.group_cols_by(tissue=tissue_expr).aggregate(
        transcript_expression=summary_agg_func(rsem_expr)
    )

    if tissue_as_row:
        transcript_ht = rsem_mt.rename({"transcript_expression": ""}).make_table()
        gene_ht = transcript_ht.key_by("gene_id").drop("transcript_id")
        tissues = list(gene_ht.row)
        tissues.remove("gene_id")
        gene_ht = gene_ht.group_by(*gene_ht.key).aggregate(
            **{tissue: hl.agg.sum(gene_ht[tissue]) for tissue in tissues}
        )
    else:
        transcript_ht = rsem_mt.localize_entries(
            columns_array_field_name="tissues",
            entries_array_field_name="transcript_expression",
        )
        transcript_ht = transcript_ht.annotate(
            transcript_expression=transcript_ht.transcript_expression.map(
                lambda x: x.transcript_expression
            )
        )
        transcript_ht = transcript_ht.annotate_globals(
            tissues=transcript_ht.tissues.map(lambda x: x.tissue)
        )
        gene_ht = transcript_ht.group_by(transcript_ht.gene_id).aggregate(
            gene_expression=hl.agg.array_sum(transcript_ht.transcript_expression)
        )

    return transcript_ht.key_by("transcript_id", "gene_id"), gene_ht.key_by("gene_id")
