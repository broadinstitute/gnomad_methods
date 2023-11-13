"""Utils module containing generic functions that are useful for adding transcript expression-aware annotations."""
from typing import Callable, Optional

import hail as hl


def summarize_rsem_mt(
    rsem_mt: hl.MatrixTable,
    rsem_expr: hl.expr.NumericExpression,
    tissue_expr: hl.expr.StringExpression,
    summary_agg_func: Optional[Callable] = None,
    tissue_as_row: bool = False,
) -> hl.Table:
    """
    Summarize an RSEM table with ENSTs and ENSGs as rows and samples as columns by tissue.

    The `summary_agg_func` argument allows the user to specify a Hail aggregation
    function to use to summarize the expression by tissue. By default, the median is
    used.

    The output can be returned in one of the following formats (both keyed by
    "transcript_id" and "gene_id"):
      - A Table with an 'rsem' field containing an array of summarized expression
        values by tissue, where the order of tissues in the array is indicated by the
        "tissues" global annotation (`tissue_as_row` set to False).
      - A Table with a row annotation for each tissue containing the summarized tissue
        expression value (`tissue_as_row` set to True).

    :param rsem_mt: MatrixTable of RSEM quantifications.
    :param tissue_expr: Column expression indicating tissue type.
    :param rsem_expr: Entry expression indicating RSEM quantification.
    :param summary_agg_func: Optional aggregation function to use to summarize the RSEM
        values by tissue. Default is None, which will use a median aggregation.
    :param tissue_as_row: If True, return a Table with a row annotation for each tissue
        instead of an array of RSEM values. Default is False.
    :return: Table of summarized transcript expression.
    """
    if summary_agg_func is None:
        summary_agg_func = lambda x: hl.median(hl.agg.collect(x))

    rsem_mt = rsem_mt.group_cols_by(tissue=tissue_expr).aggregate(
        rsem=summary_agg_func(rsem_expr)
    )

    if tissue_as_row:
        rsem_ht = rsem_mt.rename({"rsem": ""}).make_table()
    else:
        rsem_ht = rsem_mt.localize_entries(
            columns_array_field_name="tissues", entries_array_field_name="rsem"
        )
        rsem_ht = rsem_ht.annotate(rsem=rsem_ht.rsem.map(lambda x: x.rsem))
        rsem_ht = rsem_ht.annotate_globals(
            tissues=rsem_ht.tissues.map(lambda x: x.tissue)
        )

    return rsem_ht.key_by("transcript_id", "gene_id")
