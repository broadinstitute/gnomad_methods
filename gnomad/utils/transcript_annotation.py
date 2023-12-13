"""Utils module containing generic functions that are useful for adding transcript expression-aware annotations."""
from typing import Callable, List, Optional, Tuple

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


def get_expression_proportion(
    transcript_ht: hl.Table,
    gene_ht: hl.Table,
    tissues_to_filter: Optional[List[str]] = None,
) -> hl.Table:
    """
    Calculate the proportion of expression of transcript to gene per tissue.

    :param transcript_ht: Table of summarized transcript expression by tissue.
    :param gene_ht: Table of summarized gene expression by tissue.
    :param tissues_to_filter: Optional list of tissues to filter out
    :return: Table with expression proportion of transcript to gene per tissue
        and mean expression proportion across tissues.
    """
    transcript_ht = tissue_expression_ht_to_array(
        transcript_ht, tissues_to_filter=tissues_to_filter
    )
    gene_ht = tissue_expression_ht_to_array(
        gene_ht, tissues=hl.eval(transcript_ht.tissues)
    )

    # Join the transcript expression table and gene expression table.
    transcript_ht = transcript_ht.annotate(
        gene_expression=gene_ht[transcript_ht.gene_id].tissue_expression
    )

    # Calculate the proportion of expression of transcript to gene per tissue.
    transcript_ht = transcript_ht.annotate(
        exp_prop=hl.or_else(
            transcript_ht.transcript_expression / transcript_ht.gene_expression,
            hl.empty_array(hl.tfloat64),
        ),
    )
    # Calculate the mean expression proportion across tissues.
    transcript_ht = transcript_ht.annotate(
        exp_prop_mean=hl.mean(
            hl.filter(lambda e: ~hl.is_nan(e), transcript_ht.exp_prop),
        )
    )

    return transcript_ht
