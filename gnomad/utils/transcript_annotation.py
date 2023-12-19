"""Utils module containing generic functions that are useful for adding transcript expression-aware annotations."""
import logging
from typing import Callable, List, Optional, Tuple, Union

import hail as hl

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("transcript_annotation_utils")
logger.setLevel(logging.INFO)


def summarize_transcript_expression(
    mt: hl.MatrixTable,
    transcript_expression_expr: Union[
        hl.expr.NumericExpression, str
    ] = "transcript_tpm",
    tissue_expr: Union[hl.expr.StringExpression, str] = "tissue",
    summary_agg_func: Optional[Callable] = None,
) -> hl.Table:
    """
    Summarize a transcript expression MatrixTable by transcript, gene, and tissue.

    The `summary_agg_func` argument allows the user to specify a Hail aggregation
    function to use to summarize the expression by tissue. By default, the median is
    used.

    The returned Table has a row annotation for each tissue containing the summarized
    tissue expression value.

    :param mt: MatrixTable of transcript (rows) expression quantifications (entry) by
        sample (columns).
    :param transcript_expression_expr: Entry expression indicating transcript expression
        quantification. Default is 'transcript_tpm'.
    :param tissue_expr: Column expression indicating tissue type. Default is 'tissue'.
    :param summary_agg_func: Optional aggregation function to use to summarize the
        transcript expression quantification by tissue. Example: `hl.agg.mean`. Default
        is None, which will use a median aggregation.
    :return: A Table of summarized transcript expression by tissue
    """
    if summary_agg_func is None:
        summary_agg_func = lambda x: hl.median(hl.agg.collect(x))

    if isinstance(transcript_expression_expr, str):
        transcript_expression_expr = mt[transcript_expression_expr]

    if isinstance(tissue_expr, str):
        tissue_expr = mt[tissue_expr]

    mt = mt.group_cols_by(tissue=tissue_expr).aggregate(
        tx=summary_agg_func(transcript_expression_expr)
    )

    transcript_ht = mt.rename({"tx": ""}).make_table()
    transcript_ht = transcript_ht.key_by("transcript_id", "gene_id")

    return transcript_ht


def tissue_expression_ht_to_array(
    ht: hl.Table,
    tissues: Optional[List[str]] = None,
    tissues_to_filter: Optional[List[str]] = None,
) -> hl.Table:
    """
    Convert a Table with a row annotation for each tissue to a Table with tissues as an array.

    The output is a Table with a field 'tissue_expression' containing an array of
    summarized expression values by tissue, where the order of tissues in the array is
    indicated by the "tissues" global annotation.

    :param ht: Table with a row annotation for each tissue.
    :param tissues: Optional list of tissues to keep in the 'tissue_expression' array.
        Default is all non-key rows in the Table.
    :param tissues_to_filter: Optional list of tissues to exclude from the tissue
        expression array.
    :return: Table with a field 'tissue_expression' containing an array of summarized
        expression values by tissue.
    """
    if tissues is None:
        tissues = list(ht.row_value)

    if tissues_to_filter is not None:
        logger.info("Filtering tissues: %s", tissues_to_filter)
        tissues = [t for t in tissues if t not in tissues_to_filter]

    ht = ht.select_globals(tissues=tissues)
    ht = ht.select(tissue_expression=[ht[t] for t in tissues])

    return ht


def get_expression_proportion(
    transcript_ht: hl.Table,
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

    # Calculate the sum of transcript expression by gene per tissue.
    gene_ht = transcript_ht.group_by("gene_id").aggregate(
        expression_sum=hl.agg.array_sum(transcript_ht.tissue_expression)
    )

    # Calculate the mean expression proportion across tissues.
    transcript_ht = transcript_ht.annotate(
        gene_exp_sum=gene_ht[transcript_ht.gene_id].expression_sum,
        exp_prop=hl.map(
            lambda x, y: x / y,
            transcript_ht.tissue_expression,
            gene_ht[transcript_ht.gene_id].expression_sum,
        ),
        exp_prop_mean=hl.mean(
            hl.filter(
                lambda e: ~hl.is_nan(e),
                hl.map(
                    lambda x, y: x / y,
                    transcript_ht.tissue_expression,
                    gene_ht[transcript_ht.gene_id].expression_sum,
                ),
            )
        ),
    )

    return transcript_ht
