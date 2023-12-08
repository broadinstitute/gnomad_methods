"""Utils module containing generic functions that are useful for adding transcript expression-aware annotations."""
from typing import Callable, List, Optional, Tuple

import hail as hl

from gnomad.utils.vep import add_most_severe_consequence_to_consequence


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
    tissue_as_row: bool = False,
) -> hl.Table:
    """
    Calculate the proportion of expression of transcript to gene per tissue.

    :param transcript_ht: Table of summarized transcript expression by tissue.
    :param gene_ht: Table of summarized gene expression by tissue.
    :param tissues_to_filter: Optional list of tissues to filter out.
    :param tissue_as_row: If True, the input Table is with a row annotation for each
        tissue instead of an array of values. Default is False.
    :return: Table with expression proportion of transcript to gene per tissue
        and mean expression proportion across tissues.
    """
    if tissue_as_row:
        tissues1 = list(gene_ht.row)
        tissues1.remove("gene_id")
        gene_ht = gene_ht.select(
            gene_expression=hl.array([gene_ht[tissue] for tissue in tissues1])
        ).annotate_globals(tissues=tissues1)

        tissues2 = list(transcript_ht.row)
        tissues2.remove("transcript_id")
        tissues2.remove("gene_id")
        transcript_ht = transcript_ht.select(
            transcript_expression=hl.array(
                [transcript_ht[tissue] for tissue in tissues2]
            )
        ).annotate_globals(tissues=tissues2)

    # Join the transcript expression table and gene expression table
    transcript_ht = transcript_ht.annotate(
        gene_expression=gene_ht[transcript_ht.gene_id].gene_expression
    )

    if tissues_to_filter is not None:
        tissues_indices = []
        for i, t in enumerate(hl.eval(transcript_ht.tissues)):
            if t not in tissues_to_filter:
                tissues_indices.append(i)

        transcript_ht = transcript_ht.annotate(
            transcript_expression=hl.array(
                [transcript_ht.transcript_expression[i] for i in tissues_indices]
            ),
            gene_expression=hl.array(
                [transcript_ht.gene_expression[i] for i in tissues_indices]
            ),
        )

        transcript_ht = transcript_ht.select_globals(
            tissues=hl.array([transcript_ht.tissues[i] for i in tissues_indices])
        )

    # Calculate the proportion of expression of transcript to gene per tissue
    transcript_ht = transcript_ht.annotate(
        exp_prop=hl.or_else(
            transcript_ht.transcript_expression / transcript_ht.gene_expression,
            hl.empty_array(hl.tfloat64),
        ),
    )
    # Calculate the mean expression proportion across tissues
    transcript_ht = transcript_ht.annotate(
        exp_prop_mean=hl.mean(
            hl.filter(lambda e: ~hl.is_nan(e), transcript_ht.exp_prop),
        )
    )

    return transcript_ht


def tx_annotate_variants(
    variant_ht: hl.Table,
    tx_ht: hl.Table,
    tx_annotation_type: str = "expression",
    filter_to_protein_coding: bool = True,
    filter_to_genes: List[str] = None,
    filter_to_csqs: List[str] = None,
    filter_to_homs: bool = False,
) -> hl.Table:
    """
    Annotate variants with transcript-based expression values or expression proportion from GTEx.

    :param variant_ht: Table of variants to annotate, it should contain at
        least the following nested fields: `vep.transcript_consequences`,
        `freq`.
    :param tx_ht: Table of transcript expression information.
    :param tx_annotation_type: Type of transcript annotation to add. Options are
        'expression' or 'expression_proportion'. Default is 'expression'.
    :param filter_to_protein_coding: If True, filter to protein coding
        transcripts. Default is True.
    :param filter_to_genes: List of genes to filter to.
    :param filter_to_csqs: List of consequence terms to filter to.
    :param filter_to_homs: If True, filter to variants with at least one
        homozygote in `freq` field. Default is False.
    :return: MatrixTable with transcript expression information annotated
    """
    # GTEx data has transcript IDs without version numbers, so we need to
    # remove the version numbers from the transcript IDs in the variant table
    tx_ht = tx_ht.key_by()
    tx_ht = tx_ht.annotate(transcript_id=tx_ht.transcript_id.split("\\.")[0])
    tx_ht = tx_ht.key_by(tx_ht.transcript_id)

    variant_ht = variant_ht.annotate(
        vep=variant_ht.vep.annotate(
            transcript_consequences=variant_ht.vep.transcript_consequences.map(
                add_most_severe_consequence_to_consequence
            )
        )
    )

    # Explode the transcript consequences to be able to key by transcript ID
    variant_ht = variant_ht.explode(variant_ht.vep.transcript_consequences)

    if filter_to_protein_coding:
        variant_ht = variant_ht.filter(
            variant_ht.vep.transcript_consequences.biotype == "protein_coding"
        )

    if filter_to_genes:
        variant_ht = variant_ht.filter(
            hl.literal(filter_to_genes).contains(
                variant_ht.vep.transcript_consequences.gene_id
            )
        )

    if filter_to_csqs:
        variant_ht = variant_ht.filter(
            hl.literal(filter_to_csqs).contains(
                variant_ht.vep.transcript_consequences.most_severe_consequence
            )
        )

    if filter_to_homs:
        variant_ht = variant_ht.filter(variant_ht.freq[0].homozygote_count > 0)

    # Annotate the variant table with the transcript expression information
    variant_ht = variant_ht.annotate(
        tx_data=tx_ht[variant_ht.vep.transcript_consequences.transcript_id]
    )

    # Aggregate the transcript expression information by gene, csq, etc.
    tx_annot_ht = variant_ht.group_by(
        locus=variant_ht.locus,
        alleles=variant_ht.alleles,
        gene_id=variant_ht.vep.transcript_consequences.gene_id,
        gene_symbol=variant_ht.vep.transcript_consequences.gene_symbol,
        csq=variant_ht.vep.transcript_consequences.most_severe_consequence,
        lof=variant_ht.vep.transcript_consequences.lof,
        lof_flags=variant_ht.vep.transcript_consequences.lof_flags,
    ).aggregate(
        tx_annotation=hl.agg.array_sum(
            variant_ht.tx_data.transcript_expression
            if tx_annotation_type == "expression"
            else variant_ht.tx_data.exp_prop
        )
    )

    # Remove unnecessary global annotations and add tissue and
    # expression_type in global annotations
    tx_annot_ht = tx_annot_ht.select_globals().annotate_globals(
        tissues=tx_ht.index_globals().tissues,
        tx_annotation_type=tx_annotation_type,
    )

    tx_annot_ht = tx_annot_ht.key_by(tx_annot_ht.locus, tx_annot_ht.alleles)

    return tx_annot_ht
