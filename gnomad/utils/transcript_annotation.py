"""Utils module containing generic functions that are useful for adding transcript expression-aware annotations."""
import logging
from typing import Callable, List, Optional, Tuple, Union

import hail as hl

from gnomad.resources.grch37 import gencode
from gnomad.utils.filtering import filter_gencode_to_cds
from gnomad.utils.vep import (
    CSQ_CODING_ALL_IMPACT,
    CSQ_CODING_HIGH_IMPACT,
    CSQ_ORDER,
    SPLICE_CSQS,
    explode_by_vep_annotation,
    filter_vep_transcript_csqs,
    process_consequences,
)

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

    The returned Table has a row annotation for each tissue containing a struct with the
    summarized tissue expression value ('transcript_expression') and the proportion of
    expression of transcript to gene per tissue ('expression_proportion').

    Returned Table Schema example::

        Row fields:
            'transcript_id': str
            'gene_id': str
            'tissue_1': struct {
              transcript_expression: float64,
              expression_proportion: float64
            }
            'tissue_2': struct {
              transcript_expression: float64,
              expression_proportion: float64
            }

        Key: ['transcript_id', 'gene_id']

    :param mt: MatrixTable of transcript (rows) expression quantifications (entry) by
        sample (columns).
    :param transcript_expression_expr: Entry expression indicating transcript expression
        quantification. Default is 'transcript_tpm'.
    :param tissue_expr: Column expression indicating tissue type. Default is 'tissue'.
    :param summary_agg_func: Optional aggregation function to use to summarize the
        transcript expression quantification by tissue. Example: `hl.agg.mean`. Default
        is None, which will use a median aggregation.
    :return: A Table of summarized transcript expression by tissue.
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
    ht = mt.rename({"tx": ""}).make_table().key_by("transcript_id", "gene_id")

    # Annotate with the proportion of expression of transcript to gene per tissue.
    ht = ht.annotate(expression_proportion=get_expression_proportion(ht))
    ht = ht.select(
        **{
            t: hl.struct(
                transcript_expression=ht[t],
                expression_proportion=ht.expression_proportion[t],
            )
            for t in ht.expression_proportion
        }
    )

    return ht


def get_expression_proportion(ht: hl.Table) -> hl.expr.StructExpression:
    """
    Calculate the proportion of expression of transcript to gene per tissue.

    :param ht: Table of summarized transcript expression by tissue.
    :return: StructExpression containing the proportion of expression of transcript to
        gene per tissue.
    """
    tissues = list(ht.row_value)

    # Calculate the sum of transcript expression by gene per tissue.
    gene_ht = ht.group_by("gene_id").aggregate(
        **{tissue: hl.agg.sum(ht[tissue]) for tissue in tissues}
    )

    # Return the proportion of expression of transcript to gene per tissue.
    gene = gene_ht[ht.gene_id]
    return hl.struct(
        **{
            tissue: hl.utils.misc.divide_null(ht[tissue], gene[tissue])
            for tissue in tissues
        }
    )


def filter_expression_ht_by_tissues(
    ht: hl.Table,
    tissues_to_keep: Optional[List[str]] = None,
    tissues_to_filter: Optional[List[str]] = None,
) -> hl.Table:
    """
    Filter a Table with a row annotation for each tissue to only include specified tissues.

    :param ht: Table with a row annotation for each tissue.
    :param tissues_to_keep: Optional list of tissues to keep in the Table. Default is
        all non-key rows in the Table.
    :param tissues_to_filter: Optional list of tissues to exclude from the Table.
    :return: Table with only specified tissues.
    """
    if tissues_to_keep is None and tissues_to_filter is None:
        logger.info(
            "No tissues_to_keep or tissues_to_filter specified. Returning input Table."
        )
        return ht

    if tissues_to_keep is None:
        tissues = list(ht.row_value)

    if tissues_to_filter is not None:
        logger.info("Filtering tissues: %s", tissues_to_filter)
        tissues = [t for t in tissues if t not in tissues_to_filter]

    ht = ht.select(*tissues)

    return ht


def tissue_expression_ht_to_array(
    ht: hl.Table,
    tissues_to_keep: Optional[List[str]] = None,
    tissues_to_filter: Optional[List[str]] = None,
    annotations_to_extract: Optional[Union[Tuple[str], List[str]]] = (
        "transcript_expression",
        "expression_proportion",
    ),
) -> hl.Table:
    """
    Convert a Table with a row annotation for each tissue to a Table with tissues as an array.

    The output is a Table with one of the two formats:
        - An annotation of 'tissue_expression' containing an array of structs by
          tissue, where each element of the array is the Table's row value for a given
          tissue.

            Example::

                tissue_expression': array<struct {
                    transcript_expression: float64,
                    expression_proportion: float64
                }>

        - One array annotation for each field defined in the 'annotations_to_extract'
          argument, where each array is an array of the given field values by tissue.

            Example::

                'transcript_expression': array<float64>
                'expression_proportion': array<float64>

    The order of tissues in the array is indicated by the "tissues" global annotation.

    :param ht: Table with a row annotation for each tissue.
    :param tissues_to_keep: Optional list of tissues to keep in the tissue expression
        array. Default is all non-key rows in the Table.
    :param tissues_to_filter: Optional list of tissues to exclude from the
        tissue expression array.
    :param annotations_to_extract: Optional list of tissue struct fields to extract
        into top level array annotations. If None, the returned Table will contain a
        single top level annotation 'tissue_expression' that contains an array of
        structs by tissue. Default is ('transcript_expression', 'expression_proportion').
    :return: Table with requested tissue struct annotations pulled into arrays of
        tissue values and a 'tissues' global annotation indicating the order of tissues
        in the arrays.
    """
    ht = filter_expression_ht_by_tissues(ht, tissues_to_keep, tissues_to_filter)

    tissues = list(ht.row_value)
    ht = ht.select_globals(tissues=tissues)
    ht = ht.select(tissue_expression=[ht[t] for t in tissues])

    if annotations_to_extract is not None:
        ht = ht.select(
            **{
                a: ht.tissue_expression.map(lambda x: x[a])
                for a in annotations_to_extract
            }
        )

    return ht


def preprocess_variants_for_tx(
    ht: hl.Table,
    filter_to_cds: bool = True,
    filter_to_genes: Optional[List[str]] = None,
    match_by_gene_symbol: bool = False,
    filter_to_csqs: Optional[List[str]] = None,
    ignore_splicing: bool = True,
    filter_to_protein_coding: bool = True,
    vep_root: str = "vep",
) -> hl.Table:
    """
    Prepare a Table of variants with vep transcript consequences for annotation.

    :param ht: Table of variants with 'vep' annotations.
    :param cds_intervals: Optional Table of CDS intervals. Default is None.
    :param filter_to_cds: Whether to filter to CDS regions. Default is True.
    :param filter_to_genes: Optional list of genes to filter to. Default is None.
    :param match_by_gene_symbol: Whether to match by gene symbol instead of gene ID.
        Default is False.
    :param filter_to_csqs: Optional list of consequences to filter to. Default is None.
    :param ignore_splicing: If True, ignore splice variants. Default is True.
    :param filter_to_protein_coding: Whether to filter to protein coding transcripts.
        Default is True.
    :param vep_root: Name used for root VEP annotation. Default is 'vep'.
    :return: Table of variants with preprocessed/filtered transcript consequences
        prepared for annotation.
    """
    if filter_to_cds:
        logger.info("Filtering to CDS regions...")
        cds = filter_gencode_to_cds(gencode)
        ht = ht.filter(hl.is_defined(cds[ht.locus]))

    keep_csqs = True
    if ignore_splicing:
        if filter_to_csqs is not None:
            filter_to_csqs = [csq for csq in filter_to_csqs if csq not in SPLICE_CSQS]
        else:
            filter_to_csqs = SPLICE_CSQS
            keep_csqs = False

    if filter_to_csqs is not None:
        logger.info("Adding most severe consequence to VEP transcript consequences...")
        ht = process_consequences(ht, vep_root=vep_root)

    return filter_vep_transcript_csqs(
        ht,
        vep_root=vep_root,
        synonymous=False,
        canonical=False,
        protein_coding=filter_to_protein_coding,
        csqs=filter_to_csqs,
        keep_csqs=keep_csqs,
        genes=filter_to_genes,
        match_by_gene_symbol=match_by_gene_symbol,
    )


def tx_annotate_variants(
    ht: hl.Table,
    tx_ht: hl.Table,
    tissues_to_filter: Optional[List[str]] = None,
    vep_root: str = "vep",
    vep_annotation: str = "transcript_consequences",
) -> hl.Table:
    """
    Annotate variants with transcript-based expression values or expression proportion from GTEx.

    :param ht: Table of variants to annotate, it should contain at least the following
        nested fields: `vep.transcript_consequences`, `freq`.
    :param tx_ht: Table of transcript expression information.
    :param tissues_to_filter: Optional list of tissues to exclude from the output.
    :param vep_root: Name used for root VEP annotation. Default is 'vep'.
    :param vep_annotation: Name of annotation in 'vep' annotation,
        one of the processed consequences: ["transcript_consequences",
        "worst_csq_by_gene", "worst_csq_for_variant",
        "worst_csq_by_gene_canonical", "worst_csq_for_variant_canonical"].
        For example, if you want to annotate each variant with the worst
        consequence in each gene it falls on and the transcript expression,
        you would use "worst_csq_by_gene". Default is "transcript_consequences".
    :return: Input Table with transcript expression information annotated.
    """
    # Filter to tissues of interest.
    tx_ht = filter_expression_ht_by_tissues(tx_ht, tissues_to_filter=tissues_to_filter)
    tissues = list(tx_ht.row_value)

    # Calculate the mean expression proportion across all tissues.
    tx_ht = tx_ht.annotate(
        exp_prop_mean=hl.mean([tx_ht[t].expression_proportion for t in tissues])
    )

    # Explode the processed transcript consequences to be able to key by
    # transcript ID.
    ht = explode_by_vep_annotation(ht, vep_annotation=vep_annotation, vep_root=vep_root)
    ht = ht.transmute(
        **ht[vep_annotation],
        **tx_ht[ht[vep_annotation].transcript_id, ht[vep_annotation].gene_id],
    )
    ht = ht.annotate_globals(tissues=tissues)

    return ht


def tx_aggregate_variants(
    ht: hl.Table,
    additional_grouping: bool = False,
    additional_group_by: Optional[Union[Tuple[str], List[str]]] = (
        "alleles",
        "gene_symbol",
        "most_severe_consequence",
        "lof",
        "lof_flags",
    ),
) -> hl.Table:
    """
    Aggregate transcript-based expression values or expression proportion from GTEx.

    :param ht: Table of variants annotated with transcript expression information.
    :param additional_grouping: Whether to group by additional fields before sum
        aggregation. Default is False.
    :param additional_group_by: Optional list of additional fields to group by before
        sum aggregation.
    :return: Table of variants with transcript expression information aggregated.
    """
    tissues = hl.eval(ht.tissues)

    grouping = (
        ["locus", "gene_id"] + list(additional_group_by)
        if additional_grouping
        else ["locus", "gene_id"]
    )

    # Aggregate the transcript expression information by locus, gene_id and
    # annotations in additional_group_by.
    ht = ht.group_by(*grouping).aggregate(
        exp_prop_mean=hl.agg.sum(ht.exp_prop_mean),
        **{t: hl.struct(**{a: hl.agg.sum(ht[t][a]) for a in ht[t]}) for t in tissues},
    )

    ht = ht.key_by(ht.locus, ht.alleles) if additional_grouping else ht.key_by(ht.locus)

    return ht


def keep_loftee_for_high_impact(ht: hl.expr) -> hl.Table:
    """
    Keep the LOFTEE annotations for high_impact coding variants after transcript expression annotation.

    :param ht: Table of variants annotated with transcript expression
        information.
    :return: Table of variants annotated with transcript expression information
        filtered to keep LOFTEE annotations for high_impact coding variants.
    """
    csq_high = hl.literal(set(CSQ_CODING_HIGH_IMPACT))
    ht = ht.annotate(
        lof=hl.if_else(
            csq_high.contains(ht.most_severe_consequence),
            ht.lof,
            hl.missing(hl.tstr),
        ),
        lof_flag=hl.if_else(
            csq_high.contains(ht.most_severe_consequence),
            ht.lof_flags,
            hl.missing(hl.tstr),
        ),
    )
    return ht


def get_worst_csq_per_variant(
    ht: hl.Table, keep_loftee_for_csq_high: bool = False
) -> hl.Table:
    """
    Get the worst consequence of each variant (keyed by locus and alleles) or each locus (keyed by locus).

    :param ht: Table of variants annotated with transcript expression information.
    :param keep_loftee_for_csq_high: If False, include "OS" (Other Splice) annotations
        from 'lof' of LOFTEE. "OS" could be annotated in protein_coding regions on any
        variant except the ones in UTRs, stop_gained and frameshift_variant according to
        LOFTEE settings. "OS" showed to be more constraint than "LC" variants. Default is
        False.
    :return: Table of variants annotated with transcript expression information with
        annotations of the worst consequence of each variant or locus.
    """

    def _get_csq_order(keep_loftee_for_csq_high: bool) -> list:
        """
        Get the order of consequence terms.

        :param keep_loftee_for_csq_high: If False, include "OS" (Other Splice)
        :return: The order of consequence terms.
        """
        csq_order = []
        lof_values = (
            ["HC", "LC", hl.missing(hl.tstr)]
            if keep_loftee_for_csq_high
            else ["HC", "OS", "LC", hl.missing(hl.tstr)]
        )
        for lof in lof_values:
            for no_lof_flag in [True, False]:
                for consequence in CSQ_ORDER:
                    csq_order.append((lof, no_lof_flag, consequence))
        return csq_order

    def _get_sort_key(x, csq_order):
        """
        Get sort key for transcript consequence.

        :param x: The field to be sorted.
        :param csq_order: The order of consequence terms.
        :return: The sort key.
        """
        return csq_order[
            (
                x.lof,
                hl.or_else(hl.is_missing(x.lof_flags), False),
                x.most_severe_consequence,
            )
        ]

    ht = ht.collect_by_key("tx_annotation")

    csq_order = []
    # TODO: Do we add PolyPhen to prioritize missense variants?
    if keep_loftee_for_csq_high:
        ht = ht.annotate(
            tx_annotation=ht.tx_annotation.map(keep_loftee_for_high_impact)
        )

    csq_order = hl.literal(
        {x: i for i, x in enumerate(_get_csq_order(keep_loftee_for_csq_high))}
    )

    ht = ht.select(
        **hl.sorted(
            ht.tx_annotation,
            key=lambda x: _get_sort_key(x, csq_order),
        )[0]
    )

    return ht


def get_max_pext_per_gene(ht: hl.Table) -> hl.Table:
    """
    Get the maximum pext value of each gene.

    .. note::
        "For a minority of genes, when RSEM assigns higher relative expression to
        non-coding transcripts, the sum of the value of coding transcripts can be
        much smaller than the gene expression value of all transcripts, resulting in
        low pext scores for all coding variants in the gene, and thus resulting in
        possible filtering of all variants for a given gene. In many cases this
        appears to be the result of spurious non-coding transcripts with a high
        degree of exon overlap with true coding transcripts.

        To get around this artifact from affecting the analyses, you can calculate
        the maximum pext score for all variants across all protein coding genes,
        and removed any gene where the maximum pext score is below a given threshold."
        -- adapted from https://doi.org/10.1038/s41586-020-2329-2

    :param ht: Table of variants annotated with transcript expression information,
        and aggregated by additional grouping fields, including at least 'gene_id',
        'gene_symbol', 'most_severe_consequence', etc.
    :return: Table of maximum transcript expression proportion per gene and maximum
        mean transcript expression proportion across all tissues for each gene.
    """
    # TODO: Use isinstance to check if tx_annotation is already a struct.
    # TODO: Should we make it compatible with existing v2 pext?
    tissues = hl.eval(ht.tissues)
    csq_all = hl.literal(set(CSQ_CODING_ALL_IMPACT))

    # Filter to coding variants.
    ht = ht.filter(csq_all.contains(ht.most_severe_consequence))

    # Select expression proportion for each tissue and mean expression proportion.
    ht = ht.select(
        "gene_id",
        "gene_symbol",
        "exp_prop_mean",
        **{t: ht[t].expression_proportion for t in tissues},
    )

    # Group by gene and get max of expression proportion for each tissue and mean
    # expression proportion.
    max_ht = ht.group_by("gene_id", "gene_symbol").aggregate(
        max_pexts=hl.struct(
            **{t: hl.agg.max(ht[t]) for t in tissues},
            exp_prop_mean=hl.agg.max(ht.exp_prop_mean),
        )
    )

    return max_ht
