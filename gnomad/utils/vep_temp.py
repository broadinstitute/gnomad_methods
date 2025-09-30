# noqa: D100

import json
import logging
import os
import subprocess
from typing import Callable, Dict, List, Optional, Tuple, Union

import hail as hl
from deprecated import deprecated

from gnomad.resources.resource_utils import VersionedTableResource
from gnomad.utils.filtering import combine_functions

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Note that this is the current as of v105 with some included for backwards
# compatibility (VEP <= 75)
# See: https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
CSQ_CODING_HIGH_IMPACT = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",  # Considered high impact in v105, previously medium
]

CSQ_CODING_MEDIUM_IMPACT = [
    "initiator_codon_variant",  # deprecated
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",  # new in v79
]

CSQ_CODING_LOW_IMPACT = [
    "splice_region_variant",  # Considered low impact in v105, previously medium
    "incomplete_terminal_codon_variant",
    "start_retained_variant",  # new in v92
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
]

CSQ_NON_CODING = [
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "non_coding_exon_variant",  # deprecated
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "nc_transcript_variant",  # deprecated
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "feature_elongation",
    "regulatory_region_variant",
    "feature_truncation",
    "intergenic_variant",
]

CSQ_ORDER = (
    CSQ_CODING_HIGH_IMPACT
    + CSQ_CODING_MEDIUM_IMPACT
    + CSQ_CODING_LOW_IMPACT
    + CSQ_NON_CODING
)

CSQ_CODING = CSQ_CODING_HIGH_IMPACT + CSQ_CODING_MEDIUM_IMPACT + CSQ_CODING_LOW_IMPACT
"""
Constant containing all coding consequences.
"""

LOFTEE_LABELS = ["HC", "LC", "OS"]
"""
Constant that contains annotations added by LOFTEE.
"""

LOF_CSQ_SET = {
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
}
"""
Set containing loss-of-function consequence strings.
"""

POLYPHEN_ORDER = ["probably_damaging", "possibly_damaging", "benign"]
"""
Order of PolyPhen predictions from most to least severe.
"""


def filter_vep_transcript_csqs(
    t: Union[hl.Table, hl.MatrixTable],
    vep_root: str = "vep",
    synonymous: bool = True,
    canonical: bool = True,
    mane_select: bool = False,
    filter_empty_csq: bool = True,
    ensembl_only: bool = True,
    protein_coding: bool = False,
    csqs: List[str] = None,
    keep_csqs: bool = True,
    genes: Optional[List[str]] = None,
    keep_genes: bool = True,
    match_by_gene_symbol: bool = False,
    additional_filtering_criteria: Optional[List[Callable]] = None,
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Filter VEP transcript consequences based on specified criteria, and optionally filter to variants where transcript consequences is not empty after filtering.

    Transcript consequences can be filtered to those where 'most_severe_consequence' is
    'synonymous_variant' and/or the transcript is the canonical transcript, if the
    `synonymous` and `canonical` parameter are set to True, respectively.

    If `filter_empty_csq` parameter is set to True, the Table/MatrixTable is filtered
    to variants where 'transcript_consequences' within the VEP annotation is not empty
    after the specified filtering criteria is applied.

    :param t: Input Table or MatrixTable.
    :param vep_root: Name used for VEP annotation. Default is 'vep'.
    :param synonymous: Whether to filter to variants where the most severe consequence
        is 'synonymous_variant'. Default is True.
    :param canonical: Whether to filter to only canonical transcripts. Default is True.
    :param mane_select: Whether to filter to only MANE Select transcripts. Default is
        False.
    :param filter_empty_csq: Whether to filter out rows where 'transcript_consequences'
        is empty, after filtering 'transcript_consequences' to the specified criteria.
        Default is True.
    :param ensembl_only: Whether to filter to only Ensembl transcripts. This option is
        useful for deduplicating transcripts that are the same between RefSeq and
        Emsembl. Default is True.
    :param protein_coding: Whether to filter to only protein-coding transcripts.
        Default is False.
    :param csqs: Optional list of consequence terms to filter to. Transcript
        consequences are filtered to those where 'most_severe_consequence' is in the
        list of consequence terms `csqs`. Default is None.
    :param keep_csqs: Whether to keep transcript consequences that are in `csqs`. If
        set to False, transcript consequences that are in `csqs` will be removed.
        Default is True.
    :param genes: Optional list of genes to filter VEP transcript consequences to.
        Default is None.
    :param keep_genes: Whether to keep transcript consequences that are in `genes`. If
        set to False, transcript consequences that are in `genes` will be removed.
        Default is True.
    :param match_by_gene_symbol: Whether to match values in `genes` to VEP transcript
        consequences by 'gene_symbol' instead of 'gene_id'. Default is False.
    :param additional_filtering_criteria: Optional list of additional filtering
        criteria to apply to the VEP transcript consequences.
    :return: Table or MatrixTable filtered to specified criteria.
    """
    if not synonymous and not (canonical or mane_select) and not filter_empty_csq:
        logger.warning("No changes have been made to input Table/MatrixTable!")
        return t

    transcript_csqs = t[vep_root].transcript_consequences
    criteria = [lambda csq: True]
    if synonymous:
        logger.info("Filtering to most severe consequence of synonymous_variant...")
        csqs = ["synonymous_variant"]
    if csqs is not None:
        csqs = hl.literal(csqs)
        if keep_csqs:
            criteria.append(lambda csq: csqs.contains(csq.most_severe_consequence))
        else:
            criteria.append(lambda csq: ~csqs.contains(csq.most_severe_consequence))
    if canonical:
        logger.info("Filtering to canonical transcripts")
        criteria.append(lambda csq: csq.canonical == 1)
    if mane_select:
        logger.info("Filtering to MANE Select transcripts...")
        criteria.append(lambda csq: hl.is_defined(csq.mane_select))
    if ensembl_only:
        logger.info("Filtering to Ensembl transcripts...")
        criteria.append(lambda csq: csq.transcript_id.startswith("ENST"))
    if protein_coding:
        logger.info("Filtering to protein coding transcripts...")
        criteria.append(lambda csq: csq.biotype == "protein_coding")
    if genes is not None:
        logger.info("Filtering to genes of interest...")
        genes = hl.literal(genes)
        gene_field = "gene_symbol" if match_by_gene_symbol else "gene_id"
        if keep_genes:
            criteria.append(lambda csq: genes.contains(csq[gene_field]))
        else:
            criteria.append(lambda csq: ~genes.contains(csq[gene_field]))
    if additional_filtering_criteria is not None:
        logger.info("Filtering to variants with additional criteria...")
        criteria = criteria + additional_filtering_criteria

    transcript_csqs = transcript_csqs.filter(lambda x: combine_functions(criteria, x))
    is_mt = isinstance(t, hl.MatrixTable)
    vep_data = {vep_root: t[vep_root].annotate(transcript_consequences=transcript_csqs)}
    t = t.annotate_rows(**vep_data) if is_mt else t.annotate(**vep_data)

    if filter_empty_csq:
        transcript_csq_expr = t[vep_root].transcript_consequences
        filter_expr = hl.is_defined(transcript_csq_expr) & (
            hl.len(transcript_csq_expr) > 0
        )
        t = t.filter_rows(filter_expr) if is_mt else t.filter(filter_expr)

    return t


########################################################################################
# Functions for handling most severe consequences.
########################################################################################


########################################################################################
# Functions for determining and annotating most severe consequences.
#
# get_most_severe_consequence_expr finds the most severe consequence from an array of
# consequences.
#
# add_most_severe_consequence_to_consequence uses get_most_severe_consequence_expr to
# add a most_severe_consequence field to a transcript consequence or array of transcript
# consequences.
#
# add_most_severe_csq_to_tc_within_vep_root uses
# add_most_severe_consequence_to_consequence to add a most_severe_consequence field to
# a transcript consequence or array of transcript consequences within a VEP root
# annotation.
########################################################################################
def get_most_severe_consequence_expr(
    csq_expr: hl.expr.ArrayExpression,
    csq_order: Optional[List[str]] = None,
) -> hl.expr.StringExpression:
    """
    Get the most severe consequence from a collection of consequences.

    This is for a given transcript, as there are often multiple annotations for a single
    transcript: e.g. splice_region_variant&intron_variant -> splice_region_variant

    :param csq_expr: ArrayExpression of consequences.
    :param csq_order: Optional list indicating the order of VEP consequences, sorted
        from high to low impact. Default is None, which uses the value of the
        `CSQ_ORDER` global.
    :return: Most severe consequence in `csq_expr`.
    """
    if csq_order is None:
        csq_order = CSQ_ORDER
    csqs = hl.literal(csq_order)

    return csqs.find(lambda c: csq_expr.contains(c))


def add_most_severe_consequence_to_consequence(
    tc: Union[hl.expr.StructExpression, hl.expr.ArrayExpression],
    csq_order: Optional[List[str]] = None,
    most_severe_csq_field: str = "most_severe_consequence",
) -> Union[hl.expr.StructExpression, hl.expr.ArrayExpression]:
    """
    Add a `most_severe_consequence` field to a transcript consequence or array of transcript consequences.

    For a single transcript consequence, `tc` should be a StructExpression with a
    `consequence_terms` field, e.g. Struct(consequence_terms=['missense_variant']).
    For an array of transcript consequences, `tc` should be an ArrayExpression of
    StructExpressions with a `consequence_terms` field.

    :param tc: Transcript consequence or array of transcript consequences to annotate.
    :param csq_order: Optional list indicating the order of VEP consequences, sorted
        from high to low impact. Default is None, which uses the value of the
        `CSQ_ORDER` global.
    :param most_severe_csq_field: Field name to use for most severe consequence. Default
        is 'most_severe_consequence'.
    :return: Transcript consequence or array of transcript consequences annotated with
        the most severe consequence.
    """
    csq = lambda x: get_most_severe_consequence_expr(x.consequence_terms, csq_order)
    if isinstance(tc, hl.expr.StructExpression):
        return tc.annotate(**{most_severe_csq_field: csq(tc)})
    else:
        return tc.map(lambda x: x.annotate(**{most_severe_csq_field: csq(x)}))


def add_most_severe_csq_to_tc_within_vep_root(
    t: Union[hl.Table, hl.MatrixTable],
    vep_root: str = "vep",
    csq_field: str = "transcript_consequences",
    most_severe_csq_field: str = "most_severe_consequence",
    csq_order: Optional[List[str]] = None,
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Add `most_severe_csq_field` annotation to `csq_field` within the `vep_root` annotation.

    :param t: Input Table or MatrixTable.
    :param vep_root: Root for vep annotation (probably vep).
    :param csq_field: Field name of VEP consequences ArrayExpression within `vep_root`
        to add most severe consequence to. Default is 'transcript_consequences'.
    :param most_severe_csq_field: Field name to use for most severe consequence. Default
        is 'most_severe_consequence'.
    :param csq_order: Optional list indicating the order of VEP consequences, sorted
        from high to low impact. Default is None, which uses the value of the
        `CSQ_ORDER` global.
    :return: Table or MatrixTable with most_severe_consequence annotation added.
    """
    vep_expr = t[vep_root]
    csq_expr = vep_expr[csq_field]
    vep_expr = vep_expr.annotate(
        **{
            csq_field: add_most_severe_consequence_to_consequence(
                csq_expr,
                csq_order=csq_order,
                most_severe_csq_field=most_severe_csq_field,
            )
        }
    )
    return (
        t.annotate_rows(**{vep_root: vep_expr})
        if isinstance(t, hl.MatrixTable)
        else t.annotate(**{vep_root: vep_expr})
    )


def filter_to_most_severe_consequences(
    csq_expr: hl.expr.ArrayExpression,
    csq_order: List[str] = CSQ_ORDER,
    most_severe_csq_field: str = "most_severe_consequence",
) -> hl.expr.ArrayExpression:
    """
    Filter an array of VEP consequences to all entries that have the most severe consequence.

    This function expects that all entries in the `csq_list` are already annotated with
    the most severe consequence using `add_most_severe_consequence_to_consequence` or
    `add_most_severe_csq_to_tc_within_vep_root`.

    .. note::

        - If you have multiple lists of consequences and want to determine the most
          severe consequence across all lists, consider using
          `get_most_severe_consequence_for_summary`.

        - If you want to group consequences by gene and determine the most severe
          consequence for each gene, consider using `process_consequences`.

    :param csq_expr: ArrayExpression of VEP consequences.
    :param csq_order: List indicating the order of VEP consequences, sorted from high to
        low impact. Default is the value of the `CSQ_ORDER` global.
    :param most_severe_csq_field: Field containing the most severe consequence for each
        consequence in `csq_list`. Default is 'most_severe_consequence'.
    :return: ArrayExpression with of the consequences that match the most severe
        consequence.
    """
    # Get the highest impact csq label.
    ms_csq = get_most_severe_consequence_expr(
        csq_expr.map(lambda x: x[most_severe_csq_field]), csq_order=csq_order
    )

    # Filter to only consequences with the highest impact csq label, and return missing
    # if the most severe consequence is missing.
    return hl.or_missing(
        hl.is_defined(ms_csq),
        csq_expr.filter(lambda x: x[most_severe_csq_field] == ms_csq),
    )


# TODO: Can add this to `filter_vep_transcript_csqs`, but need to also make that
#  function work with just an array of transcript consequences.
def filter_vep_consequences_by_loftee(
    csq_expr: hl.expr.ArrayExpression,
    loftee_labels: Optional[List[str]] = None,
    no_lof_flags: bool = False,
    keep: bool = True,
) -> hl.expr.StructExpression:
    """
    Filter VEP transcript consequences by LOFTEE.

    :param csq_expr: ArrayExpression of VEP consequences with LOFTEE annotations.
    :param loftee_labels: List of LOFTEE labels to filter to. Default is None, which
        filters to all LOFTEE labels.
    :param no_lof_flags: Whether to filter to consequences with no LOFTEE flags.
        Default is False.
    :param keep: Whether to keep the consequences that match the filter criteria.
        Default is True.
    :return: StructExpression with the filtered consequences.
    """
    filter_criteria = [lambda csq: True]

    if loftee_labels:
        logger.info("Filtering to LOFTEE labels: %s...", loftee_labels)
        filter_criteria.append(lambda x: hl.set(loftee_labels).contains(x.lof))

    if no_lof_flags:
        logger.info("Filtering to consequences with no LOFTEE flags...")
        filter_criteria.append(
            lambda x: hl.is_missing(x.lof_flags) | (x.lof_flags == "")
        )

    return csq_expr.filter(lambda x: combine_functions(filter_criteria, x), keep=keep)


def _prioritize_loftee_hc_no_flags(
    most_severe_csq: Union[hl.Table, hl.expr.StructExpression],
) -> hl.expr.StructExpression:
    """
    Prioritize LOFTEE HC LOF consequences with no LOF flags.

    Given the result of `get_most_severe_csq_from_multiple_csq_lists`, this function
    will filter the transcript consequences to only include those with no LOF flags if
    the most severe consequence is a LOFTEE HC LOF and there are transcript consequences
    with no LOF flags.

    :param most_severe_csq: Table or StructExpression with most severe consequence
        information. This should be the result of
        `get_most_severe_csq_from_multiple_csq_lists`.
    :return: StructExpression with HC LOF consequences with no LOF flags if they exist,
        otherwise all transcript consequences.
    """
    tcl = most_severe_csq.transcript_consequences

    # Filter transcript consequences to only consequences that have no LOF flags.
    no_flags = filter_vep_consequences_by_loftee(tcl, no_lof_flags=True)

    # If the most severe consequence is a LOFTEE HC LOF and there are transcript
    # consequences with no LOF flags, return only those transcripts.
    return hl.if_else(
        (most_severe_csq.lof == "HC") & (hl.len(no_flags) > 0),
        no_flags,
        tcl,
    )


def get_most_severe_csq_from_multiple_csq_lists(
    vep_expr: hl.expr.StructExpression,
    csq_order: List[str] = CSQ_ORDER,
    loftee_labels: List[str] = LOFTEE_LABELS,
    include_transcript_csqs: bool = False,
    prioritize_protein_coding: bool = True,
    csq_list_order: Union[List[str], Tuple[str]] = (
        "transcript_consequences",
        "regulatory_feature_consequences",
        "motif_feature_consequences",
        "intergenic_consequences",
    ),
    add_order_by_csq_list: Dict[str, List[str]] = None,
) -> hl.Table:
    """
    Process multiple VEP consequences lists to determine the most severe consequence.

    Useful for generating summary annotations for VEP consequences.

    Adds the following annotations:
        - most_severe_csq: Most severe consequence for variant.
        - protein_coding: Whether the variant is present on a protein-coding transcript.
        - lof: Whether the variant is a loss-of-function variant.
        - no_lof_flags: Whether the variant has any LOFTEE flags (True if no flags).

    If `include_transcript_csqs` is True, an additional annotation is added:
        - transcript_consequences: All transcript consequences for the most severe
          consequence.

    If `prioritize_protein_coding` is True and "transcript_consequences" is in
    `csq_list_order`, protein-coding transcripts are prioritized by filtering to only
    protein-coding transcripts and determining the most severe consequence. If no
    protein-coding transcripts are present, determine the most severe consequence for
    all transcripts. If additional VEP consequences lists are requested, process those
    lists in the order they appear in `csq_list_order`.

    If `add_order_by_csq_list` is provided, additional ordering is applied to the
    consequences in the list. The key is the name of the consequences list and the value
    is the order of consequences, sorted from high to low impact. An example use of this
    parameter is to prioritize PolyPhen consequences for protein-coding transcripts.

    .. note::

        Assumes input Table is annotated with VEP and that VEP annotations have been
        filtered to canonical transcripts if wanted.

    :param vep_expr: StructExpression of VEP consequences to get the most severe
        consequence from.
    :param csq_order: Order of VEP consequences, sorted from high to low impact.
        Default is CSQ_ORDER.
    :param loftee_labels: Annotations added by LOFTEE, sorted from high to low impact.
        Default is LOFTEE_LABELS.
    :param include_transcript_csqs: Whether to include all transcript consequences for
        the most severe consequence. Default is False.
    :param prioritize_protein_coding: Whether to prioritize protein-coding transcripts
        when determining the worst consequence. Default is True.
    :param csq_list_order: Order of VEP consequences lists to be processed. Default is
        ('transcript_consequences', 'regulatory_feature_consequences',
        'motif_feature_consequences', 'intergenic_consequences').
    :param add_order_by_csq_list: Dictionary of additional ordering for VEP consequences
        lists. The key is the name of the consequences list and the value is the order
        of consequences, sorted from high to low impact. Default is None.
    :return: Table annotated with VEP summary annotations.
    """
    if add_order_by_csq_list is None:
        add_order_by_csq_list = {}

    def _get_most_severe_csq(
        csq_list: hl.ArrayExpression,
        protein_coding: bool = False,
        include_csqs: bool = False,
        prioritize_loftee: bool = False,
        additional_order: Optional[List[str]] = None,
    ) -> hl.StructExpression:
        """
        Filter a list of consequences to those that have the most severe consequence.

        If `protein_coding` is True, filter to only protein-coding transcripts before
        determining the most severe consequence. If `prioritize_loftee` is True,
        prioritize LOFTEE consequences by filtering to only LOFTEE consequences and
        determining the most severe consequence. If `additional_order` is provided,
        additional ordering is applied to the consequences in the list.

        :param csq_list: ArrayExpression of VEP consequences.
        :param protein_coding: Whether to filter to only protein-coding transcripts
            before determining the most severe consequence. Default is False.
        :param include_csqs: Whether to include all transcript consequences for the most
            severe consequence. Default is False.
        :param prioritize_loftee: Whether to prioritize LOFTEE consequences. Default is
            False.
        :param additional_order: Tuple indicating the additional ordering to apply to
            the consequences in the list. The first element is the name of the
            consequences list and the second element is the order of consequences,
            sorted from high to low impact. Default is None.
        :return: StructExpression with the most severe consequence and the list of
            consequences that match the most severe consequence.
        """
        if protein_coding:
            csq_list = csq_list.filter(lambda x: x.biotype == "protein_coding")

        lof = hl.missing(hl.tstr)
        no_lof_flags = hl.missing(hl.tbool)
        if prioritize_loftee:
            lof_csq = filter_to_most_severe_consequences(csq_list, loftee_labels, "lof")
            lof = lof_csq[0].lof

            # Check if any of the lof consequences have no lof_flags.
            no_lof_flags = (
                hl.len(filter_vep_consequences_by_loftee(lof_csq, no_lof_flags=True))
                > 0
            )

            # If there are no lof consequences, set the consequence list to the original
            # list.
            csq_list = hl.coalesce(lof_csq, csq_list)

        # Add most_severe_consequence to each consequence.
        csq_list = add_most_severe_consequence_to_consequence(csq_list, csq_order)

        # Get the most severe consequence of all consequences in the list.
        csq_list = filter_to_most_severe_consequences(csq_list, csq_order=csq_order)
        ms_csq = csq_list[0].most_severe_consequence
        result = hl.struct(
            most_severe_consequence=ms_csq,
            protein_coding=protein_coding,
            lof=lof,
            no_lof_flags=no_lof_flags,
        )

        if additional_order is not None:
            # Get the highest impact consequences from the additional ordering.
            add_csq_expr = filter_to_most_severe_consequences(
                csq_list,
                most_severe_csq_field=additional_order[0],
                csq_order=additional_order[1],
            )
            # If there are consequences from the additional ordering, set the
            # consequence list to the additional ordering, otherwise keep the original
            # list.
            csq_list = hl.coalesce(add_csq_expr, csq_list)

        result = result.annotate(consequences=csq_list)
        csq_list = hl.or_missing(hl.len(csq_list) > 0, result)

        # Drop the consequences field if not requested.
        if not include_csqs:
            csq_list = csq_list.drop("consequences")

        return csq_list

    # Get type of transcript_consequences field for use with hl.missing for other
    # consequence lists.
    tc_dtype = None
    if include_transcript_csqs and "transcript_consequences" in csq_list_order:
        tc_dtype = (
            vep_expr["transcript_consequences"]
            .map(lambda x: x.annotate(most_severe_consequence=hl.missing(hl.tstr)))
            .dtype
        )

    # Get the most severe consequence for each VEP consequences list.
    ms_csqs_list = []
    for c in csq_list_order:
        if c not in vep_expr:
            logger.warning(f"VEP consequences list %s not found in input!", c)
            continue
        csqs = vep_expr[c]
        is_tc = c == "transcript_consequences"
        ms_csqs = _get_most_severe_csq(
            csqs,
            # Only include transcript consequences if requested and the current list is
            # for transcript consequences.
            prioritize_loftee=True if is_tc else False,
            include_csqs=include_transcript_csqs and is_tc,
            additional_order=add_order_by_csq_list.get(c),
        )

        # If prioritizing protein-coding transcripts, get the most severe consequence
        # for protein-coding transcripts and coalesce with the current most severe
        # transcript consequence.
        if is_tc and prioritize_protein_coding:
            ms_csqs = hl.coalesce(
                _get_most_severe_csq(
                    csqs,
                    protein_coding=True,
                    prioritize_loftee=True,
                    include_csqs=include_transcript_csqs,
                    additional_order=add_order_by_csq_list.get(c),
                ),
                ms_csqs,
            )

        # If the current list is not for transcript consequences, annotate with missing
        # for transcript_consequences if transcript consequences are requested.
        if tc_dtype is not None:
            tc_expr = ms_csqs.consequences if is_tc else hl.missing(tc_dtype)
            ms_csqs = ms_csqs.annotate(consequences=tc_expr)

        ms_csqs_list.append(ms_csqs)

    ms_csqs = hl.coalesce(*ms_csqs_list)

    # Rename most_severe_consequence to most_severe_csq for consistency with older
    # version of code.
    rename_map = {"most_severe_consequence": "most_severe_csq"}
    if tc_dtype is not None:
        rename_map["consequences"] = "transcript_consequences"

    return ms_csqs.rename(rename_map)


@deprecated(reason="Replaced by get_most_severe_csq_from_multiple_csq_lists")
def get_most_severe_consequence_for_summary(
    ht: hl.Table,
    csq_order: List[str] = CSQ_ORDER,
    loftee_labels: List[str] = LOFTEE_LABELS,
) -> hl.Table:
    """
    Use `get_most_severe_csq_from_multiple_csq_lists` instead, this function is deprecated.

    Prepare a hail Table for summary statistics generation.

    Adds the following annotations:
        - most_severe_csq: Most severe consequence for variant
        - protein_coding: Whether the variant is present on a protein-coding transcript
        - lof: Whether the variant is a loss-of-function variant
        - no_lof_flags: Whether the variant has any LOFTEE flags (True if no flags)

    Assumes input Table is annotated with VEP and that VEP annotations have been filtered to canonical transcripts.

    :param ht: Input Table.
    :param csq_order: Order of VEP consequences, sorted from high to low impact. Default is CSQ_ORDER.
    :param loftee_labels: Annotations added by LOFTEE. Default is LOFTEE_LABELS.
    :return: Table annotated with VEP summary annotations.
    """
    csq_expr = get_most_severe_csq_from_multiple_csq_lists(
        ht.vep, csq_order=csq_order, loftee_labels=loftee_labels
    )

    # Rename most_severe_consequence to most_severe_csq for consistency with older
    # version of code.
    csq_expr = csq_expr.rename({"most_severe_consequence": "most_severe_csq"})

    return ht.annotate(**csq_expr)


def process_consequences(
    t: Union[hl.MatrixTable, hl.Table],
    vep_root: str = "vep",
    penalize_flags: bool = True,
    csq_order: Optional[List[str]] = None,
    has_polyphen: bool = True,
    prioritize_protein_coding: bool = False,
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Add most_severe_consequence into [vep_root].transcript_consequences, and worst_csq_by_gene, any_lof into [vep_root].

    `most_severe_consequence` is the worst consequence for a transcript.

    .. note::

        From gnomAD v4.0 on, the PolyPhen annotation was removed from the VEP Struct
        in the release HTs. When using this function with gnomAD v4.0 or later,
        set `has_polyphen` to False.

    :param t: Input Table or MatrixTable.
    :param vep_root: Root for VEP annotation (probably "vep").
    :param penalize_flags: Whether to penalize LOFTEE flagged variants, or treat them
        as equal to HC.
    :param csq_order: Optional list indicating the order of VEP consequences, sorted
        from high to low impact. Default is None, which uses the value of the
        `CSQ_ORDER` global.
    :param has_polyphen: Whether the input VEP Struct has a PolyPhen annotation which
        will be used to modify the consequence score. Default is True.
    :param prioritize_protein_coding: Whether to prioritize protein-coding transcripts
        when determining the worst consequence. Default is False.
    :return: HT/MT with better formatted consequences.
    """
    # If has_polyphen is True, set the order of PolyPhen consequences.
    polyphen_order = POLYPHEN_ORDER if has_polyphen else None

    def _find_worst_transcript_consequence(
        tcl: hl.expr.ArrayExpression,
    ) -> hl.expr.StructExpression:
        """
        Find the worst transcript consequence in an array of transcript consequences.

        :param tcl: Array of transcript consequences.
        :return: Worst transcript consequence.
        """
        ms_csq = get_most_severe_csq_from_multiple_csq_lists(
            hl.struct(transcript_consequences=tcl),
            csq_order=csq_order,
            include_transcript_csqs=True,
            prioritize_protein_coding=prioritize_protein_coding,
            csq_list_order=["transcript_consequences"],
            add_order_by_csq_list={"transcript_consequences": polyphen_order},
        )
        tcl = ms_csq.transcript_consequences

        # Penalize LOFTEE flagged variants.
        tcl = _prioritize_loftee_hc_no_flags(ms_csq) if penalize_flags else tcl

        return hl.or_missing(hl.len(tcl) > 0, tcl[0])

    # Annotate each transcript consequence with the 'most_severe_consequence'.
    csqs = t[vep_root].transcript_consequences.map(
        lambda tc: add_most_severe_consequence_to_consequence(tc, csq_order)
    )

    # Group transcript consequences by gene and find the worst consequence for each.
    gene_csqs = (
        csqs.group_by(lambda tc: tc.gene_symbol)
        .map_values(_find_worst_transcript_consequence)
        .values()
    )

    # Filter transcript consequences to only include canonical transcripts.
    canonical = csqs.filter(lambda csq: csq.canonical == 1)
    gene_canonical = (
        canonical.group_by(lambda tc: tc.gene_symbol)
        .map_values(_find_worst_transcript_consequence)
        .values()
    )

    # Annotate the HT/MT with the worst consequence for each gene and variant.
    vep_data = t[vep_root].annotate(
        transcript_consequences=csqs,
        worst_consequence_term=get_most_severe_consequence_expr(
            csqs.map(lambda csq: csq.most_severe_consequence), csq_order
        ),
        worst_csq_by_gene=gene_csqs,
        worst_csq_for_variant=_find_worst_transcript_consequence(csqs),
        worst_csq_by_gene_canonical=gene_canonical,
        worst_csq_for_variant_canonical=_find_worst_transcript_consequence(canonical),
    )

    return (
        t.annotate_rows(**{vep_root: vep_data})
        if isinstance(t, hl.MatrixTable)
        else t.annotate(**{vep_root: vep_data})
    )
