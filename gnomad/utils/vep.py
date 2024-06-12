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

VEP_VERSIONS = ["101", "105"]
CURRENT_VEP_VERSION = VEP_VERSIONS[-1]
"""
Versions of VEP used in gnomAD data, the latest version is 105.
"""

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

CSQ_SPLICE = [
    "splice_acceptor_variant",
    "splice_donor_variant",
    "splice_region_variant",
]
"""
Constant containing all splice consequences.
"""

POSSIBLE_REFS = ("GRCh37", "GRCh38")
"""
Constant containing supported references
"""

VEP_CONFIG_PATH = "file:///vep_data/vep-gcloud.json"
"""
Constant that contains the local path to the VEP config file
"""

VEP_CSQ_FIELDS = {
    "101": "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info",
    "105": "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|UNIPROT_ISOFORM|SOURCE|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|LoF|LoF_filter|LoF_flags|LoF_info",
}
"""
Constant that defines the order of VEP annotations used in VCF export, currently stored in a dictionary with the VEP version as the key.
"""

VEP_CSQ_HEADER = (
    "Consequence annotations from Ensembl VEP. Format:"
    f" {VEP_CSQ_FIELDS[CURRENT_VEP_VERSION]}"
)
"""
Constant that contains description for VEP used in VCF export.
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


def get_vep_help(vep_config_path: Optional[str] = None):
    """
    Return the output of vep --help which includes the VEP version.

    .. warning::
        If no `vep_config_path` is supplied, this function will only work for Dataproc clusters
        created with `hailctl dataproc start --vep`. It assumes that the command is `/path/to/vep`.

    :param vep_config_path: Optional path to use as the VEP config file. If None, `VEP_CONFIG_URI` environment variable is used
    :return: VEP help string
    """
    if vep_config_path is None:
        vep_config_path = os.environ["VEP_CONFIG_URI"]

    with hl.hadoop_open(vep_config_path) as vep_config_file:
        vep_config = json.load(vep_config_file)
        vep_command = vep_config["command"]
        vep_help = subprocess.check_output([vep_command[0]]).decode("utf-8")
        return vep_help


def get_vep_context(ref: Optional[str] = None) -> VersionedTableResource:
    """
    Get VEP context resource for the genome build `ref`.

    :param ref: Genome build. If None, `hl.default_reference` is used
    :return: VEPed context resource
    """
    import gnomad.resources.grch37.reference_data as grch37
    import gnomad.resources.grch38.reference_data as grch38

    if ref is None:
        ref = hl.default_reference().name

    if ref not in POSSIBLE_REFS:
        raise ValueError(
            f'get_vep_context passed {ref}. Expected one of {", ".join(POSSIBLE_REFS)}'
        )

    vep_context = grch37.vep_context if ref == "GRCh37" else grch38.vep_context
    return vep_context


def vep_or_lookup_vep(
    ht, reference_vep_ht=None, reference=None, vep_config_path=None, vep_version=None
):
    """
    VEP a table, or lookup variants in a reference database.

    .. warning::
        If `reference_vep_ht` is supplied, no check is performed to confirm `reference_vep_ht` was
        generated with the same version of VEP / VEP configuration as the VEP referenced in `vep_config_path`.

    :param ht: Input Table
    :param reference_vep_ht: A reference database with VEP annotations (must be in top-level `vep`)
    :param reference: If reference_vep_ht is not specified, find a suitable one in reference (if None, grabs from hl.default_reference)
    :param vep_config_path: vep_config to pass to hl.vep (if None, a suitable one for `reference` is chosen)
    :param vep_version: Version of VEPed context Table to use (if None, the default `vep_context` resource will be used)
    :return: VEPed Table
    """
    if reference is None:
        reference = hl.default_reference().name

    if vep_config_path is None:
        vep_config_path = VEP_CONFIG_PATH

    vep_help = get_vep_help(vep_config_path)

    with hl.hadoop_open(vep_config_path) as vep_config_file:
        vep_config = vep_config_file.read()

    if reference_vep_ht is None:
        if reference not in POSSIBLE_REFS:
            raise ValueError(
                f"vep_or_lookup_vep got {reference}. Expected one of"
                f" {', '.join(POSSIBLE_REFS)}"
            )

        vep_context = get_vep_context(reference)
        if vep_version is None:
            vep_version = vep_context.default_version

        if vep_version not in vep_context.versions:
            logger.warning(
                "No VEPed context Table available for genome build %s and VEP"
                " version %s, all variants will be VEPed using the following"
                " VEP:\n%s",
                reference,
                vep_version,
                vep_help,
            )
            return hl.vep(ht, vep_config_path)

        logger.info(
            "Using VEPed context Table from genome build %s and VEP version %s",
            reference,
            vep_version,
        )

        reference_vep_ht = vep_context.versions[vep_version].ht()
        vep_context_help = hl.eval(reference_vep_ht.vep_help)
        vep_context_config = hl.eval(reference_vep_ht.vep_config)

        assert vep_help == vep_context_help, (
            "The VEP context HT version does not match the version referenced in the"
            f" VEP config file.\nVEP context:\n{vep_context_help}\n\n VEP"
            f" config:\n{vep_help}"
        )

        assert vep_config == vep_context_config, (
            "The VEP context HT configuration does not match the configuration in"
            f" {vep_config_path}.\nVEP context:\n{vep_context_config}\n\n Current"
            f" config:\n{vep_config}"
        )

    ht = ht.annotate(vep=reference_vep_ht[ht.key].vep)

    vep_ht = ht.filter(hl.is_defined(ht.vep))
    revep_ht = ht.filter(hl.is_missing(ht.vep))
    revep_ht = hl.vep(revep_ht, vep_config_path)
    if "vep_proc_id" in list(revep_ht.row):
        revep_ht = revep_ht.drop("vep_proc_id")
    if "vep_proc_id" in list(vep_ht.row):
        vep_ht = vep_ht.drop("vep_proc_id")

    vep_ht = vep_ht.annotate_globals(
        vep_version=f"v{vep_version}", vep_help=vep_help, vep_config=vep_config
    )

    return vep_ht.union(revep_ht)


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


def filter_vep_to_canonical_transcripts(
    mt: Union[hl.MatrixTable, hl.Table],
    vep_root: str = "vep",
    filter_empty_csq: bool = False,
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filter VEP transcript consequences to those in the canonical transcript.

    :param mt: Input Table or MatrixTable.
    :param vep_root: Name used for VEP annotation. Default is 'vep'.
    :param filter_empty_csq: Whether to filter out rows where 'transcript_consequences' is empty. Default is False.
    :return: Table or MatrixTable with VEP transcript consequences filtered.
    """
    return filter_vep_transcript_csqs(
        mt, vep_root, synonymous=False, filter_empty_csq=filter_empty_csq
    )


def filter_vep_to_mane_select_transcripts(
    mt: Union[hl.MatrixTable, hl.Table],
    vep_root: str = "vep",
    filter_empty_csq: bool = False,
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filter VEP transcript consequences to those in the MANE Select transcript.

    :param mt: Input Table or MatrixTable.
    :param vep_root: Name used for VEP annotation. Default is 'vep'.
    :param filter_empty_csq: Whether to filter out rows where 'transcript_consequences' is empty. Default is False.
    :return: Table or MatrixTable with VEP transcript consequences filtered.
    """
    return filter_vep_transcript_csqs(
        mt,
        vep_root,
        synonymous=False,
        canonical=False,
        mane_select=True,
        filter_empty_csq=filter_empty_csq,
    )


def filter_vep_to_synonymous_variants(
    mt: Union[hl.MatrixTable, hl.Table],
    vep_root: str = "vep",
    filter_empty_csq: bool = False,
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filter VEP transcript consequences to those with a most severe consequence of 'synonymous_variant'.

    :param mt: Input Table or MatrixTable.
    :param vep_root: Name used for VEP annotation. Default is 'vep'.
    :param filter_empty_csq: Whether to filter out rows where 'transcript_consequences' is empty. Default is False.
    :return: Table or MatrixTable with VEP transcript consequences filtered.
    """
    return filter_vep_transcript_csqs(
        mt, vep_root, canonical=False, filter_empty_csq=filter_empty_csq
    )


def filter_vep_to_gene_list(
    t: Union[hl.MatrixTable, hl.Table],
    genes: List[str],
    match_by_gene_symbol: bool = False,
    vep_root: str = "vep",
    filter_empty_csq: bool = False,
):
    """
    Filter VEP transcript consequences to those in a set of genes.

    .. note::

       Filtering to a list of genes by their 'gene_id' or 'gene_symbol' will filter to
       all variants that are annotated to the gene, including
       ['upstream_gene_variant', 'downstream_gene_variant'], which will not be the
       same as if you filter to a gene interval. If you only want variants inside
       certain gene boundaries and a faster filter, you can first filter `t` to an
       interval list and then apply this filter.

    :param t: Input Table or MatrixTable.
    :param genes: Genes of interest to filter VEP transcript consequences to.
    :param match_by_gene_symbol: Whether to match values in `genes` to VEP transcript
        consequences by 'gene_symbol' instead of 'gene_id'. Default is False.
    :param vep_root: Name used for VEP annotation. Default is 'vep'.
    :param filter_empty_csq: Whether to filter out rows where 'transcript_consequences'
        is empty. Default is False.
    :return: Table or MatrixTable with VEP transcript consequences filtered.
    """
    return filter_vep_transcript_csqs(
        t,
        vep_root,
        synonymous=False,
        canonical=False,
        filter_empty_csq=filter_empty_csq,
        genes=genes,
        match_by_gene_symbol=match_by_gene_symbol,
    )


def vep_struct_to_csq(
    vep_expr: hl.expr.StructExpression,
    csq_fields: str = VEP_CSQ_FIELDS[CURRENT_VEP_VERSION],
    has_polyphen_sift: bool = True,
) -> hl.expr.ArrayExpression:
    """
    Given a VEP Struct, returns and array of VEP VCF CSQ strings (one per consequence in the struct).

    The fields and their order will correspond to those passed in `csq_fields`, which corresponds to the
    VCF header that is required to interpret the VCF CSQ INFO field.

    Note that the order is flexible and that all fields that are in the default value are supported.
    These fields will be formatted in the same way that their VEP CSQ counterparts are.

    While other fields can be added if their name are the same as those in the struct. Their value will be the result of calling
    hl.str(), so it may differ from their usual VEP CSQ representation.

    :param vep_expr: The input VEP Struct
    :param csq_fields: The | delimited list of fields to include in the CSQ (in that order), default is the CSQ fields of the CURRENT_VEP_VERSION.
    :param has_polyphen_sift: Whether the input VEP Struct has PolyPhen and SIFT annotations. Default is True.
    :return: The corresponding CSQ strings
    """
    _csq_fields = [f.lower() for f in csq_fields.split("|")]

    def get_csq_from_struct(
        element: hl.expr.StructExpression, feature_type: str
    ) -> hl.expr.StringExpression:
        # Most fields are 1-1, just lowercase
        fields = dict(element)

        # Add general exceptions
        fields.update(
            {
                "allele": element.variant_allele,
                "consequence": hl.delimit(element.consequence_terms, delimiter="&"),
                "feature_type": feature_type,
                "feature": (
                    element.transcript_id
                    if "transcript_id" in element
                    else (
                        element.regulatory_feature_id
                        if "regulatory_feature_id" in element
                        else (
                            element.motif_feature_id
                            if "motif_feature_id" in element
                            else ""
                        )
                    )
                ),
                "variant_class": vep_expr.variant_class,
            }
        )

        # Add exception for transcripts
        if feature_type == "Transcript":
            transcript_dict = {
                "canonical": hl.if_else(element.canonical == 1, "YES", ""),
                "ensp": element.protein_id,
                "gene": element.gene_id,
                "symbol": element.gene_symbol,
                "symbol_source": element.gene_symbol_source,
                "cdna_position": hl.str(element.cdna_start)
                + hl.if_else(
                    element.cdna_start == element.cdna_end,
                    "",
                    "-" + hl.str(element.cdna_end),
                ),
                "cds_position": hl.str(element.cds_start)
                + hl.if_else(
                    element.cds_start == element.cds_end,
                    "",
                    "-" + hl.str(element.cds_end),
                ),
                "mirna": hl.delimit(element.mirna, "&") if "mirna" in element else None,
                "protein_position": hl.str(element.protein_start)
                + hl.if_else(
                    element.protein_start == element.protein_end,
                    "",
                    "-" + hl.str(element.protein_end),
                ),
                "uniprot_isoform": (
                    hl.delimit(element.uniprot_isoform, "&")
                    if "uniprot_isoform" in element
                    else None
                ),
            }
            # Retain transcript dict updates only for fields that exist in the csq
            # fields.
            transcript_dict = {
                k: v
                for k, v in transcript_dict.items()
                if k in [x.lower() for x in csq_fields.split("|")]
            }
            fields.update(transcript_dict)

            if has_polyphen_sift:
                fields.update(
                    {
                        "sift": (
                            element.sift_prediction
                            + "("
                            + hl.format("%.3f", element.sift_score)
                            + ")"
                        ),
                        "polyphen": (
                            element.polyphen_prediction
                            + "("
                            + hl.format("%.3f", element.polyphen_score)
                            + ")"
                        ),
                    }
                )
            fields.update(
                {
                    "domains": hl.delimit(
                        element.domains.map(lambda d: d.db + ":" + d.name), "&"
                    ),
                }
            )
        elif feature_type == "MotifFeature":
            fields["motif_score_change"] = hl.format("%.3f", element.motif_score_change)
            if "transcription_factors" in element:
                fields["transcription_factors"] = hl.delimit(
                    element.transcription_factors, "&"
                )

        return hl.delimit(
            [hl.or_else(hl.str(fields.get(f, "")), "") for f in _csq_fields], "|"
        )

    csq = hl.empty_array(hl.tstr)
    for feature_field, feature_type in [
        ("transcript_consequences", "Transcript"),
        ("regulatory_feature_consequences", "RegulatoryFeature"),
        ("motif_feature_consequences", "MotifFeature"),
        ("intergenic_consequences", "Intergenic"),
    ]:
        csq = csq.extend(
            hl.or_else(
                vep_expr[feature_field].map(
                    lambda x: get_csq_from_struct(x, feature_type=feature_type)
                ),
                hl.empty_array(hl.tstr),
            )
        )

    return hl.or_missing(hl.len(csq) > 0, csq)


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


def explode_by_vep_annotation(
    t: Union[hl.Table, hl.MatrixTable],
    vep_annotation: str = "transcript_consequences",
    vep_root: str = "vep",
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Explode the specified VEP annotation on the input Table/MatrixTable.

    :param t: Input Table or MatrixTable.
    :param vep_annotation: Name of annotation in `vep_root` to explode.
    :param vep_root: Name used for root VEP annotation. Default is 'vep'.
    :return: Table or MatrixTable with exploded VEP annotation.
    """
    if vep_annotation not in t[vep_root].keys():
        raise ValueError(
            f"{vep_annotation} is not a row field of the {vep_root} annotation in"
            " Table/MatrixTable!"
        )
    # Create top-level annotation for `vep_annotation` and explode it.
    if isinstance(t, hl.Table):
        t = t.transmute(**{vep_annotation: t[vep_root][vep_annotation]})
        t = t.explode(t[vep_annotation])
    else:
        t = t.transmute_rows(**{vep_annotation: t[vep_root][vep_annotation]})
        t = t.explode_rows(t[vep_annotation])

    return t
