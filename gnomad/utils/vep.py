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

# Note that these terms are current as of v105 with some included for backwards
# compatibility (VEP <= 75). The impact groupings are loosely based on VEP's categories
# but have been adjusted to better serve gnomAD's use cases.
# See: https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
CSQ_CODING_HIGH_IMPACT = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
]

CSQ_CODING_MEDIUM_IMPACT = [
    "start_lost",  # considered high impact in v105, previously medium
    "initiator_codon_variant",  # deprecated
    "transcript_amplification",  # considered high impact in v105, previously medium
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",  # new in v79
]

CSQ_CODING_LOW_IMPACT = [
    "splice_donor_5th_base_variant",  # new in v105
    "splice_region_variant",  # Considered low impact in v105, previously medium
    "splice_donor_region_variant",  # new in v105
    "splice_polypyrimidine_tract_variant",  # new in v105
    "incomplete_terminal_codon_variant",
    "start_retained_variant",  # new in v92
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "coding_transcript_variant",
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
    "sequence_variant",
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
    csq_order = csq_order or CSQ_ORDER
    return hl.literal(csq_order).find(lambda c: csq_expr.contains(c))


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
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Add most_severe_consequence into [vep_root].transcript_consequences, and worst_csq_by_gene, any_lof into [vep_root].

    `most_severe_consequence` is the worst consequence for a transcript.

    Each transcript consequence is annotated with a `csq_score` which is a combination
    of the index of the consequence's `most_severe_consequence` in `csq_order` and an
    extra deduction for loss-of-function consequences, and polyphen predictions if
    `has_polyphen` is True. Lower scores translate to higher severity.

    The score adjustment is as follows:
        - lof == 'HC' & NO lof_flags (-1000 if penalize_flags, -500 if not)
        - lof == 'HC' & lof_flags (-500)
        - lof == 'OS' (-20)
        - lof == 'LC' (-10)
        - everything else (0)

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
    :return: MT with better formatted consequences.
    """
    if csq_order is None:
        csq_order = CSQ_ORDER
    csqs = hl.literal(csq_order)

    # Assign a score to each consequence based on the order in csq_order.
    csq_dict = hl.literal(dict(zip(csq_order, range(len(csq_order)))))

    def _find_worst_transcript_consequence(
        tcl: hl.expr.ArrayExpression,
    ) -> hl.expr.StructExpression:
        """
        Find the worst transcript consequence in an array of transcript consequences.

        :param tcl: Array of transcript consequences.
        :return: Worst transcript consequence.
        """
        flag = 500
        no_flag = flag * (1 + penalize_flags)

        # Score each consequence based on the order in csq_order.
        score_expr = tcl.map(
            lambda tc: csq_dict[csqs.find(lambda x: x == tc.most_severe_consequence)]
        )

        # Determine the score adjustment based on the consequence's LOF and LOF flags.
        sub_expr = tcl.map(
            lambda tc: (
                hl.case(missing_false=True)
                .when((tc.lof == "HC") & hl.or_else(tc.lof_flags == "", True), no_flag)
                .when((tc.lof == "HC") & (tc.lof_flags != ""), flag)
                .when(tc.lof == "OS", 20)
                .when(tc.lof == "LC", 10)
                .default(0)
            )
        )

        # If requested, determine the score adjustment based on the consequence's
        # PolyPhen prediction.
        if has_polyphen:
            polyphen_sub_expr = tcl.map(
                lambda tc: (
                    hl.case(missing_false=True)
                    .when(tc.polyphen_prediction == "probably_damaging", 0.5)
                    .when(tc.polyphen_prediction == "possibly_damaging", 0.25)
                    .when(tc.polyphen_prediction == "benign", 0.1)
                    .default(0)
                )
            )
            sub_expr = hl.map(lambda s, ps: s + ps, sub_expr, polyphen_sub_expr)

        # Calculate the final consequence score.
        tcl = hl.map(
            lambda tc, s, ss: tc.annotate(csq_score=s - ss), tcl, score_expr, sub_expr
        )

        # Return the worst consequence based on the calculated score.
        return hl.or_missing(hl.len(tcl) > 0, hl.sorted(tcl, lambda x: x.csq_score)[0])

    # Annotate each transcript consequence with the 'most_severe_consequence'.
    transcript_csqs = add_most_severe_consequence_to_consequence(
        t[vep_root].transcript_consequences, csq_order
    )

    # Group transcript consequences by gene and find the worst consequence for each.
    gene_dict = transcript_csqs.group_by(lambda tc: tc.gene_symbol)
    worst_csq_gene = gene_dict.map_values(_find_worst_transcript_consequence).values()
    sorted_scores = hl.sorted(worst_csq_gene, key=lambda tc: tc.csq_score)

    # Filter transcript consequences to only include canonical transcripts and find the
    # worst consequence for each gene.
    canonical = transcript_csqs.filter(lambda csq: csq.canonical == 1)
    gene_canonical_dict = canonical.group_by(lambda tc: tc.gene_symbol)
    worst_csq_gene_canonical = gene_canonical_dict.map_values(
        _find_worst_transcript_consequence
    ).values()
    sorted_canonical_scores = hl.sorted(
        worst_csq_gene_canonical, key=lambda tc: tc.csq_score
    )

    # Annotate the HT/MT with the worst consequence for each gene and variant.
    vep_data = t[vep_root].annotate(
        transcript_consequences=transcript_csqs,
        worst_consequence_term=csqs.find(
            lambda c: transcript_csqs.map(
                lambda csq: csq.most_severe_consequence
            ).contains(c)
        ),
        worst_csq_by_gene=sorted_scores,
        worst_csq_for_variant=hl.or_missing(
            hl.len(sorted_scores) > 0, sorted_scores[0]
        ),
        worst_csq_by_gene_canonical=sorted_canonical_scores,
        worst_csq_for_variant_canonical=hl.or_missing(
            hl.len(sorted_canonical_scores) > 0, sorted_canonical_scores[0]
        ),
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
    csq_order: Optional[List[str]] = None,
    loftee_labels: Optional[List[str]] = None,
    prioritize_protein_coding: bool = False,
    prioritize_loftee: bool = False,
    prioritize_loftee_no_flags: bool = False,
    additional_order_field: Optional[str] = None,
    additional_order: Optional[List[str]] = None,
) -> hl.StructExpression:
    """
    Filter an array of VEP consequences to all entries that have the most severe consequence.

    Returns a struct with the following annotations:

        - most_severe_consequence: Most severe consequence for variant.
        - lof: Whether the variant is a loss-of-function variant.
        - no_lof_flags: Whether the variant has any LOFTEE flags (True if no flags).
        - consequences: Array of consequences that match the most severe consequence.

    .. note::

        - If you have multiple lists of consequences (such as lists of both
          'transcript_consequences'  and 'intergenic_consequences') and want to
          determine the most severe consequence across all lists, consider using
          `get_most_severe_csq_from_multiple_csq_lists`.

        - If you want to group consequences by gene and determine the most severe
          consequence for each gene, consider using `process_consequences`.

    If `prioritize_protein_coding` is True, protein-coding transcripts are prioritized
    by filtering to only protein-coding transcripts and determining the
    most severe consequence. If no protein-coding transcripts are present, determine
    the most severe consequence for all transcripts.

    If `prioritize_loftee` is True, prioritize consequences with LOFTEE annotations, in
    the order of `loftee_labels`, over those without LOFTEE annotations. If
    `prioritize_loftee_no_flags` is True, prioritize LOFTEE consequences with no flags
    over those with flags.

    If `additional_order` is provided, additional ordering is applied to the
    consequences in the list after any of the above prioritization. An example use of
    this parameter is to prioritize by PolyPhen predictions.

    :param csq_expr: ArrayExpression of VEP consequences to filter.
    :param csq_order: List indicating the order of VEP consequences, sorted from high to
        low impact. Default is None, which uses the value of the `CSQ_ORDER` global.
    :param loftee_labels: Annotations added by LOFTEE, sorted from high to low impact.
        Default is None, which uses the value of the `LOFTEE_LABELS` global.
    :param prioritize_protein_coding: Whether to prioritize protein-coding transcripts
        when determining the worst consequence. Default is False.
    :param prioritize_loftee: Whether to prioritize LOFTEE consequences. Default is
        False.
    :param prioritize_loftee_no_flags: Whether to prioritize LOFTEE consequences with no
        flags over those with flags. Default is False.
    :param additional_order_field: Field name of the consequence annotation to use for
        additional ordering to apply to the consequences in the list. Default is None.
    :param additional_order: The ordering to use for prioritizing consequences in the
        `additional_order_field`. Default is None.
    :return: ArrayExpression with of the consequences that match the most severe
        consequence.
    """
    # Get the dtype of the csq_expr ArrayExpression elements
    csq_type = csq_expr.dtype.element_type

    if ((additional_order_field is None) + (additional_order is None)) == 1:
        raise ValueError(
            "If `additional_order_field` is provided, `additional_order` must also be"
            " provided and vice versa."
        )

    if additional_order_field and additional_order_field not in csq_type.fields:
        raise ValueError("Additional order field not found in consequence type.")

    # Define the order of fields to prioritize by, based on specified parameters.
    priority_order = (
        (["protein_coding"] if prioritize_protein_coding else [])
        + (["lof"] if prioritize_loftee else [])
        + (["no_lof_flags"] if prioritize_loftee_no_flags else [])
        + ["most_severe_consequence"]
    )

    # Define the impact ordering of VEP consequences and LOFTEE labels. If not provided,
    # use the globals CSQ_ORDER (set in get_most_severe_consequence_expr) and
    # LOFTEE_LABELS.
    loftee_labels = loftee_labels or LOFTEE_LABELS
    term_order = {"most_severe_consequence": csq_order, "lof": loftee_labels}

    # Add the additional order field to the priority order and term order if provided.
    if additional_order_field:
        priority_order.append(additional_order_field)
        term_order[additional_order_field] = additional_order

    # Define initial result and current expression.
    result = {}
    order_result_fields = ["most_severe_consequence", "lof"]
    curr_expr = hl.or_missing(hl.len(csq_expr) > 0, csq_expr)

    for curr_field in priority_order:
        order = term_order.get(curr_field)
        if order is None and curr_field != "most_severe_consequence":
            # If there is no order specified for the current field, then the field is
            # used as a parameter to filter_vep_transcript_csqs_expr and if there are
            # any consequences remaining, the result is set to True.
            curr_expr = filter_vep_transcript_csqs_expr(curr_expr, **{curr_field: True})
            result[curr_field] = hl.len(curr_expr) > 0
        else:
            # Handle the case where the current field is a collection of consequences
            # each with a 'consequence_terms' field (e.g. transcript_consequences) that
            # need to be flattened to determine the most severe consequence.
            f = curr_field if curr_field in csq_type.fields else "consequence_terms"
            if isinstance(csq_type[f], hl.tarray) or isinstance(csq_type[f], hl.tset):
                f_map = csq_expr.flatmap
                f_func = lambda x, csq: x.contains(csq)
            else:
                f_map = csq_expr.map
                f_func = lambda x, csq: x == csq

            # Get the most severe (highest impact) consequence for the current field.
            ms_csq_expr = get_most_severe_consequence_expr(f_map(lambda x: x[f]), order)

            # Filter to only elements that contain the most severe (highest impact)
            # consequence for the current field, and return missing if the most severe
            # consequence is missing.
            curr_expr = hl.or_missing(
                hl.is_defined(ms_csq_expr),
                csq_expr.filter(lambda x: f_func(x[f], ms_csq_expr)),
            )

            # Add the most severe consequence to the result if the field is in the order
            # result fields. When there is no most severe consequence and the current
            # field has a result expression, the result is kept as the existing result
            # value.
            if curr_field in order_result_fields:
                if curr_field in result:
                    ms_csq_expr = hl.or_else(ms_csq_expr, result[curr_field])
                result[curr_field] = ms_csq_expr

            if curr_field == "lof":
                result["no_lof_flags"] = hl.any(
                    curr_expr.map(
                        lambda x: hl.is_missing(x.lof_flags) | (x.lof_flags == "")
                    )
                )

        curr_expr = hl.or_missing(hl.len(curr_expr) > 0, curr_expr)
        csq_expr = hl.or_else(curr_expr, csq_expr)

    return hl.struct(**result, consequences=csq_expr)


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
    csq_order: Optional[List[str]] = None,
    loftee_labels: Optional[List[str]] = None,
    include_csqs: bool = False,
    prioritize_protein_coding: bool = True,
    prioritize_loftee: bool = True,
    prioritize_loftee_no_flags: bool = False,
    csq_list_order: Union[List[str], Tuple[str]] = (
        "transcript_consequences",
        "regulatory_feature_consequences",
        "motif_feature_consequences",
        "intergenic_consequences",
    ),
    add_order_by_csq_list: Dict[str, Tuple[str, List[str]]] = None,
) -> hl.Table:
    """
    Process multiple VEP consequences lists to determine the most severe consequence.

    Adds the following annotations:
        - most_severe_consequence: Most severe consequence for variant.
        - protein_coding: Whether the variant is present on a protein-coding transcript.
        - lof: Whether the variant is a loss-of-function variant.
        - no_lof_flags: Whether the variant has any LOFTEE flags (True if no flags).

    .. note::

        Assumes input Table is annotated with VEP and that VEP annotations have been
        filtered to canonical or MANE Select transcripts if wanted.

    If `include_csqs` is True, additional annotations are added for each VEP
    consequences list in `csq_list_order`, with the consequences that match the most
    severe consequence term.

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

    If `prioritize_loftee` is True, prioritize consequences with LOFTEE annotations, in
    the order of `loftee_labels`, over those without LOFTEE annotations. If
    `prioritize_loftee_no_flags` is True, prioritize LOFTEE consequences with no flags
    over those with flags.

    :param vep_expr: StructExpression of VEP consequences to get the most severe
        consequence from.
    :param csq_order: Order of VEP consequences, sorted from high to low impact. Default
        is None, which uses the value of the `CSQ_ORDER` global.
    :param loftee_labels: Annotations added by LOFTEE, sorted from high to low impact.
        Default is None, which uses the value of the `LOFTEE_LABELS` global.
    :param include_csqs: Whether to include all consequences for the most severe
        consequence. Default is False.
    :param prioritize_protein_coding: Whether to prioritize protein-coding transcripts
        when determining the worst consequence. Default is True.
    :param prioritize_loftee: Whether to prioritize consequences with LOFTEE annotations
        over those without. Default is True.
    :param prioritize_loftee_no_flags: Whether to prioritize LOFTEE annotated
        consequences with no flags over those with flags. Default is False.
    :param csq_list_order: Order of VEP consequences lists to be processed. Default is
        ('transcript_consequences', 'regulatory_feature_consequences',
        'motif_feature_consequences', 'intergenic_consequences').
    :param add_order_by_csq_list: Dictionary of additional ordering for VEP consequences
        lists. The key is the name of the consequences list and the value is the order
        of consequences, sorted from high to low impact. Default is None.
    :return: Table annotated with VEP summary annotations.
    """
    add_order_by_csq_list = add_order_by_csq_list or {}
    loftee_labels = loftee_labels or LOFTEE_LABELS

    result = {
        **({"protein_coding": hl.tbool} if prioritize_protein_coding else {}),
        **({"lof": hl.tstr} if prioritize_loftee else {}),
        **(
            {"no_lof_flags": hl.tbool}
            if prioritize_loftee or prioritize_protein_coding
            else {}
        ),
        "most_severe_consequence": hl.tstr,
    }
    result = hl.struct(**{k: hl.missing(v) for k, v in result.items()})

    # Create a struct with missing values for each VEP consequences list.
    ms_csq_list_expr = hl.struct(
        **{c: hl.missing(vep_expr[c].dtype) for c in csq_list_order if c in vep_expr}
    )

    # Build the case expression to determine the most severe consequence.
    ms_csq_expr = hl.case(missing_false=True)
    for csq_list in csq_list_order:
        if csq_list not in vep_expr:
            logger.warning("VEP consequences list %s not found in input!", csq_list)
            continue

        is_tc = csq_list == "transcript_consequences"
        csq_expr = vep_expr[csq_list]

        # Set the base arguments for filtering to the most severe consequence using
        # filter_to_most_severe_consequences.
        add_order = add_order_by_csq_list.get(csq_list)
        base_args = {
            "csq_order": csq_order,
            "loftee_labels": loftee_labels,
            "prioritize_protein_coding": (
                True if (prioritize_protein_coding and is_tc) else False
            ),
            "prioritize_loftee": True if (prioritize_loftee and is_tc) else False,
            "prioritize_loftee_no_flags": (
                True if (prioritize_loftee_no_flags and is_tc) else False
            ),
            "additional_order_field": add_order[0] if add_order else None,
            "additional_order": add_order[1] if add_order else None,
        }
        ms_expr = filter_to_most_severe_consequences(csq_expr, **base_args)
        ms_expr = result.annotate(**ms_expr)

        # Annotate the current consequence list with the consequences that match the
        # most severe consequence term.
        if include_csqs:
            ms_expr = ms_expr.annotate(
                **ms_csq_list_expr.annotate(**{csq_list: ms_expr.consequences})
            )
        ms_expr = ms_expr.drop("consequences")

        # If the length of the consequence list is not 0, set the most severe
        # consequence to the most severe consequence for the current list.
        ms_csq_expr = ms_csq_expr.when(hl.len(csq_expr) > 0, ms_expr)

    return ms_csq_expr.or_missing()


def filter_vep_transcript_csqs(
    t: Union[hl.Table, hl.MatrixTable],
    vep_root: str = "vep",
    synonymous: bool = True,
    canonical: bool = True,
    ensembl_only: bool = True,
    filter_empty_csq: bool = True,
    **kwargs,
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Filter VEP transcript consequences based on specified criteria, and optionally filter to variants where transcript consequences is not empty after filtering.

    If `filter_empty_csq` parameter is set to True, the Table/MatrixTable is filtered
    to variants where 'transcript_consequences' within the VEP annotation is not empty
    after the specified filtering criteria is applied.

    .. note::

        By default, the Table/MatrixTable is filtered to variants where
        'transcript_consequences' within the VEP annotation is not empty after filtering
        to Ensembl canonical transcripts with a most severe consequence of
        'synonymous_variant'.

    :param t: Input Table or MatrixTable.
    :param vep_root: Root for VEP annotation. Default is 'vep'.
    :param synonymous: Whether to filter to variants where the most severe consequence
        is 'synonymous_variant'. Default is True.
    :param canonical: Whether to filter to only canonical transcripts. Default is True.
    :param ensembl_only: Whether to filter to only Ensembl transcripts. This option is
        useful for deduplicating transcripts that are the same between RefSeq and
        Emsembl. Default is True.
    :param filter_empty_csq: Whether to filter out rows where 'transcript_consequences'
        is empty, after filtering 'transcript_consequences' to the specified criteria.
        Default is True.
    :param kwargs: Filtering criteria to apply to the VEP transcript consequences using
        `filter_vep_transcript_csqs_expr`. See that function for more details.
    :return: Table or MatrixTable with VEP transcript consequences filtered.
    """
    is_mt = isinstance(t, hl.MatrixTable)
    vep_data = {
        vep_root: t[vep_root].annotate(
            transcript_consequences=filter_vep_transcript_csqs_expr(
                t[vep_root].transcript_consequences,
                synonymous=synonymous,
                canonical=canonical,
                ensembl_only=ensembl_only,
                **kwargs,
            )
        )
    }
    t = t.annotate_rows(**vep_data) if is_mt else t.annotate(**vep_data)

    if filter_empty_csq:
        transcript_csq_expr = t[vep_root].transcript_consequences
        filter_expr = hl.is_defined(transcript_csq_expr) & (
            hl.len(transcript_csq_expr) > 0
        )
        t = t.filter_rows(filter_expr) if is_mt else t.filter(filter_expr)

    return t


def filter_vep_transcript_csqs_expr(
    csq_expr: hl.expr.ArrayExpression,
    synonymous: bool = False,
    canonical: bool = False,
    mane_select: bool = False,
    ensembl_only: bool = False,
    protein_coding: bool = False,
    loftee_labels: Optional[List[str]] = None,
    no_lof_flags: bool = False,
    csqs: List[str] = None,
    keep_csqs: bool = True,
    genes: Optional[List[str]] = None,
    keep_genes: bool = True,
    match_by_gene_symbol: bool = False,
    additional_filtering_criteria: Optional[List[Callable]] = None,
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Filter VEP transcript consequences based on specified criteria, and optionally filter to variants where transcript consequences is not empty after filtering.

    .. note::

        If `csqs` is not None or `synonymous` is True, and 'most_severe_consequence'
        is not already annotated on the `csq_expr` elements, the most severe
        consequence will be added to the `csq_expr` for filtering.

    :param csq_expr: ArrayExpression of VEP transcript consequences.
    :param synonymous: Whether to filter to variants where the most severe consequence
        is 'synonymous_variant'. Default is False.
    :param canonical: Whether to filter to only canonical transcripts. Default is False.
    :param mane_select: Whether to filter to only MANE Select transcripts. Default is
        False.
    :param ensembl_only: Whether to filter to only Ensembl transcripts. This option is
        useful for deduplicating transcripts that are the same between RefSeq and
        Emsembl. Default is False.
    :param protein_coding: Whether to filter to only protein-coding transcripts.
        Default is False.
    :param loftee_labels: List of LOFTEE labels to filter to. Default is None, which
        filters to all LOFTEE labels.
    :param no_lof_flags: Whether to filter to consequences with no LOFTEE flags.
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
    :return: ArrayExpression of filtered VEP transcript consequences.
    """
    csq_fields = csq_expr.dtype.element_type.fields
    criteria = [lambda csq: True]
    if synonymous:
        logger.info("Filtering to most severe consequence of synonymous_variant...")
        csqs = ["synonymous_variant"]
    if csqs is not None:
        if "most_severe_consequence" not in csq_fields:
            logger.info("Adding most_severe_consequence annotation...")
            csq_expr = add_most_severe_consequence_to_consequence(csq_expr)

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
    if loftee_labels:
        logger.info(
            "Filtering to consequences with LOFTEE labels: %s...", loftee_labels
        )
        criteria.append(lambda x: hl.set(loftee_labels).contains(x.lof))
    if no_lof_flags:
        logger.info("Filtering to consequences with no LOFTEE flags...")
        if "lof_flags" in csq_fields:
            criteria.append(lambda x: hl.is_missing(x.lof_flags) | (x.lof_flags == ""))
        else:
            logger.warning(
                "'lof_flags' not present in consequence struct, no consequences are filtered based on  LOFTEE flags"
            )
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

    if len(criteria) == 1:
        logger.warning("No changes have been made to input transcript consequences!")

    return csq_expr.filter(lambda x: combine_functions(criteria, x))


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
