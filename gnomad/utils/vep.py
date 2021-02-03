import json
import logging
import os
import subprocess
from typing import List, Optional, Union

import hail as hl

from gnomad.resources.resource_utils import VersionedTableResource

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Note that this is the current as of v81 with some included for backwards compatibility (VEP <= 75)
CSQ_CODING_HIGH_IMPACT = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
]

CSQ_CODING_MEDIUM_IMPACT = [
    "start_lost",  # new in v81
    "initiator_codon_variant",  # deprecated
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",  # new in v79
    "splice_region_variant",
]

CSQ_CODING_LOW_IMPACT = [
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

POSSIBLE_REFS = ("GRCh37", "GRCh38")
"""
Constant containing supported references
"""

VEP_CONFIG_PATH = "file:///vep_data/vep-gcloud.json"
"""
Constant that contains the local path to the VEP config file
"""

VEP_CSQ_FIELDS = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info"
"""
Constant that defines the order of VEP annotations used in VCF export.
"""

VEP_CSQ_HEADER = f"Consequence annotations from Ensembl VEP. Format: {VEP_CSQ_FIELDS}"
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
    Returns the output of vep --help which includes the VEP version.

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
    Get VEP context resource for the genome build `ref`

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
    VEP a table, or lookup variants in a reference database

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
                f'vep_or_lookup_vep got {reference}. Expected one of {", ".join(POSSIBLE_REFS)}'
            )

        vep_context = get_vep_context(reference)
        if vep_version is None:
            vep_version = vep_context.default_version

        if vep_version not in vep_context.versions:
            logger.warning(
                f"No VEPed context Table available for genome build {reference} and VEP version {vep_version}, "
                f"all variants will be VEPed using the following VEP:\n{vep_help}"
            )
            return hl.vep(ht, vep_config_path)

        logger.info(
            f"Using VEPed context Table from genome build {reference} and VEP version {vep_version}"
        )

        reference_vep_ht = vep_context.versions[vep_version].ht()
        vep_context_help = hl.eval(reference_vep_ht.vep_help)
        vep_context_config = hl.eval(reference_vep_ht.vep_config)

        assert vep_help == vep_context_help, (
            f"The VEP context HT version does not match the version referenced in the VEP config file."
            f"\nVEP context:\n{vep_context_help}\n\n VEP config:\n{vep_help}"
        )

        assert vep_config == vep_context_config, (
            f"The VEP context HT configuration does not match the configuration in {vep_config_path}."
            f"\nVEP context:\n{vep_context_config}\n\n Current config:\n{vep_config}"
        )

    ht = ht.annotate(vep=reference_vep_ht[ht.key].vep)

    vep_ht = ht.filter(hl.is_defined(ht.vep))
    revep_ht = ht.filter(hl.is_missing(ht.vep))
    revep_ht = hl.vep(revep_ht, vep_config_path)

    return vep_ht.union(revep_ht)


def add_most_severe_consequence_to_consequence(
    tc: hl.expr.StructExpression,
) -> hl.expr.StructExpression:
    """
    Add most_severe_consequence annotation to transcript consequences.

    This is for a given transcript, as there are often multiple annotations for a single transcript:
    e.g. splice_region_variant&intron_variant -> splice_region_variant
    """
    csqs = hl.literal(CSQ_ORDER)

    return tc.annotate(
        most_severe_consequence=csqs.find(lambda c: tc.consequence_terms.contains(c))
    )


def process_consequences(
    mt: Union[hl.MatrixTable, hl.Table],
    vep_root: str = "vep",
    penalize_flags: bool = True,
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Adds most_severe_consequence (worst consequence for a transcript) into [vep_root].transcript_consequences,
    and worst_csq_by_gene, any_lof into [vep_root]

    :param mt: Input MT
    :param vep_root: Root for vep annotation (probably vep)
    :param penalize_flags: Whether to penalize LOFTEE flagged variants, or treat them as equal to HC
    :return: MT with better formatted consequences
    """
    csqs = hl.literal(CSQ_ORDER)
    csq_dict = hl.literal(dict(zip(CSQ_ORDER, range(len(CSQ_ORDER)))))

    def find_worst_transcript_consequence(
        tcl: hl.expr.ArrayExpression,
    ) -> hl.expr.StructExpression:
        """
        Gets worst transcript_consequence from an array of em
        """
        flag_score = 500
        no_flag_score = flag_score * (1 + penalize_flags)

        def csq_score(tc):
            return csq_dict[csqs.find(lambda x: x == tc.most_severe_consequence)]

        tcl = tcl.map(
            lambda tc: tc.annotate(
                csq_score=hl.case(missing_false=True)
                .when(
                    (tc.lof == "HC") & (tc.lof_flags == ""),
                    csq_score(tc) - no_flag_score,
                )
                .when(
                    (tc.lof == "HC") & (tc.lof_flags != ""), csq_score(tc) - flag_score
                )
                .when(tc.lof == "OS", csq_score(tc) - 20)
                .when(tc.lof == "LC", csq_score(tc) - 10)
                .when(
                    tc.polyphen_prediction == "probably_damaging", csq_score(tc) - 0.5
                )
                .when(
                    tc.polyphen_prediction == "possibly_damaging", csq_score(tc) - 0.25
                )
                .when(tc.polyphen_prediction == "benign", csq_score(tc) - 0.1)
                .default(csq_score(tc))
            )
        )
        return hl.or_missing(hl.len(tcl) > 0, hl.sorted(tcl, lambda x: x.csq_score)[0])

    transcript_csqs = mt[vep_root].transcript_consequences.map(
        add_most_severe_consequence_to_consequence
    )

    gene_dict = transcript_csqs.group_by(lambda tc: tc.gene_symbol)
    worst_csq_gene = gene_dict.map_values(find_worst_transcript_consequence).values()
    sorted_scores = hl.sorted(worst_csq_gene, key=lambda tc: tc.csq_score)

    canonical = transcript_csqs.filter(lambda csq: csq.canonical == 1)
    gene_canonical_dict = canonical.group_by(lambda tc: tc.gene_symbol)
    worst_csq_gene_canonical = gene_canonical_dict.map_values(
        find_worst_transcript_consequence
    ).values()
    sorted_canonical_scores = hl.sorted(
        worst_csq_gene_canonical, key=lambda tc: tc.csq_score
    )

    vep_data = mt[vep_root].annotate(
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
        mt.annotate_rows(**{vep_root: vep_data})
        if isinstance(mt, hl.MatrixTable)
        else mt.annotate(**{vep_root: vep_data})
    )


def filter_vep_to_canonical_transcripts(
    mt: Union[hl.MatrixTable, hl.Table], vep_root: str = "vep"
) -> Union[hl.MatrixTable, hl.Table]:
    canonical = mt[vep_root].transcript_consequences.filter(
        lambda csq: csq.canonical == 1
    )
    vep_data = mt[vep_root].annotate(transcript_consequences=canonical)
    return (
        mt.annotate_rows(**{vep_root: vep_data})
        if isinstance(mt, hl.MatrixTable)
        else mt.annotate(**{vep_root: vep_data})
    )


def filter_vep_to_synonymous_variants(
    mt: Union[hl.MatrixTable, hl.Table], vep_root: str = "vep"
) -> Union[hl.MatrixTable, hl.Table]:
    synonymous = mt[vep_root].transcript_consequences.filter(
        lambda csq: csq.most_severe_consequence == "synonymous_variant"
    )
    vep_data = mt[vep_root].annotate(transcript_consequences=synonymous)
    return (
        mt.annotate_rows(**{vep_root: vep_data})
        if isinstance(mt, hl.MatrixTable)
        else mt.annotate(**{vep_root: vep_data})
    )


def vep_struct_to_csq(
    vep_expr: hl.expr.StructExpression, csq_fields: str = VEP_CSQ_FIELDS
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
    :param csq_fields: The | delimited list of fields to include in the CSQ (in that order)
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
                    else element.regulatory_feature_id
                    if "regulatory_feature_id" in element
                    else element.motif_feature_id
                    if "motif_feature_id" in element
                    else ""
                ),
                "variant_class": vep_expr.variant_class,
            }
        )

        # Add exception for transcripts
        if feature_type == "Transcript":
            fields.update(
                {
                    "canonical": hl.cond(element.canonical == 1, "YES", ""),
                    "ensp": element.protein_id,
                    "gene": element.gene_id,
                    "symbol": element.gene_symbol,
                    "symbol_source": element.gene_symbol_source,
                    "cdna_position": hl.str(element.cdna_start)
                    + hl.cond(
                        element.cdna_start == element.cdna_end,
                        "",
                        "-" + hl.str(element.cdna_end),
                    ),
                    "cds_position": hl.str(element.cds_start)
                    + hl.cond(
                        element.cds_start == element.cds_end,
                        "",
                        "-" + hl.str(element.cds_end),
                    ),
                    "protein_position": hl.str(element.protein_start)
                    + hl.cond(
                        element.protein_start == element.protein_end,
                        "",
                        "-" + hl.str(element.protein_end),
                    ),
                    "sift": element.sift_prediction
                    + "("
                    + hl.format("%.3f", element.sift_score)
                    + ")",
                    "polyphen": element.polyphen_prediction
                    + "("
                    + hl.format("%.3f", element.polyphen_score)
                    + ")",
                    "domains": hl.delimit(
                        element.domains.map(lambda d: d.db + ":" + d.name), "&"
                    ),
                }
            )
        elif feature_type == "MotifFeature":
            fields["motif_score_change"] = hl.format("%.3f", element.motif_score_change)

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


def get_most_severe_consequence_for_summary(
    ht: hl.Table,
    csq_order: List[str] = CSQ_ORDER,
    loftee_labels: List[str] = LOFTEE_LABELS,
) -> hl.Table:
    """
    Prepares hail Table for summary statistics generation 
    (number of variants by functional class, number of variants by functional class per sample).

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

    def _get_most_severe_csq(
        csq_list: hl.expr.ArrayExpression, protein_coding: bool
    ) -> hl.expr.StructExpression:
        """
        Processes VEP consequences to generate summary annotations.

        :param csq_list: VEP consequences list to be processed.
        :param protein_coding: Whether variant is in a protein-coding transcript.
        :return: Struct containing summary annotations.
        """
        lof = hl.null(hl.tstr)
        no_lof_flags = hl.null(hl.tbool)
        if protein_coding:
            all_lofs = csq_list.map(lambda x: x.lof)
            lof = hl.literal(loftee_labels).find(lambda x: all_lofs.contains(x))
            csq_list = hl.if_else(
                hl.is_defined(lof), csq_list.filter(lambda x: x.lof == lof), csq_list
            )
            no_lof_flags = hl.or_missing(
                hl.is_defined(lof),
                csq_list.any(lambda x: (x.lof == lof) & hl.is_missing(x.lof_flags)),
            )
        all_csq_terms = csq_list.flatmap(lambda x: x.consequence_terms)
        most_severe_csq = hl.literal(csq_order).find(
            lambda x: all_csq_terms.contains(x)
        )
        return hl.struct(
            most_severe_csq=most_severe_csq,
            protein_coding=protein_coding,
            lof=lof,
            no_lof_flags=no_lof_flags,
        )

    protein_coding = ht.vep.transcript_consequences.filter(
        lambda x: x.biotype == "protein_coding"
    )
    return ht.annotate(
        **hl.case(missing_false=True)
        .when(hl.len(protein_coding) > 0, _get_most_severe_csq(protein_coding, True))
        .when(
            hl.len(ht.vep.transcript_consequences) > 0,
            _get_most_severe_csq(ht.vep.transcript_consequences, False),
        )
        .when(
            hl.len(ht.vep.regulatory_feature_consequences) > 0,
            _get_most_severe_csq(ht.vep.regulatory_feature_consequences, False),
        )
        .when(
            hl.len(ht.vep.motif_feature_consequences) > 0,
            _get_most_severe_csq(ht.vep.motif_feature_consequences, False),
        )
        .default(_get_most_severe_csq(ht.vep.intergenic_consequences, False))
    )
