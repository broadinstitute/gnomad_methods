from typing import Union

import hail as hl
from gnomad.resources.resource_utils import DataException

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

VEP_REFERENCE_DATA = {
    "GRCh37": {
        "vep_config": "file:///vep_data/vep-gcloud.json",
        "all_possible": "gs://gnomad-public-requester-pays/resources/context/grch37_context_vep_annotated.ht",
    },
    "GRCh38": {
        "vep_config": "file:///vep_data/vep-gcloud.json",
        "all_possible": "gs://gnomad-public-requester-pays/resources/context/grch38_context_vep_annotated.ht",
    },
}


VEP_CSQ_FIELDS = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info"
"""
Constant that defines the order of VEP annotations used in VCF export.
"""

VEP_CSQ_HEADER = f"Consequence annotations from Ensembl VEP. Format: {VEP_CSQ_FIELDS}"
"""
Constant that contains description for VEP used in VCF export.
"""


def vep_context_ht_path(ref: str = "GRCh37"):
    if ref not in VEP_REFERENCE_DATA.keys():
        raise DataException(
            "Select reference as one of: {}".format(",".join(VEP_REFERENCE_DATA.keys()))
        )
    return VEP_REFERENCE_DATA[ref]["all_possible"]


def vep_config_path(ref: str = "GRCh37"):
    if ref not in VEP_REFERENCE_DATA.keys():
        raise DataException(
            "Select reference as one of: {}".format(",".join(VEP_REFERENCE_DATA.keys()))
        )
    return VEP_REFERENCE_DATA[ref]["vep_config"]


def vep_or_lookup_vep(ht, reference_vep_ht=None, reference=None, vep_config=None):
    """
    VEP a table, or lookup variants in a reference database

    :param ht: Input Table
    :param reference_vep_ht: A reference database with VEP annotations (must be in top-level `vep`)
    :param reference: If reference_vep_ht is not specified, find a suitable one in reference (if None, grabs from hl.default_reference)
    :param vep_config: vep_config to pass to hl.vep (if None, a suitable one for `reference` is chosen)
    :return: VEPped Table
    """
    if reference is None:
        reference = hl.default_reference().name
    if reference_vep_ht is None:

        possible_refs = ("GRCh37", "GRCh38")
        if reference not in possible_refs:
            raise ValueError(
                f'vep_or_lookup_vep got {reference}. Expected one of {", ".join(possible_refs)}'
            )

        reference_vep_ht = hl.read_table(vep_context_ht_path(reference))

    ht = ht.annotate(vep=reference_vep_ht[ht.key].vep)

    vep_ht = ht.filter(hl.is_defined(ht.vep))
    revep_ht = ht.filter(hl.is_missing(ht.vep))

    if vep_config is None:
        vep_config = vep_config_path(reference)

    revep_ht = hl.vep(revep_ht, vep_config)

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
