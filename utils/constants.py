POP_NAMES = {
    'AFR': "African/African American",
    'AMR': "Admixed American",
    'ASJ': "Ashkenazi Jewish",
    'EAS': "East Asian",
    'FIN': "Finnish",
    'NFE': "Non-Finnish European",
    'OTH': "Other (population not assigned)",
    'SAS': "South Asian"
}

SEXES = {
    'Male': 'Male',
    'Female': 'Female'
}

# Note that this is the current as of v81 with some included for backwards compatibility (VEP <= 75)
CSQ_CODING_HIGH_IMPACT = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost"]

CSQ_CODING_MEDIUM_IMPACT = [
    "start_lost",  # new in v81
    "initiator_codon_variant",  # deprecated
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",  # new in v79
    "splice_region_variant"
]

CSQ_CODING_LOW_IMPACT = [
    "incomplete_terminal_codon_variant",
    "start_retained_variant",  # new in v92
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant"]

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
    "intergenic_variant"
]

CSQ_ORDER = CSQ_CODING_HIGH_IMPACT + CSQ_CODING_MEDIUM_IMPACT + CSQ_CODING_LOW_IMPACT + CSQ_NON_CODING