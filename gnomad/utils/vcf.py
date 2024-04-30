# noqa: D100

import copy
import itertools
import logging
from typing import Dict, List, Optional, Union

import hail as hl

from gnomad.sample_qc.ancestry import POP_NAMES

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

SORT_ORDER = [
    "subset",
    "downsampling",
    "popmax",
    "grpmax",
    "pop",
    "gen_anc",
    "subpop",
    "sex",
    "group",
]
"""
Order to sort subgroupings during VCF export.
Ensures that INFO labels in VCF are in desired order (e.g., raw_AC_afr_female).
"""

GROUPS = ["adj", "raw"]
"""
Group names used to generate labels for high quality genotypes and all raw genotypes. Used in VCF export.
"""

HISTS = ["gq_hist_alt", "gq_hist_all", "dp_hist_alt", "dp_hist_all", "ab_hist_alt"]
"""
Quality histograms used in VCF export.
"""

FAF_POPS = {
    "v3": ["afr", "amr", "eas", "nfe", "sas"],
    "v4": ["afr", "amr", "eas", "mid", "nfe", "sas"],
}
"""
Global populations that are included in filtering allele frequency (faf) calculations. Used in VCF export.
"""

SEXES = ["XX", "XY"]
"""
Sample sexes used in VCF export.

Used to stratify frequency annotations (AC, AN, AF) for each sex.
Note that sample sexes in gnomAD v3 and earlier were 'male' and 'female'.
"""

AS_FIELDS = [
    "AS_FS",
    "AS_MQ",
    "AS_MQRankSum",
    "AS_pab_max",
    "AS_QUALapprox",
    "AS_QD",
    "AS_ReadPosRankSum",
    "AS_SB_TABLE",
    "AS_SOR",
    "AS_VarDP",
    "InbreedingCoeff",
]
"""
Allele-specific variant annotations.
"""

SITE_FIELDS = [
    "FS",
    "MQ",
    "MQRankSum",
    "QUALapprox",
    "QD",
    "ReadPosRankSum",
    "SB",
    "SOR",
    "VarDP",
]
"""
Site level variant annotations.
"""

ALLELE_TYPE_FIELDS = [
    "allele_type",
    "has_star",
    "n_alt_alleles",
    "original_alleles",
    "variant_type",
    "was_mixed",
]
"""
Allele-type annotations.
"""

REGION_FLAG_FIELDS = ["decoy", "lcr", "nonpar", "non_par", "segdup"]
"""
Annotations about variant region type.

.. note::
    decoy resource files do not currently exist for GRCh38/hg38.
"""

JOINT_REGION_FLAG_FIELDS = [
    "fail_interval_qc",
    "outside_broad_capture_region",
    "outside_ukb_capture_region",
    "outside_broad_calling_region",
    "outside_ukb_calling_region",
    "not_called_in_exomes",
    "not_called_in_genomes",
]
"""
Annotations about variant region type that are specifically created for joint dataset of exomes and genomes from gnomAD v4.1.
"""

RF_FIELDS = [
    "rf_positive_label",
    "rf_negative_label",
    "rf_label",
    "rf_train",
    "rf_tp_probability",
]
"""
Annotations specific to the variant QC using a random forest model.
"""

AS_VQSR_FIELDS = ["AS_culprit", "AS_VQSLOD"]
"""
Allele-specific VQSR annotations.
"""

VQSR_FIELDS = AS_VQSR_FIELDS + ["NEGATIVE_TRAIN_SITE", "POSITIVE_TRAIN_SITE"]
"""
Annotations specific to VQSR.
"""

INFO_VCF_AS_PIPE_DELIMITED_FIELDS = [
    "AS_QUALapprox",
    "AS_VarDP",
    "AS_MQ_DP",
    "AS_RAW_MQ",
    "AS_SB_TABLE",
]

INFO_DICT = {
    "FS": {
        "Description": "Phred-scaled p-value of Fisher's exact test for strand bias"
    },
    "InbreedingCoeff": {
        "Number": "A",
        "Description": (
            "Inbreeding coefficient, the excess heterozygosity at a variant site,"
            " computed as 1 - (the number of heterozygous genotypes)/(the number of"
            " heterozygous genotypes expected under Hardy-Weinberg equilibrium)"
        ),
    },
    "inbreeding_coeff": {
        "Number": "A",
        "Description": (
            "Inbreeding coefficient, the excess heterozygosity at a variant site,"
            " computed as 1 - (the number of heterozygous genotypes)/(the number of"
            " heterozygous genotypes expected under Hardy-Weinberg equilibrium)"
        ),
    },
    "MQ": {
        "Description": (
            "Root mean square of the mapping quality of reads across all samples"
        )
    },
    "MQRankSum": {
        "Description": (
            "Z-score from Wilcoxon rank sum test of alternate vs. reference read"
            " mapping qualities"
        )
    },
    "QD": {
        "Description": (
            "Variant call confidence normalized by depth of sample reads supporting a"
            " variant"
        )
    },
    "ReadPosRankSum": {
        "Description": (
            "Z-score from Wilcoxon rank sum test of alternate vs. reference read"
            " position bias"
        )
    },
    "SOR": {"Description": "Strand bias estimated by the symmetric odds ratio test"},
    "POSITIVE_TRAIN_SITE": {
        "Description": (
            "Variant was used to build the positive training set of high-quality"
            " variants for VQSR"
        )
    },
    "NEGATIVE_TRAIN_SITE": {
        "Description": (
            "Variant was used to build the negative training set of low-quality"
            " variants for VQSR"
        )
    },
    "positive_train_site": {
        "Description": (
            "Variant was used to build the positive training set of high-quality"
            " variants for VQSR"
        )
    },
    "negative_train_site": {
        "Description": (
            "Variant was used to build the negative training set of low-quality"
            " variants for VQSR"
        )
    },
    "BaseQRankSum": {
        "Description": (
            "Z-score from Wilcoxon rank sum test of alternate vs. reference base"
            " qualities"
        ),
    },
    "VarDP": {
        "Description": (
            "Depth over variant genotypes (does not include depth of reference samples)"
        )
    },
    "VQSLOD": {
        "Description": (
            "Log-odds ratio of being a true variant versus being a false positive under"
            " the trained VQSR Gaussian mixture model"
        ),
    },
    "culprit": {
        "Description": "Worst-performing annotation in the VQSR Gaussian mixture model",
    },
    "decoy": {"Description": "Variant falls within a reference decoy region"},
    "lcr": {"Description": "Variant falls within a low complexity region"},
    "nonpar": {
        "Description": (
            "Variant (on sex chromosome) falls outside a pseudoautosomal region"
        )
    },
    "non_par": {
        "Description": (
            "Variant (on sex chromosome) falls outside a pseudoautosomal region"
        )
    },
    "segdup": {"Description": "Variant falls within a segmental duplication region"},
    "fail_interval_qc": {
        "Description": (
            "Less than 85 percent of samples meet 20X coverage if variant is in"
            " autosomal or PAR regions or 10X coverage for non-PAR regions of"
            " chromosomes X and Y."
        )
    },
    "outside_ukb_capture_region": {
        "Description": "Variant falls outside of UK Biobank exome capture regions."
    },
    "outside_broad_capture_region": {
        "Description": "Variant falls outside of Broad exome capture regions."
    },
    "rf_positive_label": {
        "Description": (
            "Variant was labelled as a positive example for training of random forest"
            " model"
        )
    },
    "rf_negative_label": {
        "Description": (
            "Variant was labelled as a negative example for training of random forest"
            " model"
        )
    },
    "rf_label": {"Description": "Random forest training label"},
    "rf_train": {"Description": "Variant was used in training random forest model"},
    "rf_tp_probability": {
        "Description": (
            "Probability of a called variant being a true variant as determined by"
            " random forest model"
        )
    },
    "transmitted_singleton": {
        "Description": (
            "Variant was a callset-wide doubleton that was transmitted within a family"
            " from a parent to a child (i.e., a singleton amongst unrelated samples in"
            " cohort)"
        )
    },
    "sibling_singleton": {
        "Description": (
            "Variant was a callset-wide doubleton that was present only in two siblings"
            " (i.e., a singleton amongst unrelated samples in cohort)."
        )
    },
    "original_alleles": {"Description": "Alleles before splitting multiallelics"},
    "variant_type": {
        "Description": "Variant type (snv, indel, multi-snv, multi-indel, or mixed)"
    },
    "allele_type": {
        "Description": "Allele type (snv, insertion, deletion, or mixed)",
    },
    "n_alt_alleles": {
        "Number": "1",
        "Description": "Total number of alternate alleles observed at variant locus",
    },
    "was_mixed": {"Description": "Variant type was mixed"},
    "has_star": {
        "Description": (
            "Variant locus coincides with a spanning deletion (represented by a star)"
            " observed elsewhere in the callset"
        )
    },
    "AS_pab_max": {
        "Number": "A",
        "Description": (
            "Maximum p-value over callset for binomial test of observed allele balance"
            " for a heterozygous genotype, given expectation of 0.5"
        ),
    },
    "monoallelic": {
        "Description": "All samples are homozygous alternate for the variant"
    },
    "only_het": {"Description": "All samples are heterozygous for the variant"},
    "QUALapprox": {
        "Number": "1",
        "Description": "Sum of PL[0] values; used to approximate the QUAL score",
    },
    "AS_SB_TABLE": {
        "Number": ".",
        "Description": (
            "Allele-specific forward/reverse read counts for strand bias tests"
        ),
    },
}
"""
Dictionary used during VCF export to export row (variant) annotations.
"""

JOINT_REGION_FLAGS_INFO_DICT = {
    "fail_interval_qc": {
        "Description": (
            "Less than 85 percent of samples meet 20X coverage if variant is in"
            " autosomal or PAR regions or 10X coverage for non-PAR regions of"
            " chromosomes X and Y."
        )
    },
    "outside_ukb_capture_region": {
        "Description": "Variant falls outside of the UK Biobank exome capture regions."
    },
    "outside_broad_capture_region": {
        "Description": "Variant falls outside of the Broad exome capture regions."
    },
    "outside_ukb_calling_region": {
        "Description": (
            "Variant falls outside of the UK Biobank exome capture regions plus 150 bp"
            " padding."
        )
    },
    "outside_broad_calling_region": {
        "Description": (
            "Variant falls outside of the Broad exome capture regions plus 150 bp"
            " padding."
        )
    },
    "not_called_in_exomes": {
        "Description": "Variant was not called in the gnomAD exomes."
    },
    "not_called_in_genomes": {
        "Description": "Variant was not called in the gnomAD genomes."
    },
}


IN_SILICO_ANNOTATIONS_INFO_DICT = {
    "cadd_raw_score": {
        "Number": "1",
        "Description": (
            "Raw CADD scores are interpretable as the extent to which the annotation"
            " profile for a given variant suggests that the variant is likely to be"
            " 'observed' (negative values) vs 'simulated' (positive values). Larger"
            " values are more deleterious."
        ),
    },
    "cadd_phred": {
        "Number": "1",
        "Description": (
            "Cadd Phred-like scores ('scaled C-scores') ranging from 1 to 99, based on"
            " the rank of each variant relative to all possible 8.6 billion"
            " substitutions in the human reference genome. Larger values are more"
            " deleterious."
        ),
    },
    "revel_max": {
        "Number": "1",
        "Description": (
            "The maximum REVEL score at a site's MANE Select or canonical"
            " transcript. It's an ensemble score for predicting the pathogenicity of"
            " missense variants (based on 13 other variant predictors). Scores ranges"
            " from 0 to 1. Variants with higher scores are predicted to be more likely"
            " to be deleterious."
        ),
    },
    "spliceai_ds_max": {
        "Number": "1",
        "Description": (
            "Illumina's SpliceAI max delta score; interpreted as the probability of the"
            " variant being splice-altering."
        ),
    },
    "pangolin_largest_ds": {
        "Number": "1",
        "Description": (
            "Pangolin's largest delta score across 2 splicing consequences, which"
            " reflects the probability of the variant being splice-altering"
        ),
    },
    "phylop": {
        "Number": "1",
        "Description": (
            "Base-wise conservation score across the 241 placental mammals in the"
            " Zoonomia project. Score ranges from -20 to 9.28, and reflects"
            " acceleration (faster evolution than expected under neutral drift,"
            " assigned negative scores) as well as conservation (slower than expected"
            " evolution, assigned positive scores)."
        ),
    },
    "sift_max": {
        "Number": "1",
        "Description": (
            "Score reflecting the scaled probability of the amino acid substitution"
            " being tolerated, ranging from 0 to 1. Scores below 0.05 are predicted to"
            " impact protein function. We prioritize max scores for MANE Select"
            " transcripts where possible and otherwise report a score for the canonical"
            " transcript."
        ),
    },
    "polyphen_max": {
        "Number": "1",
        "Description": (
            "Score that predicts the possible impact of an amino acid substitution on"
            " the structure and function of a human protein, ranging from 0.0"
            " (tolerated) to 1.0 (deleterious).  We prioritize max scores for MANE"
            " Select transcripts where possible and otherwise report a score for the"
            " canonical transcript."
        ),
    },
}
"""
Dictionary with in silico score descriptions to include in the VCF INFO header.
"""


VRS_FIELDS_DICT = {
    "VRS_Allele_IDs": {
        "Number": "R",
        "Description": (
            "The computed identifiers for the GA4GH VRS Alleles corresponding to the"
            " values in the REF and ALT fields"
        ),
    },
    "VRS_Starts": {
        "Number": "R",
        "Description": (
            "Interresidue coordinates used as the location starts for the GA4GH VRS"
            " Alleles corresponding to the values in the REF and ALT fields"
        ),
    },
    "VRS_Ends": {
        "Number": "R",
        "Description": (
            "Interresidue coordinates used as the location ends for the GA4GH VRS"
            " Alleles corresponding to the values in the REF and ALT fields"
        ),
    },
    "VRS_States": {
        "Number": ".",
        "Description": (
            "The literal sequence states used for the GA4GH VRS Alleles corresponding"
            " to the values in the REF and ALT fields"
        ),
    },
}
"""
Dictionary with VRS annotations to include in the VCF INFO field and VCF header.
"""


ENTRIES = ["GT", "GQ", "DP", "AD", "MIN_DP", "PGT", "PID", "PL", "SB"]
"""
Densified entries to be selected during VCF export.
"""

SPARSE_ENTRIES = [
    "END",
    "DP",
    "GQ",
    "LA",
    "LAD",
    "LGT",
    "LPGT",
    "LPL",
    "MIN_DP",
    "PID",
    "RGQ",
    "SB",
]
"""
Sparse entries to be selected and densified during VCF export.
"""

FORMAT_DICT = {
    "GT": {"Description": "Genotype", "Number": "1", "Type": "String"},
    "AD": {
        "Description": "Allelic depths for the ref and alt alleles in the order listed",
        "Number": "R",
        "Type": "Integer",
    },
    "DP": {
        "Description": (
            "Approximate read depth (reads with MQ=255 or with bad mates are filtered)"
        ),
        "Number": "1",
        "Type": "Integer",
    },
    "GQ": {
        "Description": (
            "Phred-scaled confidence that the genotype assignment is correct. Value is"
            " the difference between the second lowest PL and the lowest PL (always"
            " normalized to 0)."
        ),
        "Number": "1",
        "Type": "Integer",
    },
    "MIN_DP": {
        "Description": "Minimum DP observed within the GVCF block",
        "Number": "1",
        "Type": "Integer",
    },
    "PGT": {
        "Description": (
            "Physical phasing haplotype information, describing how the alternate"
            " alleles are phased in relation to one another"
        ),
        "Number": "1",
        "Type": "String",
    },
    "PID": {
        "Description": (
            "Physical phasing ID information, where each unique ID within a given"
            " sample (but not across samples) connects records within a phasing group"
        ),
        "Number": "1",
        "Type": "String",
    },
    "PL": {
        "Description": (
            "Normalized, phred-scaled likelihoods for genotypes as defined in the VCF"
            " specification"
        ),
        "Number": "G",
        "Type": "Integer",
    },
    "SB": {
        "Description": (
            "Per-sample component statistics which comprise the Fisher's exact test to"
            " detect strand bias. Values are: depth of reference allele on forward"
            " strand, depth of reference allele on reverse strand, depth of alternate"
            " allele on forward strand, depth of alternate allele on reverse strand."
        ),
        "Number": "4",
        "Type": "Integer",
    },
}
"""
Dictionary used during VCF export to export MatrixTable entries.
"""

JOINT_FILTERS_INFO_DICT = {
    "exomes_filters": {"Description": "Filters' values from the exomes dataset."},
    "genomes_filters": {"Description": "Filters' values from the genomes dataset."},
}


def adjust_vcf_incompatible_types(
    ht: hl.Table,
    pipe_delimited_annotations: List[str] = INFO_VCF_AS_PIPE_DELIMITED_FIELDS,
) -> hl.Table:
    """
    Create a Table ready for vcf export.

    In particular, the following conversions are done:
        - All int64 are coerced to int32
        - Fields specified by `pipe_delimited_annotations` are converted from arrays to pipe-delimited strings

    :param ht: Input Table.
    :param pipe_delimited_annotations: List of info fields (they must be fields of the ht.info Struct).
    :return: Table ready for VCF export.
    """

    def get_pipe_expr(array_expr: hl.expr.ArrayExpression) -> hl.expr.StringExpression:
        return hl.delimit(array_expr.map(lambda x: hl.or_else(hl.str(x), "")), "|")

    # Make sure the HT is keyed by locus, alleles
    ht = ht.key_by("locus", "alleles")

    info_type_convert_expr = {}
    # Convert int64 fields to int32 (int64 isn't supported by VCF)
    for f, ft in ht.info.dtype.items():
        if ft == hl.dtype("int64"):
            logger.warning(
                "Coercing field info.%s from int64 to int32 for VCF output. Value"
                " will be capped at int32 max value.",
                f,
            )
            info_type_convert_expr.update(
                {
                    f: hl.or_missing(
                        hl.is_defined(ht.info[f]),
                        hl.int32(hl.min(2**31 - 1, ht.info[f])),
                    )
                }
            )
        elif ft == hl.dtype("array<int64>"):
            logger.warning(
                "Coercing field info.%s from array<int64> to array<int32> for VCF"
                " output. Array values will be capped at int32 max value.",
                f,
            )
            info_type_convert_expr.update(
                {
                    f: ht.info[f].map(
                        lambda x: hl.or_missing(
                            hl.is_defined(x), hl.int32(hl.min(2**31 - 1, x))
                        )
                    )
                }
            )

    ht = ht.annotate(info=ht.info.annotate(**info_type_convert_expr))

    info_expr = {}

    # Make sure to pipe-delimit fields that need to.
    # Note: the expr needs to be prefixed by "|" because GATK expect one value for the ref (always empty)
    # Note2: this doesn't produce the correct annotation for AS_SB_TABLE, it
    # is handled below
    for f in pipe_delimited_annotations:
        if f in ht.info and f != "AS_SB_TABLE":
            info_expr[f] = "|" + get_pipe_expr(ht.info[f])

    # Flatten SB if it is an array of arrays
    if "SB" in ht.info and not isinstance(ht.info.SB, hl.expr.ArrayNumericExpression):
        info_expr["SB"] = ht.info.SB[0].extend(ht.info.SB[1])

    if "AS_SB_TABLE" in ht.info:
        info_expr["AS_SB_TABLE"] = get_pipe_expr(
            ht.info.AS_SB_TABLE.map(lambda x: hl.delimit(x, ","))
        )

    # Annotate with new expression
    ht = ht.annotate(info=ht.info.annotate(**info_expr))

    return ht


def make_label_combos(
    label_groups: Dict[str, List[str]],
    sort_order: List[str] = SORT_ORDER,
    label_delimiter: str = "_",
) -> List[str]:
    """
    Make combinations of all possible labels for a supplied dictionary of label groups.

    For example, if label_groups is `{"sex": ["male", "female"], "pop": ["afr", "nfe", "amr"]}`,
    this function will return `["afr_male", "afr_female", "nfe_male", "nfe_female", "amr_male", "amr_female']`

    :param label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"]).
    :param sort_order: List containing order to sort label group combinations. Default is SORT_ORDER.
    :param label_delimiter: String to use as delimiter when making group label combinations.
    :return: list of all possible combinations of values for the supplied label groupings.
    """
    copy_label_groups = copy.deepcopy(label_groups)
    if len(copy_label_groups) == 1:
        return [item for sublist in copy_label_groups.values() for item in sublist]
    anchor_group = sorted(copy_label_groups.keys(), key=lambda x: sort_order.index(x))[
        0
    ]
    anchor_val = copy_label_groups.pop(anchor_group)
    combos = []
    for x, y in itertools.product(
        anchor_val,
        make_label_combos(copy_label_groups, label_delimiter=label_delimiter),
    ):
        combos.append(f"{x}{label_delimiter}{y}")
    return combos


def index_globals(
    globals_array: List[Dict[str, str]],
    label_groups: Dict[str, List[str]],
    label_delimiter: str = "_",
) -> Dict[str, int]:
    """
    Create a dictionary keyed by the specified label groupings with values describing the corresponding index of each grouping entry in the meta_array annotation.

    :param globals_array: Ordered list containing dictionary entries describing all the grouping combinations contained in the globals_array annotation.
       Keys are the grouping type (e.g., 'group', 'pop', 'sex') and values are the grouping attribute (e.g., 'adj', 'eas', 'XY').
    :param label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"])
    :param label_delimiter: String used as delimiter when making group label combinations.
    :return: Dictionary keyed by specified label grouping combinations, with values describing the corresponding index
        of each grouping entry in the globals
    """
    combos = make_label_combos(label_groups, label_delimiter=label_delimiter)
    index_dict = {}

    for combo in combos:
        combo_fields = combo.split(label_delimiter)
        for i, v in enumerate(globals_array):
            if set(v.values()) == set(combo_fields):
                index_dict.update({f"{combo}": i})
    return index_dict


def make_combo_header_text(
    preposition: str,
    combo_dict: Dict[str, str],
    pop_names: Dict[str, str],
) -> str:
    """
    Programmatically generate text to populate the VCF header description for a given variant annotation with specific groupings and subset.

    For example, if preposition is "for", group_types is ["group", "pop", "sex"], and combo_fields is ["adj", "afr", "female"],
    this function will return the string " for female samples in the African-American/African genetic ancestry group".

    :param preposition: Relevant preposition to precede automatically generated text.
    :param combo_dict: Dict with grouping types as keys and values for grouping type as values. This function generates text for these values.
        Possible grouping types are: "group", "pop", "sex", and "subpop".
        Example input: {"pop": "afr", "sex": "female"}
    :param pop_names: Dict with global population names (keys) and population descriptions (values).
    :return: String with automatically generated description text for a given set of combo fields.
    """
    header_text = " " + preposition

    if len(combo_dict) == 1:
        if combo_dict["group"] == "adj":
            return ""

    if "sex" in combo_dict:
        header_text = header_text + " " + combo_dict["sex"]

    header_text = header_text + " samples"

    if "subpop" in combo_dict or "pop" in combo_dict:
        if "subpop" in combo_dict:
            header_text = (
                header_text
                + f" in the {pop_names[combo_dict['subpop']]} genetic ancestry subgroup"
            )

        else:
            header_text = (
                header_text
                + f" in the {pop_names[combo_dict['pop']]} genetic ancestry group"
            )

    if "group" in combo_dict:
        if combo_dict["group"] == "raw":
            header_text = header_text + ", before removing low-confidence genotypes"

    return header_text


def create_label_groups(
    pops: List[str],
    sexes: List[str] = SEXES,
    all_groups: List[str] = GROUPS,
    pop_sex_groups: List[str] = ["adj"],
) -> List[Dict[str, List[str]]]:
    """
    Generate a list of label group dictionaries needed to populate info dictionary.

    Label dictionaries are passed as input to `make_info_dict`.

    :param pops: List of population names.
    :param sexes: List of sample sexes.
    :param all_groups: List of data types (raw, adj). Default is `GROUPS`, which is ["raw", "adj"].
    :param pop_sex_groups: List of data types (raw, adj) to populate with pops and sexes. Default is ["adj"].
    :return: List of label group dictionaries.
    """
    return [
        # This is to capture raw frequency fields, which are
        # not stratified by sex or population (e.g., only AC_raw exists, not AC_XX_raw)
        dict(group=all_groups),
        dict(group=pop_sex_groups, sex=sexes),
        dict(group=pop_sex_groups, pop=pops),
        dict(group=pop_sex_groups, pop=pops, sex=sexes),
    ]


def make_info_dict(
    prefix: str = "",
    suffix: str = "",
    prefix_before_metric: bool = True,
    pop_names: Dict[str, str] = POP_NAMES,
    label_groups: Dict[str, List[str]] = None,
    label_delimiter: str = "_",
    bin_edges: Dict[str, str] = None,
    faf: bool = False,
    popmax: bool = False,
    grpmax: bool = False,
    fafmax: bool = False,
    callstats: bool = False,
    freq_ctt: bool = False,
    freq_cmh: bool = False,
    freq_stat_union: bool = False,
    description_text: str = "",
    age_hist_distribution: str = None,
    sort_order: List[str] = SORT_ORDER,
) -> Dict[str, Dict[str, str]]:
    """
    Generate dictionary of Number and Description attributes of VCF INFO fields.

    Used to populate the INFO fields of the VCF header during export.

    Creates:
        - INFO fields for age histograms (bin freq, n_smaller, and n_larger for heterozygous and homozygous variant carriers)
        - INFO fields for popmax AC, AN, AF, nhomalt, and popmax population
        - INFO fields for AC, AN, AF, nhomalt for each combination of sample population, sex, and subpopulation, both for adj and raw data
        - INFO fields for filtering allele frequency (faf) annotations

    :param prefix: Prefix string for data, e.g. "gnomAD". Default is empty string.
    :param suffix: Suffix string for data, e.g. "gnomAD". Default is empty string.
    :param prefix_before_metric: Whether prefix should be added before the metric (AC, AN, AF, nhomalt, faf95, faf99) in INFO field. Default is True.
    :param pop_names: Dict with global population names (keys) and population descriptions (values). Default is POP_NAMES.
    :param label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"]).
    :param label_delimiter: String to use as delimiter when making group label combinations.
    :param bin_edges: Dictionary keyed by annotation type, with values that reflect the bin edges corresponding to the annotation.
    :param faf: If True, use alternate logic to auto-populate dictionary values associated with filter allele frequency annotations.
    :param popmax: If True, use alternate logic to auto-populate dictionary values associated with popmax annotations.
    :param grpmax: If True, use alternate logic to auto-populate dictionary values associated with grpmax annotations.
    :param fafmax: If True, use alternate logic to auto-populate dictionary values associated with fafmax annotations.
    :param callstats: If True, use alternate logic to auto-populate dictionary values associated with callstats annotations.
    :param freq_ctt: If True, use alternate logic to auto-populate dictionary values associated with frequency contingency table test (CTT) annotations.
    :param freq_cmh: If True, use alternate logic to auto-populate dictionary values associated with frequency Cochran-Mantel-Haenszel (CMH) annotations.
    :param freq_stat_union: If True, use alternate logic to auto-populate dictionary values associated with the union of the contingency table and Cochran-Mantel-Haenszel tests.
    :param description_text: Optional text to append to the end of descriptions. Needs to start with a space if specified.
    :param str age_hist_distribution: Pipe-delimited string of overall age distribution.
    :param sort_order: List containing order to sort label group combinations. Default is SORT_ORDER.
    :return: Dictionary keyed by VCF INFO annotations, where values are dictionaries of Number and Description attributes.
    """
    if prefix != "":
        prefix = f"{prefix}{label_delimiter}"
    if suffix != "":
        suffix = f"{label_delimiter}{suffix}"

    info_dict = dict()

    if age_hist_distribution:
        age_hist_dict = {
            f"{prefix}age_hist_het_bin_freq{suffix}": {
                "Number": "A",
                "Description": (
                    f"Histogram of ages of heterozygous individuals{description_text};"
                    f" bin edges are: {bin_edges['het']}; total number of individuals"
                    f" of any genotype bin: {age_hist_distribution}"
                ),
            },
            f"{prefix}age_hist_het_n_smaller{suffix}": {
                "Number": "A",
                "Description": (
                    "Count of age values falling below lowest histogram bin edge for"
                    f" heterozygous individuals{description_text}"
                ),
            },
            f"{prefix}age_hist_het_n_larger{suffix}": {
                "Number": "A",
                "Description": (
                    "Count of age values falling above highest histogram bin edge for"
                    f" heterozygous individuals{description_text}"
                ),
            },
            f"{prefix}age_hist_hom_bin_freq{suffix}": {
                "Number": "A",
                "Description": (
                    "Histogram of ages of homozygous alternate"
                    f" individuals{description_text}; bin edges are:"
                    f" {bin_edges['hom']}; total number of individuals of any genotype"
                    f" bin: {age_hist_distribution}"
                ),
            },
            f"{prefix}age_hist_hom_n_smaller{suffix}": {
                "Number": "A",
                "Description": (
                    "Count of age values falling below lowest histogram bin edge for"
                    f" homozygous alternate individuals{description_text}"
                ),
            },
            f"{prefix}age_hist_hom_n_larger{suffix}": {
                "Number": "A",
                "Description": (
                    "Count of age values falling above highest histogram bin edge for"
                    f" homozygous alternate individuals{description_text}"
                ),
            },
        }
        info_dict.update(age_hist_dict)

    if popmax:
        popmax_dict = {
            f"{prefix}popmax{suffix}": {
                "Number": "A",
                "Description": (
                    f"Population with the maximum allele frequency{description_text}"
                ),
            },
            f"{prefix}AC{label_delimiter}popmax{suffix}": {
                "Number": "A",
                "Description": (
                    "Allele count in the population with the maximum allele"
                    f" frequency{description_text}"
                ),
            },
            f"{prefix}AN{label_delimiter}popmax{suffix}": {
                "Number": "A",
                "Description": (
                    "Total number of alleles in the population with the maximum allele"
                    f" frequency{description_text}"
                ),
            },
            f"{prefix}AF{label_delimiter}popmax{suffix}": {
                "Number": "A",
                "Description": (
                    f"Maximum allele frequency across populations{description_text}"
                ),
            },
            f"{prefix}nhomalt{label_delimiter}popmax{suffix}": {
                "Number": "A",
                "Description": (
                    "Count of homozygous individuals in the population with the"
                    f" maximum allele frequency{description_text}"
                ),
            },
            f"{prefix}faf95{label_delimiter}popmax{suffix}": {
                "Number": "A",
                "Description": (
                    "Filtering allele frequency (using Poisson 95% CI) for the"
                    f" population with the maximum allele frequency{description_text}"
                ),
            },
        }
        info_dict.update(popmax_dict)
    if grpmax:
        grpmax_dict = {
            f"{prefix}grpmax{suffix}": {
                "Number": "A",
                "Description": (
                    "Genetic ancestry group with the maximum allele"
                    f" frequency{description_text}"
                ),
            },
            f"{prefix}AC{label_delimiter}grpmax{suffix}": {
                "Number": "A",
                "Description": (
                    "Allele count in the genetic ancestry group with the maximum allele"
                    f" frequency{description_text}"
                ),
            },
            f"{prefix}AN{label_delimiter}grpmax{suffix}": {
                "Number": "A",
                "Description": (
                    "Total number of alleles in the genetic ancestry group with the"
                    f" maximum allele frequency{description_text}"
                ),
            },
            f"{prefix}AF{label_delimiter}grpmax{suffix}": {
                "Number": "A",
                "Description": (
                    "Maximum allele frequency across genetic ancestry"
                    f" groups{description_text}"
                ),
            },
            f"{prefix}nhomalt{label_delimiter}grpmax{suffix}": {
                "Number": "A",
                "Description": (
                    "Count of homozygous individuals in the genetic ancestry group"
                    f" with the maximum allele frequency{description_text}"
                ),
            },
        }
        info_dict.update(grpmax_dict)

    if fafmax:
        fafmax_dict = {
            f"{prefix}fafmax{label_delimiter}faf95{label_delimiter}max{suffix}": {
                "Number": "A",
                "Description": (
                    "Maximum filtering allele frequency (using Poisson 95% CI)"
                    f" across genetic ancestry groups{description_text}"
                ),
            },
            f"{prefix}fafmax{label_delimiter}faf95{label_delimiter}max{label_delimiter}gen{label_delimiter}anc{suffix}": {
                "Number": "A",
                "Description": (
                    "Genetic ancestry group with maximum filtering allele"
                    f" frequency (using Poisson 95% CI){description_text}"
                ),
            },
            f"{prefix}fafmax{label_delimiter}faf99{label_delimiter}max{suffix}": {
                "Number": "A",
                "Description": (
                    "Maximum filtering allele frequency (using Poisson 99% CI)"
                    f" across genetic ancestry groups{description_text}"
                ),
            },
            f"{prefix}fafmax{label_delimiter}faf99{label_delimiter}max{label_delimiter}gen{label_delimiter}anc{suffix}": {
                "Number": "A",
                "Description": (
                    "Genetic ancestry group with maximum filtering allele"
                    f" frequency (using Poisson 99% CI){description_text}"
                ),
            },
        }

        info_dict.update(fafmax_dict)

    if callstats or faf or freq_ctt:
        group_types = sorted(label_groups.keys(), key=lambda x: sort_order.index(x))
        combos = make_label_combos(label_groups, label_delimiter=label_delimiter)

        for combo in combos:
            combo_fields = combo.split(label_delimiter)
            group_dict = dict(zip(group_types, combo_fields))

            for_combo = make_combo_header_text("for", group_dict, pop_names)
            in_combo = make_combo_header_text("in", group_dict, pop_names)

            metrics = ["AC", "AN", "AF", "nhomalt", "faf95", "faf99"]
            if freq_ctt:
                metrics += ["CTT_odds_ratio", "CTT_p_value"]
            if prefix_before_metric:
                metric_label_dict = {
                    metric: f"{prefix}{metric}{label_delimiter}{combo}{suffix}"
                    for metric in metrics
                }
            else:
                metric_label_dict = {
                    metric: f"{metric}{label_delimiter}{prefix}{combo}{suffix}"
                    for metric in metrics
                }

            if callstats:
                combo_dict = {
                    metric_label_dict["AC"]: {
                        "Number": "A",
                        "Description": (
                            f"Alternate allele count{for_combo}{description_text}"
                        ),
                    },
                    metric_label_dict["AN"]: {
                        "Number": "1",
                        "Description": (
                            f"Total number of alleles{in_combo}{description_text}"
                        ),
                    },
                    metric_label_dict["AF"]: {
                        "Number": "A",
                        "Description": (
                            f"Alternate allele frequency{in_combo}{description_text}"
                        ),
                    },
                    metric_label_dict["nhomalt"]: {
                        "Number": "A",
                        "Description": (
                            "Count of homozygous"
                            f" individuals{in_combo}{description_text}"
                        ),
                    },
                }
            elif faf:
                if ("XX" in combo_fields) | ("XY" in combo_fields):
                    faf_description_text = (
                        description_text + " in non-PAR regions of sex chromosomes only"
                    )
                else:
                    faf_description_text = description_text
                combo_dict = {
                    metric_label_dict["faf95"]: {
                        "Number": "A",
                        "Description": (
                            "Filtering allele frequency (using Poisson 95%"
                            f" CI){for_combo}{faf_description_text}"
                        ),
                    },
                    metric_label_dict["faf99"]: {
                        "Number": "A",
                        "Description": (
                            "Filtering allele frequency (using Poisson 99%"
                            f" CI){for_combo}{faf_description_text}"
                        ),
                    },
                }
            else:
                combo_dict = {
                    metric_label_dict["CTT_odds_ratio"]: {
                        "Number": "A",
                        "Description": (
                            "Odds ratio from from Hail's contingency_table_test with"
                            " `min_cell_count=100` comparing allele frequencies"
                            f" between exomes and genomes{for_combo}{description_text}"
                        ),
                    },
                    metric_label_dict["CTT_p_value"]: {
                        "Number": "A",
                        "Description": (
                            "P-value from Hail's contingency_table_test with"
                            " `min_cell_count=100` comparing allele frequencies"
                            f" between exomes and genomes{for_combo}{description_text}"
                        ),
                    },
                }
            info_dict.update(combo_dict)
    if freq_cmh:
        cmh_dict = {
            f"{prefix}CMH_chisq{suffix}": {
                "Number": "A",
                "Description": (
                    "Chi-squared test statistic from the Cochran-Mantel-Haenszel test"
                    " comparing allele frequencies between exomes and genomes"
                    f" stratified by genetic ancestry group{description_text}"
                ),
            },
            f"{prefix}CMH_p_value{suffix}": {
                "Number": "A",
                "Description": (
                    "Odds ratio from Cochran-Mantel-Haenszel test comparing allele"
                    " frequencies between exomes and genomes stratified by genetic"
                    f" ancestry group{description_text}"
                ),
            },
        }
        info_dict.update(cmh_dict)
    if freq_stat_union:
        freq_stat_union_dict = {
            f"{prefix}stat_union_p_value{suffix}": {
                "Number": "A",
                "Description": (
                    f"p-value from the contingency table or Cochran-Mantel-Haenszel tests{description_text}"
                ),
            },
            f"{prefix}stat_union_test_name{suffix}": {
                "Number": "A",
                "Description": (
                    f"Name of the test, either contingency_table_test or cochran_mantel_haenszel_test, used to compare allele frequencies between exomes and genomes{description_text}"
                ),
            },
            f"{prefix}stat_union_gen_ancs{suffix}": {
                "Number": ".",
                "Description": (
                    f"List of genetic ancestry groups included in the test. If stat_union_test_name is contingency_table_test, the length of gen_ancs is one and if stat_union_test_name is 'cochran_mantel_haenszel_test', the length of 'gen_ancs' is greater than one{description_text}"
                ),
            },
        }
        info_dict.update(freq_stat_union_dict)

    return info_dict


def add_as_info_dict(
    info_dict: Dict[str, Dict[str, str]] = INFO_DICT, as_fields: List[str] = AS_FIELDS
) -> Dict[str, Dict[str, str]]:
    """
    Update info dictionary with allele-specific terms and their descriptions.

    Used in VCF export.

    :param info_dict: Dictionary containing site-level annotations and their descriptions. Default is INFO_DICT.
    :param as_fields: List containing allele-specific fields to be added to info_dict. Default is AS_FIELDS.
    :return: Dictionary with allele specific annotations, their descriptions, and their VCF number field.
    """
    as_dict = {}
    for field in as_fields:
        try:
            # Strip AS_ from field name
            site_field = field[3:]

            # Get site description from info dictionary and make first letter lower case
            first_letter = info_dict[site_field]["Description"][0].lower()
            rest_of_description = info_dict[site_field]["Description"][1:]

            as_dict[field] = {}
            as_dict[field]["Number"] = "A"
            as_dict[field][
                "Description"
            ] = f"Allele-specific {first_letter}{rest_of_description}"

        except KeyError:
            logger.warning("%s is not present in input info dictionary!", field)

    return as_dict


def make_vcf_filter_dict(
    snp_cutoff: Optional[float] = None,
    indel_cutoff: Optional[float] = None,
    inbreeding_cutoff: Optional[float] = None,
    variant_qc_filter: str = "RF",
    joint: bool = False,
) -> Dict[str, str]:
    """
    Generate dictionary of Number and Description attributes to be used in the VCF header, specifically for FILTER annotations.

    Generates descriptions for:
        - AC0 filter
        - InbreedingCoeff filter
        - Variant QC filter (RF or AS_VQSR)
        - PASS (passed all variant filters)

    :param snp_cutoff: Minimum SNP cutoff score from random forest model.
    :param indel_cutoff: Minimum indel cutoff score from random forest model.
    :param inbreeding_cutoff: Inbreeding coefficient hard cutoff.
    :param variant_qc_filter: Method used for variant QC filter. One of 'RF' or 'AS_VQSR'. Default is 'RF'.
    :param joint: Whether the filter dictionary is for the joint release. Default is False.
    :return: Dictionary keyed by VCF FILTER annotations, where values are Dictionaries of Number and Description attributes.
    """
    variant_qc_filter_dict = {
        "RF": {
            "Description": (
                f"Failed random forest filtering thresholds of {snp_cutoff} for SNPs"
                f" and {indel_cutoff} for indels (probabilities of being a true"
                " positive variant)"
            )
        },
        "AS_VQSR": {
            "Description": (
                f"Failed VQSR filtering thresholds of {snp_cutoff} for SNPs and"
                f" {indel_cutoff} for indels"
            )
        },
    }

    if variant_qc_filter not in variant_qc_filter_dict:
        raise ValueError(
            f"{variant_qc_filter} is not a valid value for 'variant_qc_filter'. It must"
            " be 'RF' or 'AS_VQSR'"
        )
    if not joint and (
        snp_cutoff is None or indel_cutoff is None or inbreeding_cutoff is None
    ):
        raise ValueError(
            "snp_cutoff, indel_cutoff, and inbreeding_cutoff must be specified to generate filter descriptions."
        )

    if joint:
        filter_dict = {
            "PASS": {
                "Description": "Either passed all variant filters in both exomes and "
                "genomes, or passed all variant filters in either "
                "exomes or genomes while being absent from the other "
                "dataset"
            },
            "EXOMES_FILTERED": {
                "Description": "Failed variant filters in the exomes dataset and either "
                "passed all variant filters in the genomes dataset or the variant was "
                "not present in the genomes dataset. Refer to 'exomes_filters' within "
                "INFO for more information"
            },
            "GENOMES_FILTERED": {
                "Description": "Failed variant filters in the genomes dataset and either "
                "passed all variant filters in the exomes dataset or the variant was "
                "not present in the exomes dataset. Refer to 'genomes_filters' within "
                "INFO for more information"
            },
            "BOTH_FILTERED": {
                "Description": "Failed variant filters in both exomes and genomes datasets. "
                "Refer to 'exomes_filters' and 'genomes_filters' within INFO for more information"
            },
        }
    else:
        filter_dict = {
            "AC0": {
                "Description": (
                    "Allele count is zero after filtering out low-confidence genotypes (GQ"
                    " < 20; DP < 10; and AB < 0.2 for het calls)"
                )
            },
            "InbreedingCoeff": {
                "Description": f"Inbreeding coefficient < {inbreeding_cutoff}"
            },
            "PASS": {"Description": "Passed all variant filters"},
            variant_qc_filter: variant_qc_filter_dict[variant_qc_filter],
        }

    return filter_dict


def make_hist_bin_edges_expr(
    ht: hl.Table,
    hists: List[str] = HISTS,
    ann_with_hists: Optional[str] = None,
    prefix: str = "",
    label_delimiter: str = "_",
    include_age_hists: bool = True,
) -> Dict[str, str]:
    """
    Create dictionaries containing variant histogram annotations and their associated bin edges, formatted into a string separated by pipe delimiters.

    :param ht: Table containing histogram variant annotations.
    :param hists: List of variant histogram annotations. Default is HISTS.
    :param ann_with_hists: Name of row annotation containing histogram data. In exomes or
        genomes release HT, `histograms` is a row, but in the joint release HT, it's
        under the row of `exomes`, `genomes`, or `joint`.
    :param prefix: Prefix text for age histogram bin edges.  Default is empty string.
    :param label_delimiter: String used as delimiter between prefix and histogram annotation.
    :param include_age_hists: Include age histogram annotations.
    :return: Dictionary keyed by histogram annotation name, with corresponding
        reformatted bin edges for values.
    """
    # Add underscore to prefix if it isn't empty
    if prefix:
        prefix += label_delimiter

    edges_dict = {}

    if include_age_hists:
        for call_type in ["het", "hom"]:
            if ann_with_hists:
                bin_edges = (
                    ht.filter(
                        hl.is_defined(
                            ht[ann_with_hists]
                            .histograms.age_hists[f"age_hist_{call_type}"]
                            .bin_edges
                        )
                    )[ann_with_hists]
                    .histograms.age_hists[f"age_hist_{call_type}"]
                    .bin_edges.take(1)[0]
                )
            else:
                bin_edges = (
                    ht.filter(
                        hl.is_defined(
                            ht.histograms.age_hists[f"age_hist_{call_type}"].bin_edges
                        )
                    )
                    .histograms.age_hists[f"age_hist_{call_type}"]
                    .bin_edges.take(1)[0]
                )

            if bin_edges:
                edges_dict[f"{prefix}{call_type}"] = "|".join(
                    map(lambda x: f"{x:.1f}", bin_edges)
                )

    for hist in hists:
        # Parse hists calculated on both raw and adj-filtered data
        for hist_type in [f"{prefix}raw_qual_hists", f"{prefix}qual_hists"]:
            hist_name = hist if "raw" not in hist_type else f"{prefix}{hist}_raw"

            if ann_with_hists:
                bin_edges = (
                    ht.filter(
                        hl.is_defined(
                            ht[ann_with_hists].histograms[hist_type][hist].bin_edges
                        )
                    )[ann_with_hists]
                    .histograms[hist_type][hist]
                    .bin_edges.take(1)[0]
                )
            else:
                bin_edges = (
                    ht.filter(hl.is_defined(ht.histograms[hist_type][hist].bin_edges))
                    .histograms[hist_type][hist]
                    .bin_edges.take(1)[0]
                )
            if bin_edges:
                edges_dict[hist_name] = "|".join(
                    map(
                        lambda x: f"{x:.2f}" if "ab" in hist else str(int(x)), bin_edges
                    )
                )

    return edges_dict


def make_hist_dict(
    bin_edges: Dict[str, Dict[str, str]],
    adj: bool,
    hist_metric_list: List[str] = HISTS,
    label_delimiter: str = "_",
    drop_n_smaller_larger: bool = False,
    prefix: str = "",
    suffix: str = "",
    description_text: str = "",
) -> Dict[str, str]:
    """
    Generate dictionary of Number and Description attributes to be used in the VCF header, specifically for histogram annotations.

    :param bin_edges: Dictionary keyed by histogram annotation name, with corresponding string-reformatted bin edges for values.
    :param adj: Whether to create a header dict for raw or adj quality histograms.
    :param hist_metric_list: List of hists for which to build hist info dict
    :param label_delimiter: String used as delimiter in values stored in hist_metric_list.
    :param drop_n_smaller_larger: Whether to drop n_smaller and n_larger annotations from header dict. Default is False.
    :param prefix: Prefix text for histogram annotations. Default is empty string.
    :param suffix: Suffix text for histogram annotations. Default is empty string.
    :param description_text: Optional text to append to the end of descriptions. Needs to start with a space if specified.
    :return: Dictionary keyed by VCF INFO annotations, where values are Dictionaries of Number and Description attributes.
    """
    if prefix != "":
        prefix = f"{prefix}{label_delimiter}"
    if suffix != "":
        suffix = f"{label_delimiter}{suffix}"

    header_hist_dict = {}
    for hist in hist_metric_list:
        # Get hists for both raw and adj data
        # Add "_raw" to quality histograms calculated on raw data
        if not adj:
            hist = f"{hist}_raw"

        edges = bin_edges[hist]
        hist_fields = hist.split(label_delimiter)
        hist_text = hist_fields[0].upper()

        if hist_fields[2] == "alt":
            hist_text = hist_text + " in heterozygous individuals"
        if adj:
            hist_text = hist_text + " calculated on high quality genotypes"

        hist_dict = {
            f"{prefix}{hist}_bin_freq{suffix}": {
                "Number": "A",
                "Description": (
                    f"Histogram for {hist_text}{description_text}; bin edges are:"
                    f" {edges}"
                ),
            },
        }
        # These annotations are frequently zero and are dropped from gnomad
        # releases for most histograms.
        if not drop_n_smaller_larger:
            hist_dict.update(
                {
                    f"{prefix}{hist}_n_smaller{suffix}": {
                        "Number": "A",
                        "Description": (
                            f"Count of {hist_fields[0].upper()} values falling below"
                            f" lowest histogram bin edge {hist_text}{description_text}"
                        ),
                    },
                    f"{prefix}{hist}_n_larger{suffix}": {
                        "Number": "A",
                        "Description": (
                            f"Count of {hist_fields[0].upper()} values falling above"
                            f" highest histogram bin edge {hist_text}{description_text}"
                        ),
                    },
                }
            )
        # Only add n_larger for dp qual histograms.
        if "dp" in hist:
            hist_dict.update(
                {
                    f"{prefix}{hist}_n_larger{suffix}": {
                        "Number": "A",
                        "Description": (
                            f"Count of {hist_fields[0].upper()} values falling above"
                            f" highest histogram bin edge {hist_text}{description_text}"
                        ),
                    },
                }
            )

        header_hist_dict.update(hist_dict)

    return header_hist_dict


def set_female_y_metrics_to_na(
    t: Union[hl.Table, hl.MatrixTable],
) -> Dict[str, hl.expr.Int32Expression]:
    """
    Set AC, AN, and nhomalt chrY variant annotations for females to NA (instead of 0).

    :param t: Table/MatrixTable containing female variant annotations.
    :return: Dictionary with reset annotations
    """
    metrics = list(t.row.info)
    female_metrics = [x for x in metrics if "_female" in x or "_XX" in x]

    female_metrics_dict = {}
    for metric in female_metrics:
        female_metrics_dict.update(
            {
                f"{metric}": hl.or_missing(
                    (~t.locus.in_y_nonpar() & ~t.locus.in_y_par()),
                    t.info[f"{metric}"],
                )
            }
        )
    return female_metrics_dict


def build_vcf_export_reference(
    name: str,
    build: str = "GRCh38",
    keep_contigs: List[str] = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"],
    keep_chrM: bool = True,
) -> hl.ReferenceGenome:
    """
    Create export reference based on reference genome defined by `build`.

    By default this will return a new reference with all non-standard contigs eliminated. Keeps chr 1-22, Y, X, and M.

    An example of a non-standard contig is: ##contig=<ID=chr3_GL000221v1_random,length=155397,assembly=GRCh38>

    :param name: Name to use for new reference.
    :param build: Reference genome build to use as starting reference genome.
    :param keep_contigs: Contigs to keep from reference genome defined by `build`. Default is autosomes and sex chromosomes.
    :param keep_chrM: Whether to keep chrM. Default is True.
    :return: Reference genome for VCF export containing only contigs in `keep_contigs`.
    """
    ref = hl.get_reference(build)
    ref_args = {}

    if keep_chrM:
        keep_contigs.extend(ref.mt_contigs)
        ref_args.update({"mt_contigs": ref.mt_contigs})

    ref_args.update(
        {
            "name": name,
            "contigs": keep_contigs,
            "lengths": {contig: ref.lengths[contig] for contig in keep_contigs},
            "x_contigs": ref.x_contigs,
            "y_contigs": ref.y_contigs,
            "par": [
                (interval.start.contig, interval.start.position, interval.end.position)
                for interval in ref.par
            ],
        }
    )

    export_reference = hl.ReferenceGenome(**ref_args)

    return export_reference


def rekey_new_reference(
    t: Union[hl.Table, hl.MatrixTable], reference: hl.ReferenceGenome
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Re-key Table or MatrixTable with a new reference genome.

    :param t: Input Table/MatrixTable.
    :param reference: Reference genome to re-key with.
    :return: Re-keyed Table/MatrixTable
    """
    t = t.rename({"locus": "locus_original"})
    locus_expr = hl.locus(
        t.locus_original.contig,
        t.locus_original.position,
        reference_genome=reference,
    )

    if isinstance(t, hl.MatrixTable):
        t = t.annotate_rows(locus=locus_expr)
        t = t.key_rows_by("locus", "alleles").drop("locus_original")
    else:
        t = t.annotate(locus=locus_expr)
        t = t.key_by("locus", "alleles").drop("locus_original")

    return t
