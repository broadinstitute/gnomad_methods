import hail as hl
from typing import List
import logging

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


GROUPS = ["adj", "raw"]
"""
Group names used to generate labels for high quality genotypes and all raw genotypes. Used in VCF export.
"""

HISTS = ["gq_hist_alt", "gq_hist_all", "dp_hist_alt", "dp_hist_all", "ab_hist_alt"]
"""
Quality histograms used in VCF export.
"""

POPS = ["afr", "amr", "asj", "eas", "fin", "nfe", "oth", "sas"]
"""
Continental population names used in VCF export.
"""

NFE_SUBPOPS = ["onf", "bgr", "swe", "nwe", "seu", "est"]
"""
gnomAD subpopulations for NFE population. Used in VCF export.
"""

EAS_SUBPOPS = ["kor", "oea", "jpn"]
"""
gnomAD subpopulations for EAS population. Used in VCF export.
"""

FAF_POPS = ["afr", "amr", "eas", "nfe", "sas"]
"""
Populations that are included in filtering allele frequency calculations. Used in VCF export.
"""

SEXES = ["male", "female"]
"""
Sample sexes used in VCF export.
"""

AS_FIELDS = [
    "AS_BaseQRankSum",
    "AS_FS",
    "AS_MQ",
    "AS_MQRankSum",
    "AS_pab_max",
    "AS_QD",
    "AS_ReadPosRankSum",
    "AS_SOR",
    "AS_VarDP",
]
"""
Allele-specific variant annotations.
"""

SITE_FIELDS = [
    "FS",
    "InbreedingCoeff",
    "MQ",
    "MQRankSum",
    "QD",
    "ReadPosRankSum",
    "sibling_singleton",
    "SOR",
    "transmitted_singleton",
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

REGION_TYPE_FIELDS = ["lcr", "nonpar"]
"""
Annotations about variant region type.
"""

RF_FIELDS = [
    "fail_hard_filters",
    "filters",
    "rf_label",
    "rf_train",
    "rf_probability",
    "tp",
]
"""
Annotations specific to the random forest model.
"""

VQSR_FIELDS = ["AS_VQSLOD", "AS_culprit", "NEGATIVE_TRAIN_SITE", "POSITIVE_TRAIN_SITE"]
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
        "Description": "Inbreeding coefficient, the excess heterozygosity at a variant site, computed as 1 - (the number of heterozygous genotypes)/(the number of heterozygous genotypes expected under Hardy-Weinberg equilibrium)"
    },
    "MQ": {
        "Description": "Root mean square of the mapping quality of reads across all samples"
    },
    "MQRankSum": {
        "Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities"
    },
    "QD": {
        "Description": "Variant call confidence normalized by depth of sample reads supporting a variant"
    },
    "ReadPosRankSum": {
        "Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference read position bias"
    },
    "SOR": {"Description": "Strand bias estimated by the symmetric odds ratio test"},
    "VQSR_POSITIVE_TRAIN_SITE": {
        "Description": "Variant was used to build the positive training set of high-quality variants for VQSR"
    },
    "VQSR_NEGATIVE_TRAIN_SITE": {
        "Description": "Variant was used to build the negative training set of low-quality variants for VQSR"
    },
    "BaseQRankSum": {
        "Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference base qualities",
    },
    "VarDP": {
        "Description": "Depth over variant genotypes (does not include depth of reference samples)"
    },
    "AS_VQSLOD": {
        "Number": "A",
        "Description": "Allele-specific log-odds ratio of being a true variant versus being a false positive under the trained VQSR Gaussian mixture model",
    },
    "AS_VQSR_culprit": {
        "Number": "A",
        "Description": "Allele-specific worst-performing annotation in the VQSR Gaussian mixture model",
    },
    "lcr": {"Description": "Variant falls within a low complexity region"},
    "nonpar": {
        "Description": "Variant (on sex chromosome) falls outside a pseudoautosomal region"
    },
    "rf_positive_label": {
        "Description": "Variant was labelled as a positive example for training of random forest model"
    },
    "rf_negative_label": {
        "Description": "Variant was labelled as a negative example for training of random forest model"
    },
    "rf_label": {"Description": "Random forest training label"},
    "rf_train": {"Description": "Variant was used in training random forest model"},
    "rf_tp_probability": {
        "Description": "Probability of a called variant being a true variant as determined by random forest model"
    },
    "transmitted_singleton": {
        "Description": "Variant was a callset-wide doubleton that was transmitted within a family from a parent to a child (i.e., a singleton amongst unrelated samples in cohort)"
    },
    "original_alleles": {"Description": "Alleles before splitting multiallelics"},
    "variant_type": {
        "Description": "Variant type (snv, indel, multi-snv, multi-indel, or mixed)"
    },
    "allele_type": {
        "Number": "A",
        "Description": "Allele type (snv, insertion, deletion, or mixed)",
    },
    "n_alt_alleles": {
        "Number": "A",
        "Description": "Total number of alternate alleles observed at variant locus",
    },
    "was_mixed": {"Description": "Variant type was mixed"},
    "has_star": {
        "Description": "Variant locus coincides with a spanning deletion (represented by a star) observed elsewhere in the callset"
    },
    "pab_max": {
        "Number": "A",
        "Description": "Maximum p-value over callset for binomial test of observed allele balance for a heterozygous genotype, given expectation of 0.5",
    },
}
"""
Dictionary used during VCF export to export row (variant) annotations.
"""

ENTRIES = ["GT", "GQ", "DP", "AD", "MIN_DP", "PGT", "PID", "PL", "SB"]
"""
Densified entries to be selected during VCF export.
"""

SPARSE_ENTRIES = [
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
        "Description": "Approximate read depth (reads with MQ=255 or with bad mates are filtered)",
        "Number": "1",
        "Type": "Integer",
    },
    "GQ": {
        "Description": "Phred-scaled confidence that the genotype assignment is correct. Value is the difference between the second lowest PL and the lowest PL (always normalized to 0).",
        "Number": "1",
        "Type": "Integer",
    },
    "MIN_DP": {
        "Description": "Minimum DP observed within the GVCF block",
        "Number": "1",
        "Type": "Integer",
    },
    "PGT": {
        "Description": "Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another",
        "Number": "1",
        "Type": "String",
    },
    "PID": {
        "Description": "Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group",
        "Number": "1",
        "Type": "String",
    },
    "PL": {
        "Description": "Normalized, phred-scaled likelihoods for genotypes as defined in the VCF specification",
        "Number": "G",
        "Type": "Integer",
    },
    "SB": {
        "Description": "Per-sample component statistics which comprise the Fisher's exact test to detect strand bias. Values are: depth of reference allele on forward strand, depth of reference allele on reverse strand, depth of alternate allele on forward strand, depth of alternate allele on reverse strand.",
        "Number": "4",
        "Type": "Integer",
    },
}
"""
Dictionary used during VCF export to export MatrixTable entries.
"""


def ht_to_vcf_mt(
    info_ht: hl.Table,
    pipe_delimited_annotations: List[str] = INFO_VCF_AS_PIPE_DELIMITED_FIELDS,
) -> hl.MatrixTable:
    """
    Creates a MT ready for vcf export from a HT. In particular, the following conversions are done:
    - All int64 are coerced to int32
    - Fields specified by `pipe_delimited_annotations` will be converted from arrays to pipe-delimited strings

    .. note::

        The MT returned has no cols.

    :param info_ht: Input HT
    :param pipe_delimited_annotations: List of info fields (they must be fields of the ht.info Struct)
    :return: MatrixTable ready for VCF export
    """

    def get_pipe_expr(array_expr: hl.expr.ArrayExpression) -> hl.expr.StringExpression:
        return hl.delimit(array_expr.map(lambda x: hl.or_else(hl.str(x), "")), "|")

    # Make sure the HT is keyed by locus, alleles
    info_ht = info_ht.key_by("locus", "alleles")

    # Convert int64 fields to int32 (int64 isn't supported by VCF)
    for f, ft in info_ht.info.dtype.items():
        if ft == hl.dtype("int64"):
            logger.warning(
                f"Coercing field info.{f} from int64 to int32 for VCF output. Value will be capped at int32 max value."
            )
            info_ht = info_ht.annotate(
                info=info_ht.info.annotate(
                    **{f: hl.int32(hl.min(2 ** 31 - 1, info_ht.info[f]))}
                )
            )
        elif ft == hl.dtype("array<int64>"):
            logger.warning(
                f"Coercing field info.{f} from array<int64> to array<int32> for VCF output. Array values will be capped at int32 max value."
            )
            info_ht = info_ht.annotate(
                info=info_ht.info.annotate(
                    **{
                        f: info_ht.info[f].map(
                            lambda x: hl.int32(hl.min(2 ** 31 - 1, x))
                        )
                    }
                )
            )

    info_expr = {}

    # Make sure to pipe-delimit fields that need to.
    # Note: the expr needs to be prefixed by "|" because GATK expect one value for the ref (always empty)
    # Note2: this doesn't produce the correct annotation for AS_SB_TABLE, but it is overwritten below
    for f in pipe_delimited_annotations:
        if f in info_ht.info:
            info_expr[f] = "|" + get_pipe_expr(info_ht.info[f])

    # Flatten SB if it is an array of arrays
    if "SB" in info_ht.info and not isinstance(
        info_ht.info.SB, hl.expr.ArrayNumericExpression
    ):
        info_expr["SB"] = info_ht.info.SB[0].extend(info_ht.info.SB[1])

    if "AS_SB_TABLE" in info_ht.info:
        info_expr["AS_SB_TABLE"] = get_pipe_expr(
            info_ht.info.AS_SB_TABLE.map(lambda x: hl.delimit(x, ","))
        )

    # Annotate with new expression and add 's' empty string field required to cast HT to MT
    info_ht = info_ht.annotate(
        info=info_ht.info.annotate(**info_expr), s=hl.null(hl.tstr)
    )

    # Create an MT with no cols so that we acn export to VCF
    info_mt = info_ht.to_matrix_table_row_major(columns=["s"], entry_field_name="s")
    return info_mt.filter_cols(False)


def make_label_combos(label_groups: Dict[str, List[str]]) -> List[str]:
    """
    Make combinations of all possible labels for a supplied dictionary of label groups
    :param dict label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"])
    :return: list of all possible combinations of values for the supplied label groupings
    :rtype: list[str]
    """
    copy_label_groups = copy.deepcopy(label_groups)
    if len(copy_label_groups) == 1:
        return [item for sublist in copy_label_groups.values() for item in sublist]
    anchor_group = sorted(copy_label_groups.keys(), key=lambda x: SORT_ORDER.index(x))[
        0
    ]
    anchor_val = copy_label_groups.pop(anchor_group)
    combos = []
    for x, y in itertools.product(anchor_val, make_label_combos(copy_label_groups)):
        combos.append("{0}_{1}".format(x, y))
    return combos


def generic_field_check(
    ht: hl.Table, cond_expr, check_description, display_fields, verbose
):
    """
    Check a generic logical condition involving annotations in a Hail Table and print the results to terminal
    :param Table ht: Table containing annotations to be checked
    :param Expression cond_expr: logical expression referring to annotations in ht to be checked
    :param str check_description: String describing the condition being checked; is displayed in terminal summary message
    :param list of str display_fields: List of names of ht annotations to be displayed in case of failure (for troubleshooting purposes);
        these fields are also displayed if verbose is True
    :param bool verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks
    :rtype: None
    """
    ht_orig = ht
    ht = ht.filter(cond_expr)
    n_fail = ht.count()
    if n_fail > 0:
        logger.info(f"Found {n_fail} sites that fail {check_description} check:")
        ht = ht.flatten()
        ht.select("locus", "alleles", *display_fields).show()
    else:
        logger.info(f"PASSED {check_description} check")
        if verbose:
            ht_orig = ht_orig.flatten()
            ht_orig.select(*display_fields).show()


def make_filters_sanity_check_expr(ht: hl.Table) -> Dict[str, hl.expr.Expression]:
    """
    Make Hail Expressions to measure % variants filtered under varying conditions of interest
    :param Table ht: Hail Table containing 'filter' annotation to be examined
    :return: Dictionary containing Hail aggregation expressions to measure filter flags
    :rtype: Dict of str: Expression
    """
    filters_dict = {
        "n": hl.agg.count(),
        "frac_any_filter": hl.agg.count_where(ht.is_filtered) / hl.agg.count(),
        "frac_inbreed_coeff": hl.agg.count_where(
            ht.filters.contains("inbreeding_coeff")
        )
        / hl.agg.count(),
        "frac_ac0": hl.agg.count_where(ht.filters.contains("AC0")) / hl.agg.count(),
        "frac_rf": hl.agg.count_where(ht.filters.contains("rf")) / hl.agg.count(),
        "frac_inbreed_coeff_only": hl.agg.count_where(
            ht.filters.contains("inbreeding_coeff") & (ht.filters.length() == 1)
        )
        / hl.agg.count(),
        "frac_ac0_only": hl.agg.count_where(
            ht.filters.contains("AC0") & (ht.filters.length() == 1)
        )
        / hl.agg.count(),
        "frac_rf_only": hl.agg.count_where(
            ht.filters.contains("rf") & (ht.filters.length() == 1)
        )
        / hl.agg.count(),
    }
    return filters_dict


def add_as_info_dict(INFO_DICT: Dict[str, Dict[str, str]]) -> Dict[str, Dict[str, str]]:
    """
    Updates info dictionary with allele-specific terms and their descriptions.

    Used in VCF export.

    :param INFO_DICT: Dictionary containing site-level annotations and their descriptions.
    :return: Dictionary with allele specific annotations, their descriptions, and their VCF number field.
    """
    prefix_text = "Allele-specific"
    AS_DICT = {}

    for field in AS_FIELDS:
        # Strip AS_ from field name
        site_field = field[3:]

        AS_DICT[field] = {}
        AS_DICT[field]["Number"] = "A"
        AS_DICT[field][
            "Description"
        ] = f"{prefix_text} {INFO_DICT[site_field]['Description'][0].lower()}{INFO_DICT[site_field]['Description'][1:]}"

    return AS_DICT


def make_vcf_filter_dict(
    ht: hl.Table, snp_cutoff: float, indel_cutoff: float, inbreeding_cutoff: float
) -> Dict[str, str]:
    """
    Generates dictionary of Number and Description attributes to be used in the VCF header, specifically for FILTER annotations.

    :param Table ht: Table containing global annotations of the Random Forests SNP and indel cutoffs.
    :param float snp_cutoff: Minimum SNP cutoff score from random forest model.
    :param float indel_cutoff: Minimum indel cutoff score from random forest model.
    :param float inbreeding_cutoff: Inbreeding coefficient hard cutoff.
    :return: Dictionary keyed by VCF FILTER annotations, where values are Dictionaries of Number and Description attributes.
    :rtype: Dict[str, str]
    """
    filter_dict = {
        "AC0": {
            "Description": "Allele count is zero after filtering out low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)"
        },
        "InbreedingCoeff": {"Description": f"InbreedingCoeff < {inbreeding_cutoff}"},
        "RF": {
            "Description": f"Failed random forest filtering thresholds of {snp_cutoff}, {indel_cutoff} (probabilities of being a true positive variant) for SNPs, indels"
        },
        "PASS": {"Description": "Passed all variant filters"},
    }
    return filter_dict
