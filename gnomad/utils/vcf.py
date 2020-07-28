import copy
import itertools
import logging
from typing import Dict, List, Union

import hail as hl

from gnomad.sample_qc.ancestry import POP_NAMES

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

SORT_ORDER = ["popmax", "group", "pop", "subpop", "sex"]
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

FAF_POPS = ["afr", "amr", "eas", "nfe", "sas"]
"""
Global populations that are included in filtering allele frequency (faf) calculations. Used in VCF export.
"""

SEXES = ["male", "female"]
"""
Sample sexes used in VCF export.
Used to stratify frequency annotations (AC, AN, AF) for each sex.
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

REGION_FLAG_FIELDS = ["decoy", "lcr", "nonpar", "segdup"]
"""
Annotations about variant region type.

.. note::
    decoy and segdup resource files do not currently exist for GrCh38/hg38.
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
Annotations specific to the variant QC using a random forest model.
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
    "decoy": {"Description": "Variant falls within a reference decoy region"},
    "lcr": {"Description": "Variant falls within a low complexity region"},
    "nonpar": {
        "Description": "Variant (on sex chromosome) falls outside a pseudoautosomal region"
    },
    "segdup": {"Description": "Variant falls within a segmental duplication region"},
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


def make_label_combos(
    label_groups: Dict[str, List[str]], sort_order: List[str] = SORT_ORDER,
) -> List[str]:
    """
    Make combinations of all possible labels for a supplied dictionary of label groups.

    For example, if label_groups is `{"sex": ["male", "female"], "pop": ["afr", "nfe", "amr"]}`,
    this function will return `["afr_male", "afr_female", "nfe_male", "nfe_female", "amr_male", "amr_female']`

    :param label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"]).
    :param sort_order: List containing order to sort label group combinations. Default is SORT_ORDER.
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
    for x, y in itertools.product(anchor_val, make_label_combos(copy_label_groups)):
        combos.append("{0}_{1}".format(x, y))
    return combos


def generic_field_check(
    ht: hl.Table,
    cond_expr: hl.expr.BooleanExpression,
    check_description: str,
    display_fields: List[str],
    verbose: bool,
) -> None:
    """
    Check a generic logical condition involving annotations in a Hail Table and print the results to terminal.

    Displays the number of rows in the Table that return False for a condition.

    :param ht: Table containing annotations to be checked.
    :param cond_expr: Logical expression referring to annotations in ht to be checked.
    :param check_description: String describing the condition being checked; is displayed in terminal summary message.
    :param display_fields: List of names of ht annotations to be displayed in case of failure (for troubleshooting purposes);
        these fields are also displayed if verbose is True.
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
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
    Make Hail expressions to measure % variants filtered under varying conditions of interest.

    Checks for:
        - Total number of variants
        - Fraction of variants removed by any filter
        - Fraction of variants removed because of InbreedingCoefficient filter in combination with any other filter
        - Fraction of varinats removed because of AC0 filter in combination with any other filter
        - Fraction of variants removed because of random forest filtering in combination with any other filter
        - Fraction of variants removed only because of InbreedingCoefficient filter
        - Fraction of variants removed only because of AC0 filter
        - Fraction of variants removed only because of random forest filtering


    :param ht: Table containing 'filter' annotation to be examined.
    :return: Dictionary containing Hail aggregation expressions to measure filter flags.
    """
    filters_dict = {
        "n": hl.agg.count(),
        "frac_any_filter": hl.agg.fraction(hl.len(ht.filters) == 0),
        "frac_inbreed_coeff": hl.agg.fraction(ht.filters.contains("inbreeding_coeff")),
        "frac_ac0": hl.agg.fraction(ht.filters.contains("AC0")),
        "frac_rf": hl.agg.fraction(ht.filters.contains("rf")),
        "frac_inbreed_coeff_only": hl.agg.fraction(
            ht.filters.contains("inbreeding_coeff") & (ht.filters.length() == 1)
        ),
        "frac_ac0_only": hl.agg.fraction(
            ht.filters.contains("AC0") & (ht.filters.length() == 1)
        ),
        "frac_rf_only": hl.agg.fraction(
            ht.filters.contains("rf") & (ht.filters.length() == 1)
        ),
    }
    return filters_dict


def make_combo_header_text(
    preposition: str, group_types: List[str], combo_fields: List[str], prefix: str,
) -> str:
    """
    Programmatically generate text to populate the VCF header description for a given variant annotation with specific groupings and subset.

    For example, if preposition is "for", group_types is ["group", "pop", "sex"], and combo_fields is ["adj", "afr", "female"],
    this function will return the string " for female samples of African-American/African ancestry".

    :param preposition: Relevant preposition to precede automatically generated text.
    :param group_types: List of grouping types, e.g. "sex" or "pop".
    :param combo_fields: List of the specific values for each grouping type, for which the text is being generated.
    :param prefix: Prefix string indicating sample subset.
    :return: String with automatically generated description text for a given set of combo fields.
    """
    combo_dict = dict(zip(group_types, combo_fields))
    header_text = " " + preposition

    if len(combo_dict.keys()) == 1:
        if combo_dict["group"] == "adj":
            return ""

    if "sex" in combo_dict.keys():
        header_text = header_text + " " + combo_dict["sex"]

    header_text = header_text + " samples"

    if "subpop" in combo_dict.keys():
        header_text = header_text + f" of {POP_NAMES[combo_dict['subpop']]} ancestry"
        combo_dict.pop("pop")

    if "pop" in combo_dict.keys():
        header_text = header_text + f" of {POP_NAMES[combo_dict['pop']]} ancestry"

    if "gnomad" in prefix:
        header_text = header_text + " in gnomAD"

    if "group" in group_types:
        if combo_dict["group"] == "raw":
            header_text = header_text + ", before removing low-confidence genotypes"

    return header_text


def make_info_dict(
    prefix: str,
    label_groups: Dict[str, str] = None,
    bin_edges: Dict[str, str] = None,
    faf: bool = False,
    popmax: bool = False,
    description_text: str = "",
    age_hist_data: str = None,
) -> Dict[str, Dict[str, str]]:
    """
    Generate dictionary of Number and Description attributes of VCF INFO fields.

    Used to populate the VCF header during export.
    
    :param prefix: gnomAD or empty string.
    :param label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"]).
    :param bin_edges: Dictionary keyed by annotation type, with values that reflect the bin edges corresponding to the annotation.
    :param faf: If True, use alternate logic to auto-populate dictionary values associated with filter allele frequency annotations.
    :param popmax: If True, use alternate logic to auto-populate dictionary values associated with popmax annotations.
    :param description_text: Optinal text to append to the end of age and popmax descriptions. 
    :param str age_hist_data: Pipe-delimited string of age histograms, from `get_age_distributions`.
    :return: Dictionary keyed by VCF INFO annotations, where values are dictionaries of Number and Description attributes.
    """
    if prefix != "":
        prefix = f"{prefix}_"

    info_dict = dict()

    if age_hist_data:
        age_hist_dict = {
            f"{prefix}age_hist_het_bin_freq": {
                "Number": "A",
                "Description": f"Histogram of ages of heterozygous individuals{description_text}; bin edges are: {bin_edges['het']}; total number of individuals of any genotype bin: {age_hist_data}",
            },
            f"{prefix}age_hist_het_n_smaller": {
                "Number": "A",
                "Description": f"Count of age values falling below lowest histogram bin edge for heterozygous individuals{description_text}",
            },
            f"{prefix}age_hist_het_n_larger": {
                "Number": "A",
                "Description": f"Count of age values falling above highest histogram bin edge for heterozygous individuals{description_text}",
            },
            f"{prefix}age_hist_hom_bin_freq": {
                "Number": "A",
                "Description": f"Histogram of ages of homozygous alternate individuals{description_text}; bin edges are: {bin_edges['het']}; total number of individuals of any genotype bin: {age_hist_data}",
            },
            f"{prefix}age_hist_hom_n_smaller": {
                "Number": "A",
                "Description": f"Count of age values falling below lowest histogram bin edge for homozygous alternate individuals{description_text}",
            },
            f"{prefix}age_hist_hom_n_larger": {
                "Number": "A",
                "Description": f"Count of age values falling above highest histogram bin edge for homozygous alternate individuals{description_text}",
            },
        }
        info_dict.update(age_hist_dict)

    if popmax:
        popmax_dict = {
            f"{prefix}popmax": {
                "Number": "A",
                "Description": f"Population with maximum AF{description_text}",
            },
            f"{prefix}AC_popmax": {
                "Number": "A",
                "Description": f"Allele count in the population with the maximum AF{description_text}",
            },
            f"{prefix}AN_popmax": {
                "Number": "A",
                "Description": f"Total number of alleles in the population with the maximum AF{description_text}",
            },
            f"{prefix}AF_popmax": {
                "Number": "A",
                "Description": f"Maximum allele frequency across populations (excluding samples of Ashkenazi, Finnish, and indeterminate ancestry){description_text}",
            },
            f"{prefix}nhomalt_popmax": {
                "Number": "A",
                "Description": f"Count of homozygous individuals in the population with the maximum allele frequency{description_text}",
            },
        }
        info_dict.update(popmax_dict)

    else:
        group_types = sorted(label_groups.keys(), key=lambda x: SORT_ORDER.index(x))
        combos = make_label_combos(label_groups)

        for combo in combos:
            combo_fields = combo.split("_")

            if not faf:
                combo_dict = {
                    f"{prefix}AC_{combo}": {
                        "Number": "A",
                        "Description": f"Alternate allele count{make_combo_header_text('for', group_types, combo_fields, prefix)}",
                    },
                    f"{prefix}AN_{combo}": {
                        "Number": "1",
                        "Description": f"Total number of alleles{make_combo_header_text('in', group_types, combo_fields, prefix)}",
                    },
                    f"{prefix}AF_{combo}": {
                        "Number": "A",
                        "Description": f"Alternate allele frequency{make_combo_header_text('in', group_types, combo_fields, prefix)}",
                    },
                    f"{prefix}nhomalt_{combo}": {
                        "Number": "A",
                        "Description": f"Count of homozygous individuals{make_combo_header_text('in', group_types, combo_fields, prefix)}",
                    },
                }
            else:
                combo_dict = {
                    f"{prefix}faf95_{combo}": {
                        "Number": "A",
                        "Description": f"Filtering allele frequency (using Poisson 95% CI) {make_combo_header_text('for', group_types, combo_fields, prefix)}",
                    },
                    f"{prefix}faf99_{combo}": {
                        "Number": "A",
                        "Description": f"Filtering allele frequency (using Poisson 99% CI) {make_combo_header_text('for', group_types, combo_fields, prefix)}",
                    },
                }
            info_dict.update(combo_dict)
    return info_dict


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
    snp_cutoff: float, indel_cutoff: float, inbreeding_cutoff: float
) -> Dict[str, str]:
    """
    Generates dictionary of Number and Description attributes to be used in the VCF header, specifically for FILTER annotations.

    :param ht: Table containing global annotations of the Random Forests SNP and indel cutoffs.
    :param snp_cutoff: Minimum SNP cutoff score from random forest model.
    :param indel_cutoff: Minimum indel cutoff score from random forest model.
    :param inbreeding_cutoff: Inbreeding coefficient hard cutoff.
    :return: Dictionary keyed by VCF FILTER annotations, where values are Dictionaries of Number and Description attributes.
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


def make_hist_bin_edges_expr(ht: hl.Table, prefix: str) -> Dict[str, str]:
    """
    Create dictionaries containing variant histogram annotations and their associated bin edges, formatted into a string
    separated by pipe delimiters.

    :param ht: Table containing histogram variant annotations.
    :param prefix: Prefix text for age histogram bin edges. 
    :return: Dictionary keyed by histogram annotation name, with corresponding reformatted bin edges for values.
    """
    edges_dict = {
        f"{prefix}het": "|".join(
            map(lambda x: f"{x:.1f}", ht.head(1).age_hist_het.collect()[0].bin_edges)
        ),
        f"{prefix}hom": "|".join(
            map(lambda x: f"{x:.1f}", ht.head(1).age_hist_hom.collect()[0].bin_edges)
        ),
    }
    for hist in HISTS:

        # Parse hists calculated on both raw and adj-filtered data
        for hist_type in ["adj_qual_hists", "qual_hists"]:
            hist_name = hist
            if "adj" in hist_type:
                hist_name = f"{hist}_adj"
            edges_dict[hist_name] = (
                "|".join(
                    map(
                        lambda x: f"{x:.2f}",
                        ht.head(1)[hist_type][hist].collect()[0].bin_edges,
                    )
                )
                if "ab" in hist
                else "|".join(
                    map(
                        lambda x: str(int(x)),
                        ht.head(1)[hist_type][hist].collect()[0].bin_edges,
                    )
                )
            )
    return edges_dict


def make_hist_dict(bin_edges: Dict[str, Dict[str, str]], adj: bool) -> Dict[str, str]:
    """
    Generate dictionary of Number and Description attributes to be used in the VCF header, specifically for histogram annotations.

    :param bin_edges: Dictionary keyed by histogram annotation name, with corresponding string-reformatted bin edges for values.
    :param adj: Whether to create a header dict for raw or adj quality histograms.
    :return: Dictionary keyed by VCF INFO annotations, where values are Dictionaries of Number and Description attributes.
    """
    header_hist_dict = {}
    for hist in HISTS:

        # Get hists for both raw and adj data
        if adj:
            hist = f"{hist}_adj"

        edges = bin_edges[hist]
        hist_fields = hist.split("_")
        hist_text = hist_fields[0].upper()

        if hist_fields[2] == "alt":
            hist_text = hist_text + " in heterozygous individuals"
        if adj:
            hist_text = hist_text + " calculated on high quality genotypes"

        hist_dict = {
            f"{hist}_bin_freq": {
                "Number": "A",
                "Description": f"Histogram for {hist_text}; bin edges are: {edges}",
            },
            f"{hist}_n_smaller": {
                "Number": "A",
                "Description": f"Count of {hist_fields[0].upper()} values falling below lowest histogram bin edge {hist_text}",
            },
            f"{hist}_n_larger": {
                "Number": "A",
                "Description": f"Count of {hist_fields[0].upper()} values falling above highest histogram bin edge {hist_text}",
            },
        }

        header_hist_dict.update(hist_dict)

    return header_hist_dict


def sample_sum_check(
    ht: hl.Table,
    prefix: str,
    label_groups: Dict[str, List[str]],
    verbose: bool,
    subpop: bool = None,
) -> None:
    """
    Compute afresh the sum of annotations for a specified group of annotations, and compare to the annotated version;
    display results from checking the sum of the specified annotations in the terminal.

    :param Table ht: Table containing annotations to be summed.
    :param str prefix: String indicating sample subset.
    :param dict label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"]).
    :param bool verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param str subpop: Subpop abbreviation, supplied only if subpopulations are included in the annotation groups being checked.
    :rtype: None
    """
    combo_AC = [ht.info[f"{prefix}AC_{x}"] for x in make_label_combos(label_groups)]
    combo_AN = [ht.info[f"{prefix}AN_{x}"] for x in make_label_combos(label_groups)]
    combo_nhomalt = [
        ht.info[f"{prefix}nhomalt_{x}"] for x in make_label_combos(label_groups)
    ]

    group = label_groups.pop("group")[0]
    alt_groups = "_".join(
        sorted(label_groups.keys(), key=lambda x: SORT_ORDER.index(x))
    )

    annot_dict = {
        f"sum_AC_{group}_{alt_groups}": hl.sum(combo_AC),
        f"sum_AN_{group}_{alt_groups}": hl.sum(combo_AN),
        f"sum_nhomalt_{group}_{alt_groups}": hl.sum(combo_nhomalt),
    }

    ht = ht.annotate(**annot_dict)

    for subfield in ["AC", "AN", "nhomalt"]:
        if not subpop:
            generic_field_check(
                ht,
                (
                    ht.info[f"{prefix}{subfield}_{group}"]
                    != ht[f"sum_{subfield}_{group}_{alt_groups}"]
                ),
                f"{prefix}{subfield}_{group} = sum({subfield}_{group}_{alt_groups})",
                [
                    f"info.{prefix}{subfield}_{group}",
                    f"sum_{subfield}_{group}_{alt_groups}",
                ],
                verbose,
            )
        else:
            generic_field_check(
                ht,
                (
                    ht.info[f"{prefix}{subfield}_{group}_{subpop}"]
                    != ht[f"sum_{subfield}_{group}_{alt_groups}"]
                ),
                f"{prefix}{subfield}_{group}_{subpop} = sum({subfield}_{group}_{alt_groups})",
                [
                    f"info.{prefix}{subfield}_{group}_{subpop}",
                    f"sum_{subfield}_{group}_{alt_groups}",
                ],
                verbose,
            )


def set_female_y_metrics_to_na(
    t: Union[hl.Table, hl.MatrixTable], faf: bool = True,
) -> Dict[str, hl.expr.Int32Expression]:
    """
    Set AC, AN, and nhomalt chrY variant annotations for females to NA (instead of 0).

    :param t: Table/MatrixTable containing female variant annotations.
    :return: Dictionary with reset annotations
    """
    metrics = list(t.row.info)
    female_metrics = [x for x in metrics if "_female" in x]
    female_metrics = [
        x for x in female_metrics if ("nhomalt" in x) or ("AC" in x) or ("AN" in x)
    ]

    if faf:
        female_metrics.extend([x for x in female_metrics if "faf" in x])

    female_metrics_dict = {}
    for metric in female_metrics:
        female_metrics_dict.update(
            {
                f"{metric}": hl.or_missing(
                    (t.locus.contig.in_y_nonpar() | t.locus.contig.in_y_par()),
                    t.info[f"{metric}"],
                )
            }
        )
    return female_metrics_dict
