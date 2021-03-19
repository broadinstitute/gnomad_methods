import copy
import itertools
import logging
from typing import Dict, List, Union

import hail as hl

from gnomad.sample_qc.ancestry import POP_NAMES

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

SORT_ORDER = [
    "subset",
    "downsampling",
    "popmax",
    "pop",
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
    "InbreedingCoeff",
]
"""
Allele-specific variant annotations.
"""

SITE_FIELDS = [
    "BaseQRankSum",
    "FS",
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
    decoy and segdup resource files do not currently exist for GRCh38/hg38.
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
        "Number": "A",
        "Description": "Inbreeding coefficient, the excess heterozygosity at a variant site, computed as 1 - (the number of heterozygous genotypes)/(the number of heterozygous genotypes expected under Hardy-Weinberg equilibrium)",
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
    "POSITIVE_TRAIN_SITE": {
        "Description": "Variant was used to build the positive training set of high-quality variants for VQSR"
    },
    "NEGATIVE_TRAIN_SITE": {
        "Description": "Variant was used to build the negative training set of low-quality variants for VQSR"
    },
    "BaseQRankSum": {
        "Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference base qualities",
    },
    "VarDP": {
        "Description": "Depth over variant genotypes (does not include depth of reference samples)"
    },
    "VQSLOD": {
        "Description": "Log-odds ratio of being a true variant versus being a false positive under the trained VQSR Gaussian mixture model",
    },
    "culprit": {
        "Description": "Worst-performing annotation in the VQSR Gaussian mixture model",
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
    "allele_type": {"Description": "Allele type (snv, insertion, deletion, or mixed)",},
    "n_alt_alleles": {
        "Number": "1",
        "Description": "Total number of alternate alleles observed at variant locus",
    },
    "was_mixed": {"Description": "Variant type was mixed"},
    "has_star": {
        "Description": "Variant locus coincides with a spanning deletion (represented by a star) observed elsewhere in the callset"
    },
    "AS_pab_max": {
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
    Create a dictionary keyed by the specified label groupings with values describing the corresponding index of each grouping entry
    in the meta_array annotation

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
    preposition: str, combo_dict: Dict[str, str], pop_names: Dict[str, str],
) -> str:
    """
    Programmatically generate text to populate the VCF header description for a given variant annotation with specific groupings and subset.

    For example, if preposition is "for", group_types is ["group", "pop", "sex"], and combo_fields is ["adj", "afr", "female"],
    this function will return the string " for female samples of African-American/African ancestry".

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
                header_text + f" of {pop_names[combo_dict['subpop']]} ancestry"
            )

        else:
            header_text = header_text + f" of {pop_names[combo_dict['pop']]} ancestry"

    if "group" in combo_dict:
        if combo_dict["group"] == "raw":
            header_text = header_text + ", before removing low-confidence genotypes"

    return header_text


def make_info_dict(
    prefix: str = "",
    pop_names: Dict[str, str] = POP_NAMES,
    label_groups: Dict[str, str] = None,
    label_delimiter: str = "_",
    bin_edges: Dict[str, str] = None,
    faf: bool = False,
    popmax: bool = False,
    description_text: str = "",
    age_hist_data: str = None,
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
    :param pop_names: Dict with global population names (keys) and population descriptions (values). Default is POP_NAMES.
    :param label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"]).
    :param label_delimiter: String to use as delimiter when making group label combinations.
    :param bin_edges: Dictionary keyed by annotation type, with values that reflect the bin edges corresponding to the annotation.
    :param faf: If True, use alternate logic to auto-populate dictionary values associated with filter allele frequency annotations.
    :param popmax: If True, use alternate logic to auto-populate dictionary values associated with popmax annotations.
    :param description_text: Optional text to append to the end of descriptions. Needs to start with a space if specified.
    :param str age_hist_data: Pipe-delimited string of age histograms, from `get_age_distributions`.
    :param sort_order: List containing order to sort label group combinations. Default is SORT_ORDER.
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
                "Description": f"Histogram of ages of homozygous alternate individuals{description_text}; bin edges are: {bin_edges['hom']}; total number of individuals of any genotype bin: {age_hist_data}",
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
                "Description": f"Maximum allele frequency across populations{description_text}",
            },
            f"{prefix}nhomalt_popmax": {
                "Number": "A",
                "Description": f"Count of homozygous individuals in the population with the maximum allele frequency{description_text}",
            },
        }
        info_dict.update(popmax_dict)

    else:
        group_types = sorted(label_groups.keys(), key=lambda x: sort_order.index(x))
        combos = make_label_combos(label_groups, label_delimiter=label_delimiter)

        for combo in combos:
            combo_fields = combo.split(label_delimiter)
            group_dict = dict(zip(group_types, combo_fields))

            for_combo = make_combo_header_text("for", group_dict, pop_names)
            in_combo = make_combo_header_text("in", group_dict, pop_names)

            if not faf:
                combo_dict = {
                    f"{prefix}AC_{combo}": {
                        "Number": "A",
                        "Description": f"Alternate allele count{for_combo}{description_text}",
                    },
                    f"{prefix}AN_{combo}": {
                        "Number": "1",
                        "Description": f"Total number of alleles{in_combo}{description_text}",
                    },
                    f"{prefix}AF_{combo}": {
                        "Number": "A",
                        "Description": f"Alternate allele frequency{in_combo}{description_text}",
                    },
                    f"{prefix}nhomalt_{combo}": {
                        "Number": "A",
                        "Description": f"Count of homozygous individuals{in_combo}{description_text}",
                    },
                }
            else:
                combo_dict = {
                    f"{prefix}faf95_{combo}": {
                        "Number": "A",
                        "Description": f"Filtering allele frequency (using Poisson 95% CI){for_combo}{description_text}",
                    },
                    f"{prefix}faf99_{combo}": {
                        "Number": "A",
                        "Description": f"Filtering allele frequency (using Poisson 99% CI){for_combo}{description_text}",
                    },
                }
            info_dict.update(combo_dict)
    return info_dict


def add_as_info_dict(
    info_dict: Dict[str, Dict[str, str]] = INFO_DICT, as_fields: List[str] = AS_FIELDS
) -> Dict[str, Dict[str, str]]:
    """
    Updates info dictionary with allele-specific terms and their descriptions.

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
            logger.warning(f"{field} is not present in input info dictionary!")

    return as_dict


def make_vcf_filter_dict(
    snp_cutoff: float, indel_cutoff: float, inbreeding_cutoff: float
) -> Dict[str, str]:
    """
    Generates dictionary of Number and Description attributes to be used in the VCF header, specifically for FILTER annotations.

    Generates descriptions for:
        - AC0 filter
        - InbreedingCoeff filter
        - RF filter
        - PASS (passed all variant filters)

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
            "Description": f"Failed random forest filtering thresholds of {snp_cutoff} for SNPs and {indel_cutoff} for indels (probabilities of being a true positive variant)"
        },
        "PASS": {"Description": "Passed all variant filters"},
    }
    return filter_dict


def make_hist_bin_edges_expr(
    ht: hl.Table, hists: List[str] = HISTS, prefix: str = ""
) -> Dict[str, str]:
    """
    Create dictionaries containing variant histogram annotations and their associated bin edges, formatted into a string
    separated by pipe delimiters.

    :param ht: Table containing histogram variant annotations.
    :param hists: List of variant histogram annotations. Default is HISTS.
    :param prefix: Prefix text for age histogram bin edges.  Default is empty string.
    :return: Dictionary keyed by histogram annotation name, with corresponding reformatted bin edges for values.
    """
    # Add underscore to prefix if it isn't empty
    if prefix != "":
        prefix += "_"

    edges_dict = {
        f"{prefix}{call_type}": "|".join(
            map(
                lambda x: f"{x:.1f}",
                ht.head(1)[f"age_hist_{call_type}"].collect()[0].bin_edges,
            )
        )
        for call_type in ["het", "hom"]
    }

    for hist in hists:

        # Parse hists calculated on both raw and adj-filtered data
        for hist_type in ["raw_qual_hists", "qual_hists"]:

            hist_name = hist
            if "raw" in hist_type:
                hist_name = f"{hist}_raw"

            edges_dict[hist_name] = "|".join(
                map(
                    lambda x: f"{x:.2f}" if "ab" in hist else str(int(x)),
                    ht.head(1)[hist_type][hist].collect()[0].bin_edges,
                )
            )
    return edges_dict


def make_hist_dict(
    bin_edges: Dict[str, Dict[str, str]], adj: bool, label_delimiter: str = "_"
) -> Dict[str, str]:
    """
    Generate dictionary of Number and Description attributes to be used in the VCF header, specifically for histogram annotations.

    :param bin_edges: Dictionary keyed by histogram annotation name, with corresponding string-reformatted bin edges for values.
    :param adj: Whether to create a header dict for raw or adj quality histograms.
    :param label_delimiter: String used as delimiter when making group label combinations.
    :return: Dictionary keyed by VCF INFO annotations, where values are Dictionaries of Number and Description attributes.
    """
    header_hist_dict = {}
    for hist in HISTS:

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
                    (~t.locus.in_y_nonpar() & ~t.locus.in_y_par()), t.info[f"{metric}"],
                )
            }
        )
    return female_metrics_dict
