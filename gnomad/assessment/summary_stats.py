# noqa: D100
import itertools
import logging
from copy import deepcopy
from typing import Dict, List, Optional, Set, Union

import hail as hl

from gnomad.utils.filtering import filter_low_conf_regions, low_conf_regions_expr
from gnomad.utils.vep import (
    LOF_CSQ_SET,
    add_most_severe_consequence_to_consequence,
    filter_vep_to_canonical_transcripts,
    filter_vep_to_mane_select_transcripts,
    get_most_severe_consequence_for_summary,
    process_consequences,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def freq_bin_expr(
    freq_expr: hl.expr.ArrayExpression, index: int = 0
) -> hl.expr.StringExpression:
    """
    Return frequency string annotations based on input AC or AF.

    .. note::

        - Default index is 0 because function assumes freq_expr was calculated with `annotate_freq`.
        - Frequency index 0 from `annotate_freq` is frequency for all pops calculated on adj genotypes only.

    :param freq_expr: Array of structs containing frequency information.
    :param index: Which index of freq_expr to use for annotation. Default is 0.
    :return: StringExpression containing bin name based on input AC or AF.
    """
    return (
        hl.case(missing_false=True)
        .when(freq_expr[index].AC == 0, "Not found")
        .when(freq_expr[index].AC == 1, "Singleton")
        .when(freq_expr[index].AC == 2, "Doubleton")
        .when(freq_expr[index].AC <= 5, "AC 3 - 5")
        .when(freq_expr[index].AF < 1e-4, "AC 6 - 0.01%")
        .when(freq_expr[index].AF < 1e-3, "0.01% - 0.1%")
        .when(freq_expr[index].AF < 1e-2, "0.1% - 1%")
        .when(freq_expr[index].AF < 1e-1, "1% - 10%")
        .when(freq_expr[index].AF > 0.95, ">95%")
        .default("10% - 95%")
    )


def get_summary_counts_dict(
    locus_expr: hl.expr.LocusExpression,
    allele_expr: hl.expr.ArrayExpression,
    lof_expr: hl.expr.StringExpression,
    no_lof_flags_expr: hl.expr.BooleanExpression,
    most_severe_csq_expr: hl.expr.StringExpression,
    prefix_str: str = "",
) -> Dict[str, hl.expr.Int64Expression]:
    """
    Return dictionary containing containing counts of multiple variant categories.

    Categories are:
        - Number of variants
        - Number of indels
        - Number of SNVs
        - Number of LoF variants
        - Number of LoF variants that pass LOFTEE
        - Number of LoF variants that pass LOFTEE without any flgs
        - Number of LoF variants annotated as 'other splice' (OS) by LOFTEE
        - Number of LoF variants that fail LOFTEE
        - Number of missense variants
        - Number of synonymous variants
        - Number of autosomal variants
        - Number of allosomal variants

    .. warning::
        Assumes `allele_expr` contains only two variants (multi-allelics have been split).

    :param locus_expr: LocusExpression.
    :param allele_expr: ArrayExpression containing alleles.
    :param lof_expr: StringExpression containing LOFTEE annotation.
    :param no_lof_flags_expr: BooleanExpression indicating whether LoF variant has any flags.
    :param most_severe_csq_expr: StringExpression containing most severe consequence annotation.
    :param prefix_str: Desired prefix string for category names. Default is empty str.
    :return: Dict of categories and counts per category.
    """
    logger.warning("This function expects that multi-allelic variants have been split!")
    return {
        f"{prefix_str}num_variants": hl.agg.count(),
        f"{prefix_str}indels": hl.agg.count_where(
            hl.is_indel(allele_expr[0], allele_expr[1])
        ),
        f"{prefix_str}snps": hl.agg.count_where(
            hl.is_snp(allele_expr[0], allele_expr[1])
        ),
        f"{prefix_str}LOF": hl.agg.count_where(hl.is_defined(lof_expr)),
        f"{prefix_str}pass_loftee": hl.agg.count_where(lof_expr == "HC"),
        f"{prefix_str}pass_loftee_no_flag": hl.agg.count_where(
            (lof_expr == "HC") & (no_lof_flags_expr)
        ),
        f"{prefix_str}loftee_os": hl.agg.count_where(lof_expr == "OS"),
        f"{prefix_str}fail_loftee": hl.agg.count_where(lof_expr == "LC"),
        f"{prefix_str}num_missense": hl.agg.count_where(
            most_severe_csq_expr == "missense_variant"
        ),
        f"{prefix_str}num_synonymous": hl.agg.count_where(
            most_severe_csq_expr == "synonymous_variant"
        ),
        f"{prefix_str}num_autosomal_variants": hl.agg.filter(
            locus_expr.in_autosome_or_par(), hl.agg.count()
        ),
        f"{prefix_str}num_allosomal_variants": hl.agg.filter(
            locus_expr.in_x_nonpar() | locus_expr.in_y_nonpar(), hl.agg.count()
        ),
    }


def get_summary_ac_dict(
    ac_expr: hl.expr.Int64Expression,
    lof_expr: hl.expr.StringExpression,
    no_lof_flags_expr: hl.expr.BooleanExpression,
    most_severe_csq_expr: hl.expr.StringExpression,
) -> Dict[str, hl.expr.Int64Expression]:
    """
    Return dictionary containing containing total allele counts for variant categories.

    Categories are:
        - All variants
        - LoF variants
        - LoF variants that pass LOFTEE
        - LoF variants that pass LOFTEE without any flags
        - LoF variants that are annotate as 'other splice' (OS) by LOFTEE
        - LoF variants that fail LOFTEE
        - Missense variants
        - Synonymous variants

    .. warning::
        Assumes `allele_expr` contains only two variants (multi-allelics have been split).

    :param allele_expr: ArrayExpression containing alleles.
    :param lof_expr: StringExpression containing LOFTEE annotation.
    :param no_lof_flags_expr: BooleanExpression indicating whether LoF variant has any flags.
    :return: Dict of variant categories and their total allele counts.
    """
    logger.warning("This function expects that multi-allelic variants have been split!")
    return {
        "total_ac": hl.agg.sum(ac_expr),
        "total_ac_LOF": hl.agg.filter(hl.is_defined(lof_expr), hl.agg.sum(ac_expr)),
        "total_ac_pass_loftee": hl.agg.filter(lof_expr == "HC", hl.agg.sum(ac_expr)),
        "total_ac_pass_loftee_no_flag": hl.agg.filter(
            (lof_expr == "HC") & (no_lof_flags_expr), hl.agg.sum(ac_expr)
        ),
        "total_ac_loftee_os": hl.agg.filter(lof_expr == "OS", hl.agg.sum(ac_expr)),
        "total_ac_fail_loftee": hl.agg.filter(lof_expr == "LC", hl.agg.sum(ac_expr)),
        "total_ac_missense": hl.agg.filter(
            most_severe_csq_expr == "missense_variant", hl.agg.sum(ac_expr)
        ),
        "total_ac_synonymous": hl.agg.filter(
            most_severe_csq_expr == "synonymous_variant", hl.agg.sum(ac_expr)
        ),
    }


def get_summary_counts(
    ht: hl.Table,
    freq_field: str = "freq",
    filter_field: str = "filters",
    filter_decoy: bool = False,
    canonical_only: bool = True,
    mane_select_only: bool = False,
    index: int = 0,
) -> hl.Table:
    """
    Generate a struct with summary counts across variant categories.

    Summary counts:
        - Number of variants
        - Number of indels
        - Number of SNVs
        - Number of LoF variants
        - Number of LoF variants that pass LOFTEE (including with LoF flags)
        - Number of LoF variants that pass LOFTEE without LoF flags
        - Number of OS (other splice) variants annotated by LOFTEE
        - Number of LoF variants that fail LOFTEE filters

    Also annotates Table's globals with total variant counts.

    Before calculating summary counts, function:
        - Filters out low confidence regions
        - Uses the most severe consequence
        - Filters to canonical transcripts (if `canonical_only` is True) or MANE Select
          transcripts (if `mane_select_only` is True)

    Assumes that:
        - Input HT is annotated with VEP.
        - Multiallelic variants have been split and/or input HT contains bi-allelic variants only.
        - freq_expr was calculated with `annotate_freq`.
        - (Frequency index 0 from `annotate_freq` is frequency for all pops calculated on adj genotypes only.)

    :param ht: Input Table.
    :param freq_field: Name of field in HT containing frequency annotation (array of structs). Default is "freq".
    :param filter_field: Name of field in HT containing variant filter information. Default is "filters".
    :param canonical_only: Whether to filter to canonical transcripts. Default is True.
    :param mane_select_only: Whether to filter to MANE Select transcripts. Default is False.
    :param filter_decoy: Whether to filter decoy regions. Default is False.
    :param index: Which index of freq_expr to use for annotation. Default is 0.
    :return: Table grouped by frequency bin and aggregated across summary count categories.
    """
    if canonical_only and mane_select_only:
        raise ValueError(
            "Only one of `canonical_only` and `mane_select_only` can be True."
        )

    logger.info("Checking if multi-allelic variants have been split...")
    max_alleles = ht.aggregate(hl.agg.max(hl.len(ht.alleles)))
    if max_alleles > 2:
        logger.info("Splitting multi-allelics and VEP transcript consequences...")
        ht = hl.split_multi_hts(ht)

    logger.info("Filtering to PASS variants in high confidence regions...")
    ht = ht.filter((hl.len(ht[filter_field]) == 0))
    ht = filter_low_conf_regions(ht, filter_decoy=filter_decoy)

    if canonical_only:
        logger.info("Filtering to canonical transcripts...")
        ht = filter_vep_to_canonical_transcripts(ht)
    elif mane_select_only:
        logger.info("Filtering to mane select transcripts...")
        ht = filter_vep_to_mane_select_transcripts(ht)

    logger.info("Getting VEP summary annotations...")
    ht = get_most_severe_consequence_for_summary(ht)

    logger.info("Annotating with frequency bin information...")
    ht = ht.annotate(freq_bin=freq_bin_expr(ht[freq_field], index))

    logger.info(
        "Annotating HT globals with total counts/total allele counts per variant"
        " category..."
    )
    summary_counts = ht.aggregate(
        hl.struct(
            **get_summary_counts_dict(
                ht.locus,
                ht.alleles,
                ht.lof,
                ht.no_lof_flags,
                ht.most_severe_csq,
                prefix_str="total_",
            )
        )
    )
    summary_ac_counts = ht.aggregate(
        hl.struct(
            **get_summary_ac_dict(
                ht[freq_field][index].AC,
                ht.lof,
                ht.no_lof_flags,
                ht.most_severe_csq,
            )
        )
    )
    ht = ht.annotate_globals(
        summary_counts=summary_counts.annotate(**summary_ac_counts)
    )
    return ht.group_by("freq_bin").aggregate(
        **get_summary_counts_dict(
            ht.locus,
            ht.alleles,
            ht.lof,
            ht.no_lof_flags,
            ht.most_severe_csq,
        )
    )


def get_an_criteria(
    mt: hl.MatrixTable,
    samples_by_sex: Optional[Dict[str, int]] = None,
    meta_root: str = "meta",
    sex_field: str = "sex_imputation.sex_karyotype",
    xy_str: str = "XY",
    xx_str: str = "XX",
    freq_field: str = "freq",
    freq_index: int = 0,
    an_proportion_cutoff: float = 0.8,
) -> hl.expr.BooleanExpression:
    """
    Generate criteria to filter samples based on allele number (AN).

    Uses allele number as proxy for call rate.

    :param mt: Input MatrixTable.
    :param samples_by_sex: Optional Dictionary containing number of samples (value) for each sample sex (key).
    :param meta_root: Name of field in MatrixTable containing sample metadata information. Default is 'meta'.
    :param sex_field: Name of field in MatrixTable containing sample sex assignment. Defualt is 'sex_imputation.sex_karyotype'.
    :param xy_str: String marking whether a sample has XY sex. Default is 'XY'.
    :param xx_str: String marking whether a sample has XX sex. Default is 'XX'.
    :param freq_field: Name of field in MT that contains frequency information. Default is 'freq'.
    :param freq_index: Which index of frequency struct to use. Default is 0.
    :param an_proportion_cutoff: Desired allele number proportion cutoff. Default is 0.8.
    """
    if samples_by_sex is None:
        samples_by_sex = mt.aggregate_cols(hl.agg.counter(mt[meta_root][sex_field]))
    return (
        hl.case()
        .when(
            mt.locus.in_autosome_or_par(),
            mt[freq_field][freq_index].AN
            >= an_proportion_cutoff * 2 * sum(samples_by_sex.values()),
        )
        .when(
            mt.locus.in_x_nonpar(),
            mt[freq_field][freq_index].AN
            >= an_proportion_cutoff
            * (samples_by_sex[xy_str] + samples_by_sex[xx_str] * 2),
        )
        .when(
            mt.locus.in_y_nonpar(),
            mt[freq_field][freq_index].AN
            >= an_proportion_cutoff * samples_by_sex[xy_str],
        )
        .or_missing()
    )


def get_tx_expression_expr(
    key_expr: hl.expr.StructExpression,
    tx_ht: hl.Table,
    csq_expr: hl.expr.StructExpression,
    gene_field: str = "ensg",
    csq_field: str = "csq",
    tx_struct: str = "tx_annotation",
) -> hl.expr.Float64Expression:
    """
    Pull appropriate transcript expression annotation struct given a specific locus and alleles (provided in `key_expr`).

    Assumes that `key_expr` contains a locus and alleles.
    Assumes that multi-allelic variants have been split in both `tx_ht` and `key_expr`.

    :param row_key_expr: StructExpression containing locus and alleles to search in `tx_ht`.
    :param tx_ht: Input Table containing transcript expression information.
    :param csq_expr: Input StructExpression that contains VEP consequence information.
    :param gene_field: Field in `csq_expr` that contains gene ID.
    :param csq_field: Field in `csq_expr` that contains `most_severe_consequence` annotation.
    :param tx_struct: StructExpression that contains transcript expression information.
    :return: StructExpression that contains transcript expression information for given gene ID in `csq_expr`.
    """
    return hl.find(
        lambda csq: (csq[gene_field] == csq_expr.gene_id)
        & (csq[csq_field] == csq_expr.most_severe_consequence),
        tx_ht[key_expr][tx_struct],
    )


def get_summary_stats_variant_filter_expr(
    t: Union[hl.Table, hl.MatrixTable],
    filter_lcr: bool = False,
    filter_expr: hl.expr.SetExpression = None,
    freq_expr: hl.expr.SetExpression = None,
    max_af: Optional[Union[float, List[float]]] = None,
    min_an_proportion: Optional[Union[float, List[float]]] = None,
) -> Dict[str, Union[hl.expr.BooleanExpression, hl.expr.StructExpression]]:
    """
    Generate variant filtering expression for summary stats.

    The possible filtering groups are:

        - 'no_lcr' if `filter_lcr` is True.
        - 'variant_qc_pass' if `filter_expr` is provided.
        - 'max_af' as a struct with a field for each `af` in `max_af` if `max_af` is
          provided.
        - 'min_an' as a struct with a field for each `an_proportion` in
          `min_an_proportion` if `min_an_proportion` is provided.

    :param t: Input Table/MatrixTable.
    :param filter_lcr: Whether to filter out low confidence regions. Default is False.
    :param filter_expr: SetExpression containing variant filters. Default is None.
    :param freq_expr: SetExpression containing frequency information. Default is None.
    :param max_af: Maximum allele frequency cutoff(s). Can be a single float or a list
        of floats. Default is None.
    :param min_an_proportion: Minimum allele number proportion (used as a proxy for
        call rate). Default is None.
    :return: Dict of BooleanExpressions or StructExpressions for filtering variants.
    """
    if max_af is not None and freq_expr is None:
        raise ValueError("Frequency expression must be provided when filtering by AF!")

    log_list = []
    ss_filter_expr = {}
    if filter_lcr:
        log_list.append("variants not in low confidence regions")
        ss_filter_expr["no_lcr"] = low_conf_regions_expr(t.locus, filter_decoy=False)
    if filter_expr is not None:
        log_list.append("variants that pass all variant QC filters")
        ss_filter_expr["variant_qc_pass"] = hl.len(filter_expr) == 0
    if max_af is not None:
        if isinstance(max_af, float):
            max_af = [max_af]
        log_list.extend([f"variants with (AF < {af:.2e})" for af in max_af])
        ss_filter_expr["max_af"] = hl.struct(
            **{f"{af}": freq_expr < af for af in max_af}
        )
    if min_an_proportion is not None:
        log_list.append(
            "variants that meet a minimum call rate of %.2f (using AN as a call rate "
            "proxy)" % min_an_proportion,
        )
        ss_filter_expr[f"min_an_{min_an_proportion}"] = get_an_criteria(
            t, an_proportion_cutoff=min_an_proportion
        )

    logger.info("Adding filtering for:\n\t%s...", "\n\t".join(log_list))

    return ss_filter_expr


def get_summary_stats_csq_filter_expr(
    t: Union[hl.Table, hl.MatrixTable, hl.StructExpression],
    lof_csq_set: Optional[Set[str]] = None,
    lof_label_set: Optional[Set[str]] = None,
    lof_no_flags: bool = False,
    lof_any_flags: bool = False,
    additional_csq_sets: Optional[Dict[str, Set[str]]] = None,
    additional_csqs: Optional[Set[str]] = None,
) -> Dict[str, Union[hl.expr.BooleanExpression, hl.expr.StructExpression]]:
    """
    Generate consequence filtering expression for summary stats.

    .. note::

        - Assumes that input Table/MatrixTable/StructExpression contains the required
          annotations for the requested filtering groups.

            - 'lof' annotation for `lof_csq_set`.
            - 'no_lof_flags' annotation for `lof_no_flags` and `lof_any_flags`.

    The possible filtering groups are:

        - 'lof' if `lof_csq_set` is provided.
        - 'loftee_no_flags' if `lof_no_flags` is True.
        - 'loftee_with_flags' if `lof_any_flags` is True.
        - 'loftee_label' as a struct with a field for each `lof_label` in
          `lof_label_set` if provided.
        - 'csq' as a struct with a field for each consequence in `additional_csqs` and
          `lof_csq_set` if provided.
        - 'csq_set' as a struct with a field for each consequence set in
          `additional_csq_sets` if provided. This will also have an `lof` field if
          `lof_csq_set` is provided.

    :param t: Input Table/MatrixTable/StructExpression.
    :param lof_csq_set: Set of LoF consequence strings. Default is None.
    :param lof_label_set: Set of LoF consequence labels. Default is None.
    :param lof_no_flags: Whether to filter to variants with no flags. Default is
        False.
    :param lof_any_flags: Whether to filter to variants with any flags. Default
        is False.
    :param additional_csq_sets: Dictionary containing additional consequence sets.
        Default is None.
    :param additional_csqs: Set of additional consequences to keep. Default is None.
    :return: Dict of BooleanExpressions or StructExpressions for filtering by
        consequence.
    """
    if (lof_no_flags or lof_any_flags) and "no_lof_flags" not in t.row:
        raise ValueError(
            "The required 'no_lof_flags' annotation is not found in input "
            "Table/MatrixTable/StructExpression."
        )
    if lof_label_set and "lof" not in t.row:
        raise ValueError(
            "The required 'lof' annotation is not found in input "
            "Table/MatrixTable/StructExpression."
        )

    def _create_filter_by_csq(
        t: Union[hl.Table, hl.MatrixTable],
        csq_set: Union[Set, List, hl.expr.CollectionExpression],
    ) -> hl.expr.BooleanExpression:
        """
        Create filtering expression for a set of consequences.

        :param t: Input Table/MatrixTable.
        :param csq_set: Set of consequences to filter.
        :return: BooleanExpression for filtering by consequence.
        """
        if not isinstance(csq_set, hl.expr.CollectionExpression):
            csq_set = hl.set(csq_set)

        return csq_set.contains(t.most_severe_csq)

    # Set up filters for specific consequences or sets of consequences.
    csq_filters = {
        "csq": {
            **({f"{c}": {c} for c in lof_csq_set or []}),
            **({f"{c}": {c} for c in additional_csqs or []}),
        },
        "csq_set": {
            **({"lof": lof_csq_set or {}}),
            **({f"{l}": c for l, c in (additional_csq_sets or {}).items()}),
        },
    }

    # Create filtering expressions for each consequence/ consequence set.
    ss_filter_expr = {
        group_name: hl.struct(
            **{
                filter_name: _create_filter_by_csq(t, csq_set)
                for filter_name, csq_set in group.items()
            }
        )
        for group_name, group in csq_filters.items()
    }

    # Add filtering expressions for LOFTEE labels.
    ss_filter_expr["loftee_label"] = hl.struct(
        **{
            lof_label: hl.or_else(t.lof == lof_label, False)
            for lof_label in lof_label_set or []
        }
    )

    # Add filtering expressions variants with no flags or any flags.
    if lof_no_flags:
        ss_filter_expr["loftee_no_flags"] = t.no_lof_flags
    if lof_any_flags:
        ss_filter_expr["loftee_with_flags"] = ~t.no_lof_flags

    return ss_filter_expr


def generate_filter_combinations(
    combos: List[Union[List[str], Dict[str, List[str]]]],
    combo_options: Optional[Dict[str, List[str]]] = None,
) -> List[Dict[str, str]]:
    """
    Generate list of all possible filter combinations from a list of filter options.

    Example input:

    .. code-block:: python

        [
            {'pass_filters': [False, True]},
            {'pass_filters': [False, True], 'capture': ['ukb', 'broad']}
        ]

    Example output:

    .. code-block:: python

        [
            {'pass_filters': False},
            {'pass_filters': True},
            {'pass_filters': False, 'capture': 'ukb'},
            {'pass_filters': False, 'capture': 'broad'},
            {'pass_filters': True, 'capture': 'ukb'},
            {'pass_filters': True, 'capture': 'broad'},
        ]

    :param combos: List of filter groups and their options.
    :param combo_options: Dictionary of filter groups and their options that can be
        supplied if `combos` is a list of lists.
    :return: List of all possible filter combinations for each filter group.
    """
    if isinstance(combos[0], list):
        if combo_options is None:
            raise ValueError(
                "If `combos` is a list of lists, `combo_options` must be provided."
            )
        combos = [{k: combo_options[k] for k in combo} for combo in combos]

    # Create combinations of filter options.

    def _expand_combinations(filter_dict: Dict[str, List[str]]) -> List[Dict[str, str]]:
        """
        Expand filter combinations.

        :param filter_dict: Dictionary of filter options.
        :return: List of dicts of expanded filter options.
        """
        keys, values = zip(*filter_dict.items())
        return [dict(zip(keys, combo)) for combo in itertools.product(*values)]

    # Flatten list of combinations.
    expanded_meta = sum([_expand_combinations(sublist) for sublist in combos], [])

    return expanded_meta


def get_summary_stats_filter_group_meta(
    all_sum_stat_filters: Dict[str, List[str]],
    common_filter_combos: List[List[str]] = None,
    common_filter_override: Dict[str, List[str]] = None,
    lof_filter_combos: Optional[List[List[str]]] = None,
    lof_filter_override: Dict[str, List[str]] = None,
    filter_key_rename: Dict[str, str] = None,
) -> List[Dict[str, str]]:
    """
    Generate list of filter group combination metadata for summary stats.

    This function combines various filter settings for summary statistics and generates
    all possible filter combinations. It ensures that the generated combinations include
    both common filters and specific loss-of-function (LOF) filters.

    .. note::

        - The "variant_qc" filter group is removed if the value is "none", which can
          lead to a filter group of {} (no filters).
        - The `filter_key_rename` parameter can be used to rename keys in the
          `all_sum_stat_filters`, `common_filter_override`, or `lof_filter_override`
          after creating all combinations.

    Example:
        Given the following input:

        .. code-block:: python

            all_sum_stat_filters = {
                "variant_qc": ["none", "pass"],
                "capture": ["1", "2"],
                "max_af": [0.01],
                "lof_csq": ["stop_gained"],
                "lof_csq_set": ["lof"],
            }
            common_filter_combos = [["variant_qc"], ["variant_qc", "capture"]]
            common_filter_override = {"variant_qc": ["pass"], "capture": ["1"]}
            lof_filter_combos = [
                ["lof_csq_set", "loftee_HC"],
                ["lof_csq_set", "loftee_HC", "loftee_flags"],
                ["lof_csq", "loftee_HC", "loftee_flags"],
            ]
            lof_filter_override = {"loftee_HC": ["HC"], "loftee_flags": ["with_flags"]}
            filter_key_rename = {
                "lof_csq": "csq",
                "loftee_HC": "loftee_labels",
                "lof_csq_set": "csq_set",
            }

        The function will generate the following filter combinations:

        .. code-block:: python

            [
               # Combinations of all common filter keys and their possible values.
                {},
                {'capture': '1'},
                {'capture': '2'},
                {'variant_qc': 'pass'},
                {'variant_qc': 'pass', 'capture': '1'},
                {'variant_qc': 'pass', 'capture': '2'},

                # Combinations of all requested common filter combinations with all
                # possible other filter keys and values.
                {'variant_qc': 'pass', 'max_af': '0.01'},
                {'variant_qc': 'pass', 'csq': 'stop_gained'},
                {'variant_qc': 'pass', 'csq_set': 'lof'},
                {'variant_qc': 'pass', 'capture': '1', 'max_af': '0.01'},
                {'variant_qc': 'pass', 'capture': '1', 'csq': 'stop_gained'},
                {'variant_qc': 'pass', 'capture': '1', 'csq_set': 'lof'},

                # Combinations of all requested common filter combinations with all
                # requested LOF filter combination keys and their requested values.
                {'variant_qc': 'pass', 'csq_set': 'lof', 'loftee_labels': 'HC'},
                {
                    'variant_qc': 'pass', 'csq_set': 'lof', 'loftee_labels': 'HC',
                    'loftee_flags': 'with_flags'
                },
                {
                    'variant_qc': 'pass', 'csq': 'stop_gained', 'loftee_labels': 'HC',
                    'loftee_flags': 'with_flags'
                },
                {
                    'variant_qc': 'pass', 'capture': '1', 'csq_set': 'lof',
                    'loftee_labels': 'HC'
                },
                {
                    'variant_qc': 'pass', 'capture': '1', 'csq_set': 'lof',
                    'loftee_labels': 'HC', 'loftee_flags': 'with_flags'
                },
                {
                    'variant_qc': 'pass', 'capture': '1', 'csq': 'stop_gained',
                    'loftee_labels': 'HC', 'loftee_flags': 'with_flags'
                }
            ]

    :param all_sum_stat_filters: Dictionary of all possible filter types.
    :param common_filter_combos: Optional list of lists of common filter keys to use
        for creating common filter combinations.
    :param common_filter_override: Optional dictionary of filter groups and their
        options to override the values in `all_sum_stat_filters` for use with values in
        `common_filter_combos`. This is only used if `common_filter_combos` is not None.
    :param lof_filter_combos: Optional List of loss-of-function keys in all_sum_stat_filters
        to use for creating filter combinations.
    :param lof_filter_override: Optional Dictionary of filter groups and their options
        to override the values in `all_sum_stat_filters` for use with values in
        `lof_combos`. This is only used if `lof_filter_combos` is not None.
    :param filter_key_rename: Optional dictionary to rename keys in
        `all_sum_stat_filters`, `common_filter_override`, or `lof_filter_override` to
        final metadata keys.
    :return: Dictionary of filter field to metadata.
    """
    if common_filter_override is not None and common_filter_combos is None:
        raise ValueError(
            "If `common_combo_override` is provided, `common_filter_combos` must be "
            "provided."
        )
    if lof_filter_override is not None and lof_filter_combos is None:
        raise ValueError(
            "If `lof_filter_override` is provided, `lof_filter_combos` must be "
            "provided."
        )

    # Initialize dictionaries and lists to handle cases where the optional parameters
    # are not provided.
    common_filter_override = common_filter_override or {}
    lof_filter_override = lof_filter_override or {}
    filter_key_rename = filter_key_rename or {}

    # Generate all possible filter combinations for common filter combinations.
    all_sum_stat_filters = deepcopy(all_sum_stat_filters)
    filter_combinations = []
    if common_filter_combos is not None:
        filter_combinations.extend(
            generate_filter_combinations(
                [
                    {f: all_sum_stat_filters[f] for f in combo}
                    for combo in common_filter_combos
                ]
            )
        )
    common_filter_combos = common_filter_combos or [[]]

    # Update the all_sum_stat_filters dictionary with common_filter_override values and
    # remove all filters in common_filter_combos from all_sum_stat_filters and
    # into a common_filters dictionary.
    all_sum_stat_filters.update(common_filter_override)
    common_filters = {
        k: all_sum_stat_filters.pop(k) for k in set(sum(common_filter_combos, []))
    }

    # Add combinations of common filters with all other filters.
    filter_combinations.extend(
        generate_filter_combinations(
            [c + [f] for c in common_filter_combos for f in all_sum_stat_filters],
            {**all_sum_stat_filters, **common_filters},
        )
    )

    if lof_filter_combos is not None:
        # Add combinations of common filters with LOF specific filters.
        all_sum_stat_filters.update(lof_filter_override)
        filter_combinations.extend(
            generate_filter_combinations(
                [c + f for c in common_filter_combos for f in lof_filter_combos],
                {**all_sum_stat_filters, **common_filters},
            )
        )

    filter_combinations = [
        {
            filter_key_rename.get(str(k), str(k)): str(v)
            for k, v in filter_group.items()
            if not (k == "variant_qc" and v == "none")
        }
        for filter_group in filter_combinations
    ]

    return filter_combinations


def default_generate_gene_lof_matrix(
    mt: hl.MatrixTable,
    tx_ht: Optional[hl.Table],
    high_expression_cutoff: float = 0.9,
    low_expression_cutoff: float = 0.1,
    filter_field: str = "filters",
    freq_field: str = "freq",
    freq_index: int = 0,
    additional_csq_set: Set[str] = {"missense_variant", "synonymous_variant"},
    all_transcripts: bool = False,
    filter_an: bool = False,
    filter_to_rare: bool = False,
    pre_loftee: bool = False,
    lof_csq_set: Set[str] = LOF_CSQ_SET,
    remove_ultra_common: bool = False,
) -> hl.MatrixTable:
    """
    Generate loss-of-function gene matrix.

    Used to generate summary metrics on LoF variants.

    :param mt: Input MatrixTable.
    :param tx_ht: Optional Table containing expression levels per transcript.
    :param high_expression_cutoff: Minimum mean proportion expressed cutoff for a transcript to be considered highly expressed. Default is 0.9.
    :param low_expression_cutoff: Upper mean proportion expressed cutoff for a transcript to lowly expressed. Default is 0.1.
    :param filter_field: Name of field in MT that contains variant filters. Default is 'filters'.
    :param freq_field: Name of field in MT that contains frequency information. Default is 'freq'.
    :param freq_index: Which index of frequency struct to use. Default is 0.
    :param additional_csq_set: Set of additional consequences to keep. Default is {'missense_variant', 'synonymous_variant'}.
    :param all_transcripts: Whether to use all transcripts instead of just the transcript with most severe consequence. Default is False.
    :param filter_an: Whether to filter using allele number as proxy for call rate. Default is False.
    :param filter_to_rare: Whether to filter to rare (AF < 5%) variants. Default is False.
    :param pre_loftee: Whether LoF consequences have been annotated with LOFTEE. Default is False.
    :param lof_csq_set: Set of LoF consequence strings. Default is {"splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant"}.
    :param remove_ultra_common: Whether to remove ultra common (AF > 95%) variants. Default is False.
    """
    logger.info("Filtering to PASS variants...")
    filt_criteria = hl.len(mt[filter_field]) == 0
    if filter_an:
        logger.info(
            "Using AN (as a call rate proxy) to filter to variants that meet a minimum"
            " call rate..."
        )
        mt = mt.filter_rows(get_an_criteria(mt))
    if remove_ultra_common:
        logger.info("Removing ultra common (AF > 95%) variants...")
        filt_criteria &= mt[freq_field][freq_index].AF < 0.95
    if filter_to_rare:
        logger.info("Filtering to rare (AF < 5%) variants...")
        filt_criteria &= mt[freq_field][freq_index].AF < 0.05
    mt = mt.filter_rows(filt_criteria)

    if all_transcripts:
        logger.info("Exploding transcript_consequences field...")
        explode_field = "transcript_consequences"
    else:
        logger.info(
            "Adding most severe (worst) consequence and expoding worst_csq_by_gene"
            " field..."
        )
        mt = process_consequences(mt)
        explode_field = "worst_csq_by_gene"

    if additional_csq_set:
        logger.info("Including these consequences: %s", additional_csq_set)
        additional_cats = hl.literal(additional_csq_set)

    if pre_loftee:
        logger.info("Filtering to LoF consequences: %s", lof_csq_set)
        lof_cats = hl.literal(lof_csq_set)
        criteria = lambda x: lof_cats.contains(
            add_most_severe_consequence_to_consequence(x).most_severe_consequence
        )
        if additional_csq_set:
            criteria = lambda x: lof_cats.contains(
                add_most_severe_consequence_to_consequence(x).most_severe_consequence
            ) | additional_cats.contains(
                add_most_severe_consequence_to_consequence(x).most_severe_consequence
            )

    else:
        logger.info("Filtering to LoF variants that pass LOFTEE with no LoF flags...")
        criteria = lambda x: (x.lof == "HC") & hl.is_missing(x.lof_flags)
        if additional_csq_set:
            criteria = lambda x: (x.lof == "HC") & hl.is_missing(
                x.lof_flags
            ) | additional_cats.contains(
                add_most_severe_consequence_to_consequence(x).most_severe_consequence
            )

    csqs = mt.vep[explode_field].filter(criteria)
    mt = mt.select_rows(mt[freq_field], csqs=csqs)
    mt = mt.explode_rows(mt.csqs)
    annotation_expr = {
        "gene_id": mt.csqs.gene_id,
        "gene": mt.csqs.gene_symbol,
        "indel": hl.is_indel(mt.alleles[0], mt.alleles[1]),
        "most_severe_consequence": mt.csqs.most_severe_consequence,
    }

    if tx_ht:
        logger.info("Adding transcript expression annotation...")
        tx_annotation = get_tx_expression_expr(
            mt.row_key,
            tx_ht,
            mt.csqs,
        ).mean_proportion
        annotation_expr["expressed"] = (
            hl.case()
            .when(tx_annotation >= high_expression_cutoff, "high")
            .when(tx_annotation > low_expression_cutoff, "medium")
            .when(hl.is_defined(tx_annotation), "low")
            .default("missing")
        )
    else:
        annotation_expr["transcript_id"] = mt.csqs.transcript_id
        annotation_expr["canonical"] = hl.is_defined(mt.csqs.canonical)

    mt = mt.annotate_rows(**annotation_expr)
    return (
        mt.group_rows_by(*list(annotation_expr.keys()))
        .aggregate_rows(
            n_sites=hl.agg.count(),
            n_sites_array=hl.agg.array_sum(mt.freq.map(lambda x: hl.int(x.AC > 0))),
            classic_caf=hl.agg.sum(mt[freq_field][freq_index].AF),
            max_af=hl.agg.max(mt[freq_field][freq_index].AF),
            classic_caf_array=hl.agg.array_sum(mt[freq_field].map(lambda x: x.AF)),
        )
        .aggregate_entries(
            num_homs=hl.agg.count_where(mt.GT.is_hom_var()),
            num_hets=hl.agg.count_where(mt.GT.is_het()),
            defined_sites=hl.agg.count_where(hl.is_defined(mt.GT)),
        )
        .result()
    )


def get_het_hom_summary_dict(
    csq_set: Set[str],
    most_severe_csq_expr: hl.expr.StringExpression,
    defined_sites_expr: hl.expr.Int64Expression,
    num_homs_expr: hl.expr.Int64Expression,
    num_hets_expr: hl.expr.Int64Expression,
    pop_expr: hl.expr.StringExpression,
) -> Dict[str, hl.expr.Int64Expression]:
    """
    Generate dictionary containing summary counts.

    Summary counts are:
        - Number of sites with defined genotype calls
        - Number of samples with heterozygous calls
        - Number of samples with homozygous calls

    Function has option to generate counts by population.

    :param csq_set: Set containing transcript consequence string(s).
    :param most_severe_csq_expr: StringExpression containing most severe consequence.
    :param defined_sites_expr: Int64Expression containing number of sites with defined genotype calls.
    :param num_homs_expr: Int64Expression containing number of samples with homozygous genotype calls.
    :param num_hets_expr: Int64Expression containing number of samples with heterozygous genotype calls.
    :param pop_expr: StringExpression containing sample population labels.
    :return: Dictionary of summary annotation names and their values.
    """
    csq_filter_expr = hl.literal(csq_set).contains(most_severe_csq_expr)
    return {
        "no_alt_calls": hl.agg.count_where(
            (csq_filter_expr)
            & (defined_sites_expr > 0)
            & (num_homs_expr + num_hets_expr == 0)
        ),
        "obs_het": hl.agg.count_where(
            (csq_filter_expr) & (num_homs_expr == 0) & (num_hets_expr > 0)
        ),
        "obs_hom": hl.agg.count_where((csq_filter_expr) & (num_homs_expr > 0)),
        "defined": hl.agg.count_where((csq_filter_expr) & (defined_sites_expr > 0)),
        "pop_no_alt_calls": hl.agg.group_by(
            pop_expr,
            hl.agg.count_where(
                (csq_filter_expr)
                & (defined_sites_expr > 0)
                & (num_homs_expr + num_hets_expr == 0)
            ),
        ),
        "pop_obs_het": hl.agg.group_by(
            pop_expr,
            hl.agg.count_where(
                (csq_filter_expr) & (num_homs_expr == 0) & (num_hets_expr > 0)
            ),
        ),
        "pop_obs_hom": hl.agg.group_by(
            pop_expr,
            hl.agg.count_where((csq_filter_expr) & (num_homs_expr > 0)),
        ),
        "pop_defined": hl.agg.group_by(
            pop_expr,
            hl.agg.count_where((csq_filter_expr) & (defined_sites_expr > 0)),
        ),
    }


def default_generate_gene_lof_summary(
    mt: hl.MatrixTable,
    collapse_indels: bool = False,
    tx: bool = False,
    lof_csq_set: Set[str] = LOF_CSQ_SET,
    meta_root: str = "meta",
    pop_field: str = "pop",
    filter_loftee: bool = False,
) -> hl.Table:
    """
    Generate summary counts for loss-of-function (LoF), missense, and synonymous variants.

    Also calculates p, proportion of of haplotypes carrying a putative LoF (pLoF) variant,
    and observed/expected (OE) ratio of samples with homozygous pLoF variant calls.

    Summary counts are (all per gene):
        - Number of samples with no pLoF variants.
        - Number of samples with heterozygous pLoF variants.
        - Number of samples with homozygous pLoF variants.
        - Total number of sites with genotype calls.
        - All of the above stats grouped by population.

    Assumes MT was created using `default_generate_gene_lof_matrix`.

    .. note::
        Assumes LoF variants in MT were filtered (LOFTEE pass and no LoF flag only).
        If LoF variants have not been filtered and `filter_loftee` is True,
        expects MT has the row annotation `vep`.

    :param mt: Input MatrixTable.
    :param collapse_indels: Whether to collapse indels. Default is False.
    :param tx: Whether input MT has transcript expression data. Default is False.
    :param lof_csq_set: Set containing LoF transcript consequence strings. Default is LOF_CSQ_SET.
    :param meta_root: String indicating top level name for sample metadata. Default is 'meta'.
    :param pop_field: String indiciating field with sample population assignment information. Default is 'pop'.
    :param filter_loftee: Filters to LOFTEE pass variants (and no LoF flags) only. Default is False.
    :return: Table with het/hom summary counts.
    """
    if collapse_indels:
        grouping = ["gene_id", "gene", "most_severe_consequence"]
        if tx:
            grouping.append("expressed")
        else:
            grouping.extend(["transcript_id", "canonical"])
        mt = (
            mt.group_rows_by(*grouping)
            .aggregate_rows(
                n_sites=hl.agg.sum(mt.n_sites),
                n_sites_array=hl.agg.array_sum(mt.n_sites_array),
                classic_caf=hl.agg.sum(mt.classic_caf),
                max_af=hl.agg.max(mt.max_af),
                classic_caf_array=hl.agg.array_sum(mt.classic_caf_array),
            )
            .aggregate_entries(
                num_homs=hl.agg.sum(mt.num_homs),
                num_hets=hl.agg.sum(mt.num_hets),
                defined_sites=hl.agg.sum(mt.defined_sites),
            )
            .result()
        )

    if filter_loftee:
        lof_ht = get_most_severe_consequence_for_summary(mt.rows())
        mt = mt.filter_rows(
            hl.is_defined(lof_ht[mt.row_key].lof)
            & (lof_ht[mt.row_key].lof == "HC")
            & (lof_ht[mt.row_key].no_lof_flags)
        )

    ht = mt.annotate_rows(
        lof=hl.struct(
            **get_het_hom_summary_dict(
                csq_set=lof_csq_set,
                most_severe_csq_expr=mt.most_severe_consequence,
                defined_sites_expr=mt.defined_sites,
                num_homs_expr=mt.num_homs,
                num_hets_expr=mt.num_hets,
                pop_expr=mt[meta_root][pop_field],
            ),
        ),
        missense=hl.struct(
            **get_het_hom_summary_dict(
                csq_set={"missense_variant"},
                most_severe_csq_expr=mt.most_severe_consequence,
                defined_sites_expr=mt.defined_sites,
                num_homs_expr=mt.num_homs,
                num_hets_expr=mt.num_hets,
                pop_expr=mt[meta_root][pop_field],
            ),
        ),
        synonymous=hl.struct(
            **get_het_hom_summary_dict(
                csq_set={"synonymous_variant"},
                most_severe_csq_expr=mt.most_severe_consequence,
                defined_sites_expr=mt.defined_sites,
                num_homs_expr=mt.num_homs,
                num_hets_expr=mt.num_hets,
                pop_expr=mt[meta_root][pop_field],
            ),
        ),
    ).rows()
    ht = ht.annotate(
        p=(1 - hl.sqrt(hl.float64(ht.lof.no_alt_calls) / ht.lof.defined)),
        pop_p=hl.dict(
            hl.array(ht.lof.pop_defined).map(
                lambda x: (
                    x[0],
                    1 - hl.sqrt(hl.float64(ht.lof.pop_no_alt_calls.get(x[0])) / x[1]),
                )
            )
        ),
    )
    ht = ht.annotate(exp_hom_lof=ht.lof.defined * ht.p * ht.p)
    return ht.annotate(oe=ht.lof.obs_hom / ht.exp_hom_lof)
