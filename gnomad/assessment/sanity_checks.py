# noqa: D100

import logging
from typing import Dict, List, Optional, Tuple, Union

import hail as hl
from hail.expr.functions import missing

from gnomad.resources.grch38.gnomad import HGDP_POPS, TGP_POPS, POPS, SEXES
from gnomad.utils.vcf import HISTS, make_label_combos, SORT_ORDER


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def make_field_check_dicts(
    field_check_expr: dict,
    field_check_details: dict,
    check_description: str,
    cond_expr: hl.expr.BooleanExpression,
    display_fields: List[str],
) -> Tuple(dict, dict):
    """
    Create dictionary for aggregating each sanity check's failure count and another for each check's details.

    :param field_check_expr: Dictionary of check description and aggregated expression filtering count
    :param field_check_details: Dictionary of structs containing each check descriptions details
    :param check_description: Check to be added to the dictionary
    :param cond_expr: Logical expression referring to annotations in ht to be checked.
    :param display_fields: List of ht annotations to be displayed in case of failure (for troubleshooting purposes);
        these fields are also displayed if verbose is True.
    :return: Tuple of dictionaries
    """
    field_check_expr[check_description] = hl.agg.filter(cond_expr, hl.agg.count())
    field_check_details[check_description] = hl.struct(
        cond_expr=cond_expr, display_fields=display_fields
    )

    return field_check_expr, field_check_details


def generic_field_check_loop(
    ht: hl.Table,
    field_check_expr: dict,
    field_check_details: dict,
    verbose: bool,
    show_percent_sites: bool = False,
    ht_count: int = None,
):
    """
    Loop through all conditional checks for a given hail Table.

    This loop allows aggregation across the hail Table once, as opposed to aggregating during every conditional check. 
    :param ht: Table containing annotations to be checked.
    :param field_check_expr: Dictionary whose keys are conditions being checked and values are the expressions for filtering to condition.
    :param field_check_details: Dictionary whose keys are the the check descriptions and values are struct expressions used for check and what to display for the check in the terminal.
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param show_percent_sites: Show percentage of sites that fail checks. Default is False.
    :param ht_count: Previously computed sum of sites within hail Table. 
    """
    ht_field_check_counts = ht.aggregate(hl.struct(**field_check_expr))
    for check_description, n_fail in ht_field_check_counts.items():
        generic_field_check(
            ht,
            cond_expr=field_check_details[check_description].cond_expr,
            check_description=check_description,
            n_fail=n_fail,
            display_fields=field_check_details[check_description].display_fields,
            verbose=verbose,
            show_percent_sites=show_percent_sites,
            ht_count=ht_count,
        )


def generic_field_check(
    ht: hl.Table,
    cond_expr: hl.expr.BooleanExpression,
    check_description: str,
    display_fields: List[str],
    verbose: bool = False,
    show_percent_sites: bool = False,
    n_fail: int = None,
    ht_count: int = None,
) -> None:
    """
    Check a generic logical condition involving annotations in a Hail Table and print the results to terminal.
    Displays the number of rows (and percent of rows, if `show_percent_sites` is True) in the Table that match the `cond_expr` and fail to be the desired condition (`check_description`).
    If the number of rows that match the `cond_expr` is 0, then the Table passes that check; otherwise, it fails.
    .. note::
        `cond_expr` and `check_description` are opposites and should never be the same.
        E.g., If `cond_expr` filters for instances where the raw AC is less than adj AC,
        then it is checking sites that fail to be the desired condition (`check_description`)
        of having a raw AC greater than or equal to the adj AC.
    :param ht: Table containing annotations to be checked.
    :param cond_expr: Logical expression referring to annotations in ht to be checked.
    :param check_description: String describing the condition being checked; is displayed in terminal summary message.
    :param display_fields: List of names of ht annotations to be displayed in case of failure (for troubleshooting purposes);
        these fields are also displayed if verbose is True.
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param show_percent_sites: Show percentage of sites that fail checks. Default is False.
    :param n_fail: Previously computed number of sites that fail the conditional checks.
    :param ht_count: Previously computed sum of sites within hail Table. 
    :return: None
    """
    if show_percent_sites and (ht_count is None):
        ht_count = ht.count()

    if n_fail is None and cond_expr:
        n_fail = ht.filter(cond_expr).count()

    if n_fail > 0:
        logger.info("Found %d sites that fail %s check:", n_fail, check_description)
        if show_percent_sites:
            logger.info("Percentage of sites that fail: %d", n_fail / ht_count)
        ht.select(**display_fields).show()
    else:
        logger.info("PASSED %s check", check_description)
        if verbose:
            ht.select(**display_fields).show()


def make_filters_sanity_check_expr(
    ht: hl.Table,
    extra_filter_checks: Optional[Dict[str, hl.expr.Expression]] = None,
    rf: bool = True,
) -> Dict[str, hl.expr.Expression]:
    """
    Make Hail expressions to measure % variants filtered under varying conditions of interest.

    Checks for:
        - Total number of variants
        - Fraction of variants removed due to:
            - Any filter
            - Inbreeding coefficient filter in combination with any other filter
            - AC0 filter in combination with any other filter
            - VQSR or random forest filtering in combination with any other filter
            - Only inbreeding coefficient filter
            - Only AC0 filter
            - Only VQSR or random forest filtering

    :param ht: Table containing 'filter' annotation to be examined.
    :param extra_filter_checks: Optional dictionary containing filter condition name (key) extra filter expressions (value) to be examined.
    :param rf: True if the random forest was used for variant filtration, False if VQSR was used.
    :return: Dictionary containing Hail aggregation expressions to examine filter flags.
    """
    variant_filter_method = "RF" if rf else "VQSR"
    filters_dict = {
        "n": hl.agg.count(),
        "frac_any_filter": hl.agg.fraction(hl.len(ht.filters) != 0),
        "frac_inbreed_coeff": hl.agg.fraction(ht.filters.contains("InbreedingCoeff")),
        "frac_ac0": hl.agg.fraction(ht.filters.contains("AC0")),
        f"frac_{variant_filter_method.lower()}": hl.agg.fraction(
            ht.filters.contains(f"{variant_filter_method}")
        ),
        "frac_inbreed_coeff_only": hl.agg.fraction(
            ht.filters.contains("InbreedingCoeff") & (ht.filters.length() == 1)
        ),
        "frac_ac0_only": hl.agg.fraction(
            ht.filters.contains("AC0") & (ht.filters.length() == 1)
        ),
        f"frac_{variant_filter_method.lower()}_only": hl.agg.fraction(
            ht.filters.contains(f"{variant_filter_method}") & (ht.filters.length() == 1)
        ),
    }
    if extra_filter_checks:
        filters_dict.update(extra_filter_checks)

    return filters_dict


def filters_sanity_check(t: Union[hl.MatrixTable, hl.Table]) -> None:
    """
    Summarize variants filtered under various conditions in input MatrixTable or Table.

    Summarize counts for:
        - Total number of variants
        - Fraction of variants removed due to:
            - Any filter
            - Inbreeding coefficient filter in combination with any other filter
            - AC0 filter in combination with any other filter
            - VQSR filtering in combination with any other filter
            - Only inbreeding coefficient filter
            - Only AC0 filter
            - Only VQSR filtering

    :param t: Input MatrixTable or Table to be checked.
    :return: None
    """
    t = t.rows() if isinstance(t, hl.MatrixTable) else t

    filters = t.aggregate(hl.agg.counter(t.filters))
    logger.info("hl.agg.counter filters: %s", filters)

    filtered_expr = hl.len(t.filters) > 0
    problematic_region_expr = hl.any(
        lambda x: x,
        [
            t.info.lcr,
            t.info.segdup,
            t.info.nonpar,
        ],  # NOTE: in_problematic_region check will need to be updated if we get hg38 decoy
    )

    t = t.annotate(
        is_filtered=filtered_expr, in_problematic_region=problematic_region_expr
    )

    def _filter_agg_order(
        t: Union[hl.MatrixTable, hl.Table],
        group_exprs: Dict[str, hl.expr.Expression],
        n_rows: int = None,
        n_cols: int = None,
        extra_filter_checks: Optional[Dict[str, hl.expr.Expression]] = None,
    ) -> None:
        """
        Perform sanity checks to measure percentages of variants filtered under different conditions.

        :param t: Input MatrixTable or Table.
        :param group_exprs: Dictionary of expressions to group the Table by.
        :param n_rows: Number of rows to show.
        :param n_cols: Number of columns to show.
        :return: None
        """
        t = t.rows() if isinstance(t, hl.MatrixTable) else t
        # NOTE: make_filters_sanity_check_expr returns a dict with %ages of variants filtered
        t.group_by(**group_exprs).aggregate(
            **make_filters_sanity_check_expr(t, extra_filter_checks)
        ).order_by(hl.desc("n")).show(n_rows, n_cols)

    logger.info(
        "Checking distributions of filtered variants amongst variant filters..."
    )

    new_filters_dict = {
        "frac_vqsr": hl.agg.fraction(t.filters.contains("AS_VQSR")),
        "frac_vqsr_only": hl.agg.fraction(
            t.filters.contains("AS_VQSR") & (t.filters.length() == 1)
        ),
    }
    _filter_agg_order(
        t, {"is_filtered": t.is_filtered}, extra_filter_checks=new_filters_dict
    )

    logger.info("Checking distributions of variant type amongst variant filters...")
    _filter_agg_order(
        t, {"allele_type": t.info.allele_type}, extra_filter_checks=new_filters_dict
    )

    logger.info(
        "Checking distributions of variant type and region type amongst variant filters..."
    )
    _filter_agg_order(
        t,
        {
            "allele_type": t.info.allele_type,
            "in_problematic_region": t.in_problematic_region,
        },
        50,
        140,
        extra_filter_checks=new_filters_dict,
    )

    logger.info(
        "Checking distributions of variant type, region type, and number of alt alleles amongst variant filters..."
    )
    _filter_agg_order(
        t,
        {
            "allele_type": t.info.allele_type,
            "in_problematic_region": t.in_problematic_region,
            "n_alt_alleles": t.info.n_alt_alleles,
        },
        50,
        140,
        extra_filter_checks=new_filters_dict,
    )


def subset_freq_sanity_checks(
    t: Union[hl.MatrixTable, hl.Table],
    subsets: List[str],
    verbose: bool,
    show_percent_sites: bool = True,
    delimiter: str = "-",
    metric_first_label: bool = True,
) -> None:
    """
    Perform sanity checks on frequency data in input Table.

    Check:
        - Number of sites where callset frequency is equal to a subset frequency (raw and adj)
            - eg. t.info.AC-adj != t.info.AC-subset1-adj 
        - Total number of sites where the allele count annotation is defined (raw and adj)
        
    :param t: Input MatrixTable or Table.
    :param subsets: List of sample subsets.
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param show_percent_sites: If true, show the percentage and count of overall sites that fail; if False, only show the number of sites that fail.
    :param delimiter: String to use as delimiter when making group label combinations.
    :param metric_first_label: If True, metric precedes label group, e.g. AC-afr-male. If False, label group precedes metric, afr-male-AC.
    :return: None
    """
    t = t.rows() if isinstance(t, hl.MatrixTable) else t
    field_check_expr = {}
    field_check_details = {}
    for subset in subsets:
        if subset != "":
            subset += delimiter
            for field in ["AC", "AN", "nhomalt"]:
                for group in ["adj", "raw"]:
                    logger.info(
                        "Comparing subset %s frequencies to entire callset", subset
                    )
                    check_field_left = f"{field}{delimiter}{group}"
                    if metric_first_label:
                        check_field_right = f"{field}{delimiter}{subset}{group}"
                    else:
                        check_field_right = f"{subset}{field}{delimiter}{group}"
                    field_check_expr, field_check_details = make_field_check_dicts(
                        field_check_expr=field_check_expr,
                        field_check_details=field_check_details,
                        check_description=f"{check_field_left} != {check_field_right}",
                        cond_expr=t.info[check_field_left] == t.info[check_field_right],
                        display_fields=hl.struct(
                            **{
                                check_field_left: t.info[check_field_left],
                                check_field_right: t.info[check_field_right],
                            }
                        ),
                    )

    generic_field_check_loop(
        t,
        field_check_expr,
        field_check_details,
        verbose,
        show_percent_sites=show_percent_sites,
    )

    freq_counts = t.aggregate(
        hl.struct(
            total_defined_AC=hl.agg.count_where(
                hl.is_defined(t.info[f"AC{delimiter}adj"])
            ),
            total_defined_AC_raw=hl.agg.count_where(
                hl.is_defined(t.info[f"AC{delimiter}raw"])
            ),
        )
    )
    logger.info("Frequency spot check counts: %d", freq_counts)


def sample_sum_check(
    t: Union[hl.MatrixTable, hl.Table],
    subset: str,
    label_groups: Dict[str, List[str]],
    sort_order: List[str] = SORT_ORDER,
    delimiter: str = "-",
    metric_first_label: bool = True,
) -> Tuple(dict, dict):
    """
    Compute the sum of call stats annotations for a specified group of annotations, compare to the annotated version, and display the result in stdout.

    For example, if pop1 consists of pop1, pop2, and pop3, check that t.info.AC-subset1 == sum(t.info.AC-subset1-pop1, t.info.AC-subset1-pop2, t.info.AC-subset1-pop3).

    :param t: Input MatrixTable or Table containing annotations to be summed.
    :param subset: String indicating sample subset.
    :param label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["XY", "XX"] or ["afr", "nfe", "amr"]).
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param sort_order: List containing order to sort label group combinations. Default is SORT_ORDER.
    :param delimiter: String to use as delimiter when making group label combinations.
    :param metric_first_label: If True, metric precedes label group, e.g. AC-afr-male. If False, label group precedes metric, afr-male-AC.
    :return: Tuple of dictionaries
    """

    t = t.rows() if isinstance(t, hl.MatrixTable) else t

    if subset:  # TODO: add & subset != "":
        subset += delimiter

    label_combos = make_label_combos(label_groups, label_delimiter=delimiter)
    # Grab the adj group for checks
    group = label_groups.pop("group")[0]
    alt_groups = delimiter.join(
        sorted(label_groups.keys(), key=lambda x: sort_order.index(x))
    )
    info_fields = t.info.keys()

    annot_dict = {}
    for subfield in ["AC", "AN", "nhomalt"]:
        subfield_values = []
        for x in label_combos:
            if metric_first_label:
                field = f"{subfield}{delimiter}{subset}{x}"
            else:
                field = f"{subset}{subfield}{delimiter}{x}"

            if field in info_fields:
                subfield_values.append(t.info[field])
            else:
                logger.info("%s is not in table's info field", field)

        annot_dict[f"sum_{subfield}"] = hl.sum(subfield_values)

    field_check_expr = {}
    field_check_details = {}
    for subfield in ["AC", "AN", "nhomalt"]:
        if metric_first_label:
            check_field_left = f"{subfield}{delimiter}{subset}{group}"
        else:
            check_field_left = f"{subset}{subfield}{delimiter}{group}"

        check_field_right = f"sum_{subfield}{delimiter}{group}{delimiter}{alt_groups}"
        field_check_expr, field_check_details = make_field_check_dicts(
            field_check_expr=field_check_expr,
            field_check_details=field_check_details,
            check_description=f"{check_field_left} = {check_field_right}",
            cond_expr=t.info[check_field_left] != annot_dict[f"sum_{subfield}"],
            display_fields=hl.struct(
                **{
                    check_field_left: t.info[check_field_left],
                    f"sum_{subfield}": annot_dict[f"sum_{subfield}"],
                }
            ),
        )

    return field_check_expr, field_check_details


def sample_sum_sanity_checks(
    t: Union[hl.MatrixTable, hl.Table],
    subsets: List[str],
    pops: List[str] = POPS,
    sexes: List[str] = SEXES,
    subset_pops: Dict[str, List[str]] = {"hgdp": HGDP_POPS, "tgp": TGP_POPS},
    verbose: bool = False,
    metric_first_label: bool = True,
) -> None:
    """
    Performs sanity checks on sample sums in input Table.

    Computes afresh the sum of annotations for a specified group of annotations, and compare to the annotated version;
    displays results from checking the sum of the specified annotations in the terminal.
    Also checks that annotations for all expected sample populations are present.

    :param t: Input Table.
    :param subsets: List of sample subsets.
    :param pops: List of pops in table.
    :param sexes: List of sexes in table.
    :param subset_pops: Dict with subset (keys) and populations within subset (values).
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param metric_first_label: If True, metric precedes label group, e.g. AC-afr-male. If False, label group precedes metric, afr-male-AC.
    :return: None
    """
    t = t.rows() if isinstance(t, hl.MatrixTable) else t

    # Add "" for sum checks on entire callset
    subsets.append("")
    field_check_expr = {}
    field_check_details = {}
    for subset in subsets:
        pop_names = pops
        if subset_pops and subset in subset_pops:
            pop_names = subset_pops[subset]

        # We do not store the raw callstats for the pop or sex groupings of any subset so only check adj here.
        field_check_expr_s, field_check_details_s = sample_sum_check(
            t, subset, metric_first_label, dict(group=["adj"], pop=pop_names), verbose,
        )
        field_check_expr.update(field_check_expr_s)
        field_check_details.update(field_check_details_s)
        field_check_expr_s, field_check_details_s = sample_sum_check(
            t, subset, metric_first_label, dict(group=["adj"], sex=sexes), verbose
        )
        field_check_expr.update(field_check_expr_s)
        field_check_details.update(field_check_details_s)
        field_check_expr_s, field_check_details_s = sample_sum_check(
            t,
            subset,
            metric_first_label,
            dict(group=["adj"], pop=list(set(pop_names)), sex=sexes),
            verbose,
        )
        field_check_expr.update(field_check_expr_s)
        field_check_details.update(field_check_details_s)

        generic_field_check_loop(t, field_check_expr, field_check_details, verbose)


def compare_row_counts(ht1: hl.Table, ht2: hl.Table) -> bool:
    """
    Check if the row counts in two Tables are the same.

    :param ht1: First Table to be checked
    :param ht2: Second Table to be checked
    :return: Whether the row counts are the same
    """
    r_count1 = ht1.count()
    r_count2 = ht2.count()
    logger.info("%d rows in left table; %d rows in right table", r_count1, r_count2)
    return r_count1 == r_count2


def summarize_variants(
    t: Union[hl.MatrixTable, hl.Table],
    monoallelic_expr: Optional[hl.expr.BooleanExpression] = None,
) -> hl.Struct:
    """
    Get summary of variants in a MatrixTable or Table.

    Print the number of variants to stdout and check that each chromosome has variant calls.

    :param t: Input MatrixTable or Table to be checked.
    :param monoallelic_expr: Boolean expression to log how many monoallelic sites are in the Table.
    :return: Struct of variant summary
    """

    if isinstance(t, hl.MatrixTable):
        logger.info("Dataset has %d samples.", t.count_cols())

    var_summary = hl.summarize_variants(t, show=False)
    logger.info(
        "Dataset has %d variants distributed across the following contigs: %s",
        var_summary.n_variants,
        var_summary.contigs,
    )

    for contig in var_summary.contigs:
        if var_summary.contigs[contig] == 0:
            logger.warning("%s has no variants called", var_summary.contigs)

    if monoallelic_expr is not None:
        if isinstance(t, hl.MatrixTable):
            mono_sites = t.filter_rows(monoallelic_expr).count_rows()
        else:
            mono_sites = t.filter(monoallelic_expr).count()
        logger.info("There are %d monoallelic sites in the dataset.", mono_sites)

    return var_summary


def histograms_sanity_check(
    t: Union[hl.MatrixTable, hl.Table], verbose: bool, hists: List[str] = HISTS
) -> None:
    """
    Check the number of that variants that have nonzero values, with the excepion of DP hists, in their n_smaller and n_larger bins of quality histograms for both raw and adj.

    Histogram annotations must exist within an info struct. For example, check that t.info.dp_hist_all_n_smaller != 0. 
    All n_smaller and n_larger annotations must be within an info struct annotation. 

    :param t: Input MatrixTable or Table.
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param hists: List of variant annotation histograms.
    :return: None
    """
    t = t.rows() if isinstance(t, hl.MatrixTable) else t
    field_check_expr = {}
    field_check_details = {}

    for hist in hists:
        for suffix in ["", "raw"]:
            if suffix == "raw":
                hist = f"{hist}_{suffix}"

            check_field = f"{hist}_n_smaller"
            check_description = f"{check_field} == 0"
            field_check_expr, field_check_details = make_field_check_dicts(
                field_check_expr=field_check_expr,
                field_check_details=field_check_details,
                check_description=check_description,
                cond_expr=t.info[check_field] != 0,
                display_fields=hl.struct(**{check_field: t.info[check_field]}),
            )

            if hist not in {
                "dp_hist_alt",
                "dp_hist_all",
            }:  # NOTE: DP hists can have nonzero values in n_larger bin
                check_field = f"{hist}_n_larger"
                check_description = f"{check_field} == 0"
                field_check_expr, field_check_details = make_field_check_dicts(
                    field_check_expr=field_check_expr,
                    field_check_details=field_check_details,
                    check_description=check_description,
                    cond_expr=t.info[check_field] != 0,
                    display_fields=hl.struct(**{check_field: t.info[check_field]}),
                )
    generic_field_check_loop(
        t, field_check_expr, field_check_details, verbose,
    )


def raw_and_adj_sanity_checks(
    t: Union[hl.MatrixTable, hl.Table],
    subsets: List[str],
    verbose: bool,
    delimiter: str = "-",
    metric_first_label: bool = True,
) -> None:
    """
    Perform sanity checks on raw and adj data in input Table/MatrixTable.

    Check that:
        - Raw AC, AN, AF are not 0
        - Adj AN is not 0 and AC and AF are not negative
        - Raw values for AC, AN, nhomalt in each sample subset are greater than or equal to their corresponding adj values

    Raw and adj call stat annotations must be in an info struct annotation on the Table/MatrixTable, e.g. t.info.AC-raw.

    :param t: Input MatrixTable or Table to check.
    :param subsets: List of sample subsets.
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param delimiter: String to use as delimiter when making group label combinations.
    :param metric_first_label: If True, metric precedes label group, e.g. AC-afr-male. If False, label group precedes metric, afr-male-AC.
    :return: None
    """
    t = t.rows() if isinstance(t, hl.MatrixTable) else t

    field_check_expr = {}
    field_check_details = {}
    for subfield in ["AC", "AF"]:
        # Check raw AC, AF > 0
        check_field = f"{subfield}{delimiter}raw"
        field_check_expr, field_check_details = make_field_check_dicts(
            field_check_expr=field_check_expr,
            field_check_details=field_check_details,
            check_description=f"{check_field} > 0",
            cond_expr=t.info[check_field] <= 0,
            display_fields=hl.struct(**{check_field: t.info[check_field]}),
        )

        check_field = f"{subfield}{delimiter}adj"
        field_check_expr, field_check_details = make_field_check_dicts(
            field_check_expr=field_check_expr,
            field_check_details=field_check_details,
            check_description=f"{check_field} >= 0",
            cond_expr=t.info[check_field] < 0,
            display_fields=hl.struct(
                **{check_field: t.info[check_field], "filters": t.filters}
            ),
        )

    # Check raw AN > 0
    check_field = f"AN{delimiter}raw"
    field_check_expr, field_check_details = make_field_check_dicts(
        field_check_expr=field_check_expr,
        field_check_details=field_check_details,
        check_description=f"{check_field} > 0",
        cond_expr=t.info[check_field] <= 0,
        display_fields=hl.struct(**{check_field: t.info[check_field]}),
    )

    # Check adj AN >= 0
    check_field = f"AN{delimiter}adj"
    field_check_expr, field_check_details = make_field_check_dicts(
        field_check_expr=field_check_expr,
        field_check_details=field_check_details,
        check_description=f"{check_field} >= 0",
        cond_expr=t.info[check_field] < 0,
        display_fields=hl.struct(**{check_field: t.info[check_field]}),
    )
    # Check overall gnomad's raw subfields >= adj
    for subfield in ["AC", "AN", "nhomalt"]:
        check_field_left = f"{subfield}{delimiter}raw"
        check_field_right = f"{subfield}{delimiter}adj"
        field_check_expr, field_check_details = make_field_check_dicts(
            field_check_expr=field_check_expr,
            field_check_details=field_check_details,
            check_description=f"{check_field_left} >= {check_field_right}",
            cond_expr=t.info[check_field_left] < t.info[check_field_right],
            display_fields=hl.struct(
                **{
                    check_field_left: t.info[check_field_left],
                    check_field_right: t.info[check_field_right],
                }
            ),
        )

        for subset in subsets:
            field_check_label = (
                f"{subfield}{delimiter}{subset}{delimiter}"
                if metric_first_label
                else f"{subset}{delimiter}{subfield}{delimiter}"
            )
            check_field_left = f"{field_check_label}raw"
            check_field_right = f"{field_check_label}adj"

            field_check_expr, field_check_details = make_field_check_dicts(
                field_check_expr=field_check_expr,
                field_check_details=field_check_details,
                check_description=f"{check_field_left} >= {check_field_right}",
                cond_expr=t.info[check_field_left] < t.info[check_field_right],
                display_fields=hl.struct(
                    **{
                        check_field_left: t.info[check_field_left],
                        check_field_right: t.info[check_field_right],
                    }
                ),
            )

    generic_field_check_loop(t, field_check_expr, field_check_details, verbose)


def sex_chr_sanity_checks(
    t: Union[hl.MatrixTable, hl.Table],
    info_metrics: List[str],
    contigs: List[str],
    verbose: bool,
    delimiter: str = "-",
) -> None:
    """
    Performs sanity checks for annotations on the sex chromosomes.
    Check:
        - That metrics for chrY variants in XX samples are NA and not 0
        - That nhomalt counts are equal to XX nhomalt counts for all non-PAR chrX variants

    :param t: Input MatrixTable or Table.
    :param info_metrics: List of metrics in info struct of input Table.
    :param contigs: List of contigs present in input Table.
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param delimiter: String to use as the delimiter in XX metrics
    :return: None
    """
    t = t.rows() if isinstance(t, hl.MatrixTable) else t

    xx_metrics = [
        x for x in info_metrics if f"{delimiter}female" in x or f"{delimiter}XX" in x
    ]

    if "chrY" in contigs:
        logger.info("Check values of XX metrics for Y variants are NA:")
        t_y = hl.filter_intervals(t, [hl.parse_locus_interval("chrY")])
        metrics_values = {}
        for metric in xx_metrics:
            metrics_values[metric] = hl.agg.any(hl.is_defined(t_y.info[metric]))
        output = dict(t_y.aggregate(hl.struct(**metrics_values)))
        for metric, value in output.items():
            if value:
                values_found = t_y.aggregate(
                    hl.agg.filter(
                        hl.is_defined(t_y.info[metric]),
                        hl.agg.take(t_y.info[metric], 1),
                    )
                )
                logger.info(
                    "FAILED %s = %s check for Y variants. Values found: %s",
                    metric,
                    None,
                    values_found,
                )
            else:
                logger.info("PASSED %s = %s check for Y variants", metric, None)

    t_x = hl.filter_intervals(t, [hl.parse_locus_interval("chrX")])
    t_xnonpar = t_x.filter(t_x.locus.in_x_nonpar())
    n = t_xnonpar.count()
    logger.info("Found %d X nonpar sites", n)

    logger.info("Check (nhomalt == nhomalt_xx) for X nonpar variants:")
    xx_metrics = [x for x in xx_metrics if "nhomalt" in x]

    field_check_expr = {}
    field_check_details = {}
    for metric in xx_metrics:
        standard_field = metric.replace(f"{delimiter}female", "").replace(
            f"{delimiter}XX", ""
        )
        check_field_left = f"{metric}"
        check_field_right = f"{standard_field}"
        field_check_expr, field_check_details = make_field_check_dicts(
            field_check_expr=field_check_expr,
            field_check_details=field_check_details,
            check_description=f"{check_field_left} == {check_field_right}",
            cond_expr=t_xnonpar.info[check_field_left]
            != t_xnonpar.info[check_field_right],
            display_fields=hl.struct(
                **{
                    check_field_left: t_xnonpar.info[check_field_left],
                    check_field_right: t_xnonpar.info[check_field_right],
                }
            ),
        )

    generic_field_check_loop(t_xnonpar, field_check_expr, field_check_details, verbose)


def missingness_sanity_checks(
    t: Union[hl.MatrixTable, hl.Table],
    info_metrics: List[str],
    non_info_metrics: List[str],
    n_sites: int,
    missingness_threshold: float,
) -> None:
    """
    Check amount of missingness in all row annotations.
    
    Print metric to sdout if the metric annotations missingness exceeds the missingness_threshold.

    :param t: Input MatrixTable or Table.
    :param info_metrics: List of metrics in info struct of input Table.
    :param non_info_metrics: List of row annotations minus info struct from input Table.
    :param n_sites: Number of sites in input Table.
    :param missingness_threshold: Upper cutoff for allowed amount of missingness.
    :return: None
    """
    t = t.rows() if isinstance(t, hl.MatrixTable) else t

    logger.info(
        "Missingness threshold (upper cutoff for what is allowed for missingness checks): %d",
        missingness_threshold,
    )
    metrics_frac_missing = {}
    for x in info_metrics:
        metrics_frac_missing[x] = hl.agg.sum(hl.is_missing(t.info[x])) / n_sites
    for x in non_info_metrics:
        metrics_frac_missing[x] = hl.agg.sum(hl.is_missing(t[x])) / n_sites
    output = t.aggregate(hl.struct(**metrics_frac_missing))

    n_fail = 0
    for metric, value in dict(output).items():
        message = "missingness check for %s: %d% missing", metric, 100 * value
        if value > missingness_threshold:
            logger.info("FAILED %s", message)
            n_fail += 1
        else:
            logger.info("Passed %s", message)
    logger.info("%d missing metrics checks failed", n_fail)


def vcf_field_check(
    t: Union[hl.MatrixTable, hl.Table],
    header_dict: Dict[str, Dict[str, Dict[str, str]]],
    row_annotations: List[str],
    hists: List[str] = HISTS,
) -> bool:
    """
    Check that all VCF fields and descriptions are present in input Table and VCF header dictionary.

    :param t: Input MatrixTable or Tableto be exported to VCF.
    :param header_dict: VCF header dictionary.
    :param row_annotations: List of row annotations in MatrixTable.
    :param hists: List of variant histogram annotations. Default is HISTS.
    :return: Boolean with whether all expected fields and descriptions are present.
    """
    t = t.rows() if isinstance(t, hl.MatrixTable) else t

    # Confirm all VCF fields/descriptions are present before exporting
    hist_fields = []
    for hist in hists:
        hist_fields.extend(
            [
                f"{hist}_bin_freq",
                f"{hist}_n_smaller",
                f"{hist}_n_larger",
                f"{hist}_raw_bin_freq",
                f"{hist}_raw_n_smaller",
                f"{hist}_raw_n_larger",
            ]
        )

    missing_fields = []
    missing_descriptions = []
    for item in ["info", "filter"]:
        if item == "info":
            annots = row_annotations
        else:
            annot_t = t.explode(t.filters)
            annots = list(annot_t.aggregate(hl.agg.collect_as_set(annot_t.filters)))

        temp_missing_fields = []
        temp_missing_descriptions = []
        for field in annots:
            try:
                description = header_dict[item][field]
                if len(description) == 0:
                    logger.warning(
                        "%s in T info field has empty description in VCF header!", field
                    )
                    temp_missing_descriptions.append(field)
            except KeyError:
                logger.warning(
                    "%s in T info field does not exist in VCF header!", field
                )
                # NOTE: some hists are not exported, so ignoring here
                # END entry is also not exported (removed during densify)
                if (field not in hist_fields) and (field != "END"):
                    temp_missing_fields.append(field)

        missing_fields.extend(temp_missing_fields)
        missing_descriptions.extend(temp_missing_descriptions)

    if len(missing_fields) != 0 or len(missing_descriptions) != 0:
        logger.error(
            "Some fields are either missing or missing descriptions in the VCF header! Please reconcile."
        )
        logger.error("Missing fields: %s", missing_fields)
        logger.error("Missing descriptions: %s", missing_descriptions)
        return False

    logger.info("Passed VCF fields check!")
    return True


def sanity_check_release_t(
    t: Union[hl.MatrixTable, hl.Table],
    subsets: List[str],
    missingness_threshold: float = 0.5,
    monoallelic_check: bool = True,
    verbose: bool = True,
    show_percent_sites: bool = True,
    delimiter: str = "-",
    metric_first_label: bool = True,
    hists: List[str] = HISTS,
    pops: List[str] = POPS,
    sexes: List[str] = SEXES,
    subset_pops: Dict[str, List[str]] = {"hgdp": HGDP_POPS, "tgp": TGP_POPS},
    summarize_variants_check: bool = True,
    filters_check: bool = True,
    histograms_check: bool = True,
    raw_adj_check: bool = True,
    subset_freq_check: bool = True,
    samples_sum_check: bool = True,
    sex_chr_check: bool = True,
    missingness_check: bool = True,
) -> None:
    """
    Perform a battery of sanity checks on a specified group of subsets in a MatrixTable containing variant annotations.

    Includes:
    - Summaries of % filter status for different partitions of variants
    - Histogram outlier bin checks
    - Checks on AC, AN, and AF annotations
    - Checks that subgroup annotation values add up to the supergroup annotation values
    - Checks on sex-chromosome annotations; and summaries of % missingness in variant annotations

    All annotations must be within an info struct, e.g. t.info.AC-raw.

    :param t: Input MatrixTable or Table containing variant annotations to check.
    :param subsets: List of subsets to be checked.
    :param missingness_threshold: Upper cutoff for allowed amount of missingness. Default is 0.5.
    :param monoallelic_check: Log how many monoallelic sites are in the Table; requires a monoallelic annotation within an info struct.
    :param verbose: If True, display top values of relevant annotations being checked, regardless of whether check
        conditions are violated; if False, display only top values of relevant annotations if check conditions are violated.
    :param show_percent_sites: Show percentage of sites that fail checks. Default is False.
    :param metric_first_label: If True, metric precedes label group, e.g. AC-afr-male. If False, label group precedes metric, afr-male-AC.
    :param hists: List of variant annotation histograms.
    :param pops: List of pops in table.
    :param sexes: List of sexes in table.
    :param subset_pops: Dict with subset (keys) and populations within subset (values).
    :param summarize_variants_check: When true, runs the summarize_variants method.
    :param filters_check: When true, runs the filters_sanity_check method.
    :param histograms_check: When true, runs the histograms_sanity_check method.
    :param raw_adj_check: When true, runs the raw_and_adj_sanity_checks method.
    :param subset_freq_check: When true, runs the subset_freq_sanity_checks method.
    :param samples_sum_check: When true, runs the sample_sum_sanity_checks method.
    :param sex_chr_check: When true, runs the sex_chr_sanity_checks method.
    :param missingness_check: When true, runs the missingness_sanity_checks method.
    :return: None (terminal display of results from the battery of sanity checks).
    """

    # Perform basic checks -- number of variants, number of contigs, number of samples
    if summarize_variants_check:
        logger.info("BASIC SUMMARY OF INPUT TABLE:")
        summarize_variants(t, monoallelic_check)

    if filters_check:
        logger.info("VARIANT FILTER SUMMARIES:")
        filters_sanity_check(t)

    if histograms_check:
        logger.info("HISTOGRAM CHECKS:")
        histograms_sanity_check(t, verbose, hists)

    if raw_adj_check:
        logger.info("RAW AND ADJ CHECKS:")
        raw_and_adj_sanity_checks(t, subsets, verbose, delimiter, metric_first_label)

    if subset_freq_check:
        logger.info("SUBSET FREQUENCY CHECKS:")
        subset_freq_sanity_checks(
            t, subsets, verbose, show_percent_sites, delimiter, metric_first_label
        )

    if samples_sum_check:
        logger.info("SAMPLE SUM CHECKS:")
        sample_sum_sanity_checks(
            t, subsets, pops, sexes, subset_pops, verbose, metric_first_label
        )

    info_metrics = list(t.row.info)

    if sex_chr_check:
        logger.info("SEX CHROMOSOME ANNOTATION CHECKS:")
        contigs = t.aggregate(hl.agg.collect_as_set(t.locus.contig))
        sex_chr_sanity_checks(t, info_metrics, contigs, verbose, delimiter)

    if missingness_check:
        logger.info("MISSINGNESS CHECKS:")
        non_info_metrics = list(t.row)
        non_info_metrics.remove("info")
        n_sites = t.count()
        missingness_sanity_checks(
            t, info_metrics, non_info_metrics, n_sites, missingness_threshold
        )
