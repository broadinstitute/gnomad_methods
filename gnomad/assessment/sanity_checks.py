# noqa: D100

import logging
from typing import Dict, List, Optional, Tuple, Union

import hail as hl

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
    ht: hl.Table, field_check_expr: dict, field_check_details: dict, verbose: bool
):
    """
    Loop through all conditional checks for a given hail Table.

    This loop allows aggregation across the hail Table once, as opposed to aggregating during every conditional check. 
    :param ht: Table containing annotations to be checked.
    :param field_check_expr: Dictionary whose keys are conditions being checked and values are the expressions for filtering to condition.
    :param field_check_details: Dictionary whose keys are the the check descriptions and values are struct expressions used for check and what to display for the check in the terminal.
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
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
    ht: hl.Table, extra_filter_checks: Optional[Dict[str, hl.expr.Expression]] = None
) -> Dict[str, hl.expr.Expression]:
    """
    Make Hail expressions to measure % variants filtered under varying conditions of interest.

    Checks for:
        - Total number of variants
        - Fraction of variants removed due to:
            - Any filter
            - Inbreeding coefficient filter in combination with any other filter
            - AC0 filter in combination with any other filter
            - Random forest filtering in combination with any other filter
            - Only inbreeding coefficient filter
            - Only AC0 filter
            - Only random forest filtering

    :param ht: Table containing 'filter' annotation to be examined.
    :param extra_filter_checks: Optional dictionary containing filter condition name (key) extra filter expressions (value) to be examined.
    :return: Dictionary containing Hail aggregation expressions to examine filter flags.
    """
    filters_dict = {
        "n": hl.agg.count(),
        "frac_any_filter": hl.agg.fraction(hl.len(ht.filters) != 0),
        "frac_inbreed_coeff": hl.agg.fraction(ht.filters.contains("InbreedingCoeff")),
        "frac_ac0": hl.agg.fraction(ht.filters.contains("AC0")),
        "frac_rf": hl.agg.fraction(ht.filters.contains("RF")),
        "frac_inbreed_coeff_only": hl.agg.fraction(
            ht.filters.contains("InbreedingCoeff") & (ht.filters.length() == 1)
        ),
        "frac_ac0_only": hl.agg.fraction(
            ht.filters.contains("AC0") & (ht.filters.length() == 1)
        ),
        "frac_rf_only": hl.agg.fraction(
            ht.filters.contains("RF") & (ht.filters.length() == 1)
        ),
    }
    if extra_filter_checks:
        filters_dict.update(extra_filter_checks)

    return filters_dict


def sample_sum_check(
    t: Union[hl.MatrixTable, hl.Table],
    subset: str,
    label_groups: Dict[str, List[str]],
    verbose: bool,
    subpop: bool = None,
    sort_order: List[str] = SORT_ORDER,
    delimiter: str = "-",
    metric_first_label: bool = True,
) -> None:
    """
    Compute the sum of call stats annotations for a specified group of annotations, compare to the annotated version, and display the result in stdout.

    For example, if pop1 consists of pop1, pop2, and pop3, check that t.info.AC-subset1 == sum(t.info.AC-subset1-pop1, t.info.AC-subset1-pop2, t.info.AC-subset1-pop3).

    :param t: Input MatrixTable or Table containing annotations to be summed.
    :param subset: String indicating sample subset.
    :param label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["XY", "XX"] or ["afr", "nfe", "amr"]).
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param subpop: Subpop abbreviation, supplied only if subpopulations are included in the annotation groups being checked.
    :param sort_order: List containing order to sort label group combinations. Default is SORT_ORDER.
    :param delimiter: String to use as delimiter when making group label combinations.
    :param metric_first_label: If True, metric precedes label group, e.g. AC-afr-male. If False, label group precedes metric, afr-male-AC.
    :return: Tuple of dictionaries
    """

    t = t.rows() if isinstance(t, hl.MatrixTable) else t

    if subset:  # TODO: add & subset != "":
        subset += delimiter

    label_combos = make_label_combos(label_groups, label_delimiter=delimiter)
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
    :param hl.Table ht: Input Table.
    :param List[str] subsets: List of sample subsets.
    :param bool verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param pop_names: Dict with global population names (keys) and population descriptions (values).
    :param metric_first_label: If True, metric precedes label group, e.g. AC-afr-male. If False, label group precedes metric, afr-male-AC.
    :return: None
    :rtype: None
    """
    t = t.rows() if isinstance(t, hl.MatrixTable) else t
    # Add "" for sum checks on entire callset
    subsets.append("")
    # Perform sample sum checks per subset
    field_check_expr = {}
    field_check_details = {}
    for subset in subsets:
        pop_names = pops
        if subset_pops and subset in subset_pops:
            pop_names = subset_pops[subset]

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
    :rtype: Struct
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
    :rtype: None
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
    :return: None
    :rtype: None
    """
    t = t.rows() if isinstance(t, hl.MatrixTable) else t

    for field in ["AC", "AF"]:
        field = f"{field}{delimiter}"
        # Check raw AC, AF > 0
        generic_field_check(
            t,
            cond_expr=(t.info[f"{field}raw"] <= 0),
            check_description=f"{field}raw > 0",
            display_fields=[f"info.{field}raw"],
            verbose=verbose,
        )
        # Check adj AC, AF >=0
        generic_field_check(
            t,
            cond_expr=(t.info[f"{field}adj"] < 0),
            check_description=f"{field}adj >= 0",
            display_fields=[f"info.{field}adj", "filters"],
            verbose=verbose,
        )

    # Check raw AN > 0
    an_raw_field = f"AN{delimiter}raw"
    generic_field_check(
        t,
        cond_expr=(t.info[an_raw_field] <= 0),
        check_description=f"{an_raw_field} > 0",
        display_fields=[f"info.{an_raw_field}"],
        verbose=verbose,
    )

    an_adj_field = f"AN{delimiter}adj"
    # Check adj AN >= 0
    generic_field_check(
        t,
        cond_expr=(t.info[an_adj_field] < 0),
        check_description=f"{an_adj_field} >= 0",
        display_fields=[f"info.{an_adj_field}"],
        verbose=verbose,
    )

    # Check overall raw subfields >= adj
    for field in ["AC", "AN", "nhomalt"]:
        field = f"{field}{delimiter}"
        generic_field_check(
            t,
            cond_expr=(t.info[f"{field}raw"] < t.info[f"{field}adj"]),
            check_description=f"{field}raw >= {field}adj",
            display_fields=[f"info.{field}raw", f"info.{field}adj",],
            verbose=verbose,
        )

    for subset in subsets:
        for field in ["AC", "AN", "nhomalt"]:
            # Check AC_raw >= AC adj
            field = f"{field}{delimiter}{subset}{delimiter}"
            generic_field_check(
                t,
                cond_expr=(t.info[f"{field}raw"] < t.info[f"{field}adj"]),
                check_description=f"{field}raw >= {field}adj",
                display_fields=[f"info.{field}raw", f"info.{field}adj",],
                verbose=verbose,
            )
