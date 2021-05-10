# noqa: D100

import logging
from typing import Dict, List, Optional, Union

import hail as hl

from gnomad.utils.vcf import HISTS, make_label_combos, SORT_ORDER


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def generic_field_check(
    ht: hl.Table,
    cond_expr: hl.expr.BooleanExpression,
    check_description: str,
    display_fields: List[str],
    verbose: bool,
    show_percent_sites: bool = False,
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
    :return: None
    """
    ht_orig = ht
    ht = ht.filter(cond_expr)
    n_fail = ht.count()
    if n_fail > 0:
        logger.info("Found %d sites that fail %s check:", n_fail, check_description)
        if show_percent_sites:
            logger.info("Percentage of sites that fail: %f", n_fail / ht_orig.count())
        ht = ht.flatten()
        ht.select("locus", "alleles", *display_fields).show()
    else:
        logger.info("PASSED %s check", check_description)
        if verbose:
            ht_orig = ht_orig.flatten()
            ht_orig.select(*display_fields).show()


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
    ht: hl.Table,
    prefix: str,
    label_groups: Dict[str, List[str]],
    verbose: bool,
    subpop: bool = None,
    sort_order: List[str] = SORT_ORDER,
) -> None:
    """
    Compute the sum of annotations for a specified group of annotations, and compare to the annotated version.

    Results from checking the sum of the specified annotations are printed in the terminal.

    :param ht: Table containing annotations to be summed.
    :param prefix: String indicating sample subset.
    :param label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"]).
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param subpop: Subpop abbreviation, supplied only if subpopulations are included in the annotation groups being checked.
    :param sort_order: List containing order to sort label group combinations. Default is SORT_ORDER.
    :return: None
    """
    if prefix != "":
        prefix = f"{prefix}_"

    label_combos = make_label_combos(label_groups)
    combo_AC = [ht.info[f"{prefix}AC_{x}"] for x in label_combos]
    combo_AN = [ht.info[f"{prefix}AN_{x}"] for x in label_combos]
    combo_nhomalt = [ht.info[f"{prefix}nhomalt_{x}"] for x in label_combos]

    group = label_groups.pop("group")[0]
    alt_groups = "_".join(
        sorted(label_groups.keys(), key=lambda x: sort_order.index(x))
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
                    ht.info[f"{prefix}{subfield}_{subpop}_{group}"]
                    != ht[f"sum_{subfield}_{group}_{alt_groups}"]
                ),
                f"{prefix}{subfield}_{subpop}_{group} = sum({subfield}_{group}_{alt_groups})",
                [
                    f"info.{prefix}{subfield}_{subpop}_{group}",
                    f"sum_{subfield}_{group}_{alt_groups}",
                ],
                verbose,
            )


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
    t: Union[hl.MatrixTable, hl.Table], monoallelic_expr: Optional[hl.expr.BooleanExpression] = None
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
        "Dataset has %d variants distributed across the following contigs: %s", var_summary.n_variants, var_summary.contigs
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
    Check that variants have nonzero values, with the excepion of DP hists, in their n_smaller and n_larger bins of quality histograms for both raw and adj.

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

    for hist in hists:
        for suffix in ["", "raw"]:
            if suffix == "raw":
                logger.info("Checking raw qual hists...")
                hist = f"{hist}_{suffix}"
            else:
                logger.info("Checking adj qual hists...")

            generic_field_check(
                t,
                cond_expr=(t.info[f"{hist}_n_smaller"] != 0),
                check_description=f"{hist}_n_smaller == 0",
                display_fields=[f"info.{hist}_n_smaller"],
                verbose=verbose,
            )
            if hist not in {
                "dp_hist_alt",
                "dp_hist_all",
            }:  # NOTE: DP hists can have nonzero values in n_larger bin
                generic_field_check(
                    t,
                    cond_expr=(t.info[f"{hist}_n_larger"] != 0),
                    check_description=f"{hist}_n_larger == 0",
                    display_fields=[f"info.{hist}_n_larger"],
                    verbose=verbose,
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

