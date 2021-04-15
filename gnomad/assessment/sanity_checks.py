import logging
from typing import Dict, List, Optional

import hail as hl

from gnomad.utils.vcf import make_label_combos, SORT_ORDER


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
        logger.info(f"Found {n_fail} sites that fail {check_description} check:")
        if show_percent_sites:
            logger.info(f"Percentage of sites that fail: {n_fail / ht_orig.count()}")
        ht = ht.flatten()
        ht.select("locus", "alleles", *display_fields).show()
    else:
        logger.info(f"PASSED {check_description} check")
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
    Compute afresh the sum of annotations for a specified group of annotations, and compare to the annotated version;
    display results from checking the sum of the specified annotations in the terminal.

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
    logger.info(f"{r_count1} rows in left table; {r_count2} rows in right table")
    return r_count1 == r_count2
