"""Tests for the validity_checks module."""

import hail as hl
import pytest

from gnomad.assessment.validity_checks import (
    check_array_struct_missingness,
    check_missingness_of_struct,
    compute_and_check_summations,
    flatten_missingness_struct,
    unfurl_array_annotations,
)


@pytest.fixture
def ht_for_check_missingness_of_struct() -> hl.Table:
    """Fixture to set up a Hail Table with the desired nested structure and data."""
    # Create test data with nested structs.
    data = [
        {
            "idx": 0,
            "s": hl.struct(
                a=1,
                b="value1",
                c=hl.struct(
                    d=[hl.missing(hl.tstr), hl.missing(hl.tstr)],
                    e="test1",
                    f={"v1", "v2"},
                ),
            ),
        },
        {
            "idx": 1,
            "s": hl.struct(
                a=2,
                b="value2",
                c=hl.struct(
                    d=["not missing", hl.missing(hl.tstr)],
                    e=hl.missing(hl.tstr),
                    f={"v3", hl.missing(hl.tstr)},
                ),
            ),
        },
        {
            "idx": 2,
            "s": hl.struct(
                a=3,
                b=hl.missing(hl.tstr),
                c=hl.struct(
                    d=hl.missing(hl.tarray(hl.tstr)),
                    e=hl.missing(hl.tstr),
                    f=hl.empty_set(hl.tstr),
                ),
            ),
        },
        {
            "idx": 3,
            "s": hl.struct(
                a=4,
                b="value3",
                c=hl.struct(d=["foo", "bar"], e="test2", f=hl.empty_set(hl.tstr)),
            ),
        },
        {
            "idx": 4,
            "s": hl.struct(
                a=5,
                b="value4",
                c=hl.struct(
                    d=hl.empty_array(hl.tstr), e="test3", f=hl.empty_set(hl.tstr)
                ),
            ),
        },
    ]

    # Convert data into a Hail table.
    ht = hl.Table.parallelize(
        data,
        hl.tstruct(
            idx=hl.tint32,
            s=hl.tstruct(
                a=hl.tint32,
                b=hl.tstr,
                c=hl.tstruct(d=hl.tarray(hl.tstr), e=hl.tstr, f=hl.tset(hl.tstr)),
            ),
        ),
    )

    return ht


def test_check_missingness_of_struct(
    ht_for_check_missingness_of_struct: hl.Table,
) -> None:
    """Test check_missingness_of_struct and flatten results."""
    ht = ht_for_check_missingness_of_struct

    # Apply check_missingness_of_struct aggregate results.
    metric_missingness = check_missingness_of_struct(ht.s, "s")
    missingness_struct = ht.aggregate(hl.struct(**metric_missingness))

    # Flatten the result struct into a dictionary.
    missingness_dict = flatten_missingness_struct(missingness_struct)

    # Define expected missingness results.
    expected_result = {
        "s.a": 0.0,
        "s.b": 0.2,
        "s.c.d": 0.6,
        "s.c.e": 0.4,
        "s.c.f": 0.6,
    }

    # Compare the results with the expected values.
    for key, expected_value in expected_result.items():
        assert missingness_dict[key] == expected_value, (
            f"Mismatch for {key}: "
            f"expected {expected_value}, got {missingness_dict[key]}"
        )


@pytest.fixture()
def ht_for_check_array_struct_missingness() -> hl.Table:
    """Fixture to set up a sample Hail Table with array<struct> fields and global index dictionary ."""
    # Define rows of the Table.
    data = [
        {
            "idx": 0,
            "freq": [
                {"AC": 5, "AF": 0.1, "AN_eas": 20, "AN_sas": 3},
                {"AC": 10, "AF": 0.05, "AN_eas": 5, "AN_sas": None},
            ],
        },
        {
            "idx": 1,
            "freq": [
                {"AC": 6, "AF": 0.08, "AN_eas": None, "AN_sas": 4},
                {"AC": 8, "AF": 0.50, "AN_eas": None, "AN_sas": None},
            ],
        },
        {
            "idx": 2,
            "freq": [
                {"AC": 65, "AF": 0.18, "AN_eas": None, "AN_sas": 2},
                {"AC": 88, "AF": 0.20, "AN_eas": None, "AN_sas": None},
            ],
        },
        {"idx": 3, "freq": [{"AC": 8, "AF": 0.08, "AN_eas": 16, "AN_sas": 2}, None]},
    ]

    # Create Table.
    ht = hl.Table.parallelize(
        data,
        hl.tstruct(
            idx=hl.tint32,
            freq=hl.tarray(
                hl.tstruct(
                    AC=hl.tint32, AF=hl.tfloat64, AN_eas=hl.tint32, AN_sas=hl.tint32
                )
            ),
        ),
    )

    # Define global annotation for freq_index_dict.
    freq_index_dict = {"adj": 0, "raw": 1}
    ht = ht.annotate_globals(freq_index_dict=freq_index_dict)

    # Unfurl indexed array annotations.
    annotations = unfurl_array_annotations(ht, {"freq": "freq_index_dict"})
    ht = ht.annotate(**annotations)

    return ht


def test_unfurl_array_annotations(
    ht_for_check_array_struct_missingness: hl.Table,
) -> None:
    """Test the unfurl_array_annotations function for all rows."""
    ht = ht_for_check_array_struct_missingness
    indexed_array_annotations = {"freq": "freq_index_dict"}

    # Call the unfurl_array_annotations function.
    result = unfurl_array_annotations(ht, indexed_array_annotations)

    # Define names of the expected keys after unfurling.
    expected_keys = {
        "AC_adj",
        "AF_adj",
        "AN_eas_adj",
        "AN_sas_adj",
        "AC_raw",
        "AF_raw",
        "AN_eas_raw",
        "AN_sas_raw",
    }
    assert (
        set(result.keys()) == expected_keys
    ), "Unfurled keys do not match expected keys."

    # Annotate table with unfurled fields.
    ht = ht.annotate(**result)
    rows = ht.collect()

    # Define the expected values for each unfurled annotation.
    expected_values = [
        {
            "AC_adj": 5,
            "AF_adj": 0.1,
            "AN_eas_adj": 20,
            "AN_sas_adj": 3,
            "AC_raw": 10,
            "AF_raw": 0.05,
            "AN_eas_raw": 5,
            "AN_sas_raw": None,
        },
        {
            "AC_adj": 6,
            "AF_adj": 0.08,
            "AN_eas_adj": None,
            "AN_sas_adj": 4,
            "AC_raw": 8,
            "AF_raw": 0.50,
            "AN_eas_raw": None,
            "AN_sas_raw": None,
        },
        {
            "AC_adj": 65,
            "AF_adj": 0.18,
            "AN_eas_adj": None,
            "AN_sas_adj": 2,
            "AC_raw": 88,
            "AF_raw": 0.20,
            "AN_eas_raw": None,
            "AN_sas_raw": None,
        },
        {
            "AC_adj": 8,
            "AF_adj": 0.08,
            "AN_eas_adj": 16,
            "AN_sas_adj": 2,
            "AC_raw": None,
            "AF_raw": None,
            "AN_eas_raw": None,
            "AN_sas_raw": None,
        },
    ]

    # Validate each expected value.
    for row, expected in zip(rows, expected_values):
        for key, expected_value in expected.items():
            assert row[key] == expected_value, (
                f"Mismatch in row {row['idx']} for key '{key}': "
                f"expected {expected_value}, got {row[key]}"
            )


def test_check_array_struct_missingness(
    ht_for_check_array_struct_missingness: hl.Table,
) -> None:
    """Test the check_array_struct_missingness function for all fields."""
    ht = ht_for_check_array_struct_missingness
    indexed_array_annotations = {"freq": "freq_index_dict"}

    # Call the check_array_struct_missingness function.
    missingness = check_array_struct_missingness(ht, indexed_array_annotations)

    # Define the expected missingness percentages for each unfurled field.
    # All 'adj' values have no missing data.
    # AC_raw is missing only in row 4.
    # AF_raw is missing in rows 2, 3, and 4.
    # AN_raw is missing in all rows.
    expected_missingness = {
        "AC_adj": 0.0,
        "AF_adj": 0.0,
        "AC_raw": 0.25,
        "AF_raw": 0.25,
        "AN_eas_adj": 0.50,
        "AN_eas_raw": 0.75,
        "AN_sas_adj": 0.0,
        "AN_sas_raw": 1.00,
    }

    # Validate each field's missingness percentage.
    for field, expected_value in expected_missingness.items():
        assert missingness[field] == expected_value, (
            f"Mismatch in missingness for field '{field}': "
            f"expected {expected_value}, got {missingness[field]}"
        )


@pytest.fixture
def ht_for_compute_and_check_summations() -> hl.Table:
    """Fixture to set up a Hail Table with the desired structure and data for testing compute_and_check_summations."""
    data = [
        {
            "idx": 0,
            "AC_afr_adj": 5,
            "AC_amr_adj": 10,
            "AC_adj": 15,
            "AN_afr_XX_adj": 20,
            "AN_afr_XY_adj": 30,
            "AN_adj": 50,
        },
        {
            "idx": 1,
            "AC_afr_adj": 3,
            "AC_amr_adj": 7,
            "AC_adj": 10,
            "AN_afr_XX_adj": 15,
            "AN_afr_XY_adj": 25,
            "AN_adj": 40,
        },
        {
            "idx": 2,
            "AC_afr_adj": 2,
            "AC_amr_adj": 3,
            "AC_adj": 6,  # This should cause a mismatch
            "AN_afr_XX_adj": 10,
            "AN_afr_XY_adj": 20,
            "AN_adj": 35,  # This should cause a mismatch
        },
    ]

    ht = hl.Table.parallelize(
        data,
        hl.tstruct(
            idx=hl.tint32,
            AC_afr_adj=hl.tint32,
            AC_amr_adj=hl.tint32,
            AC_adj=hl.tint32,
            AN_afr_XX_adj=hl.tint32,
            AN_afr_XY_adj=hl.tint32,
            AN_adj=hl.tint32,
        ),
    )

    return ht


def test_compute_and_check_summations(
    ht_for_compute_and_check_summations: hl.Table,
) -> None:
    """Test compute_and_check_summations function."""
    ht = ht_for_compute_and_check_summations

    comparison_groups = {
        "AC_group_adj_gen_anc": {
            "values_to_sum": ["AC_afr_adj", "AC_amr_adj"],
            "expected_total": "AC_adj",
        },
        "AN_group_adj_gen_anc_sex": {
            "values_to_sum": ["AN_afr_XX_adj", "AN_afr_XY_adj"],
            "expected_total": "AN_adj",
        },
    }

    result = compute_and_check_summations(ht, comparison_groups)

    expected_result = {
        "AC_group_adj_gen_anc": 1,  # One mismatch in row 2
        "AN_group_adj_gen_anc_sex": 1,  # One mismatch in row 2
    }

    mismatches = {
        key: (result[key], expected_result[key])
        for key in expected_result
        if result.get(key) != expected_result[key]
    }

    assert (
        not mismatches
    ), f"Mismatches found: {', '.join(f'{key} (expected {exp}, got {result})' for key, (result, exp) in mismatches.items())}"
