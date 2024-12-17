"""Tests for the validity_checks module."""
import pytest
import hail as hl
from gnomad.assessment.validity_checks import check_array_struct_missingness,check_missingness_of_struct, flatten_missingness_struct, unfurl_array_annotations


@pytest.fixture
def ht_for_check_missingness_of_struct():
    """
    Fixture to set up a Hail Table with the desired nested structure and data.
    """
    # Create test data with nested structs.
    data = [
        hl.struct(idx=0, s=hl.struct(a=1, b="value1", c=hl.struct(d=[hl.missing(hl.tstr), hl.missing(hl.tstr)], e="test1"))),
        hl.struct(idx=0, s=hl.struct(a=2, b="value2", c=hl.struct(d=["not missing", hl.missing(hl.tstr)], e=hl.missing(hl.tstr)))),
        hl.struct(idx=0, s=hl.struct(a=3, b=hl.missing(hl.tstr), c=hl.struct(d=hl.missing(hl.tarray(hl.tstr)), e=hl.missing(hl.tstr)))),
        hl.struct(idx=0, s=hl.struct(a=4, b="value3", c=hl.struct(d=["foo", "bar"], e="test2"))),
        hl.struct(idx=0, s=hl.struct(a=5, b="value4", c=hl.struct(d=hl.empty_array(hl.tstr), e="test3"))),
    ]

    # Convert data into a Hail table.
    ht = hl.Table.parallelize(
        data,
        hl.tstruct(
            idx=hl.tint32,
            s=hl.tstruct(
                a=hl.tint32,
                b=hl.tstr,
                c=hl.tstruct(
                    d=hl.tarray(hl.tstr),
                    e=hl.tstr
                )
            )
        )
    )

    return ht


def test_check_missingness_of_struct(ht_for_check_missingness_of_struct):
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
    }

    # Compare the results with the expected values.
    for key, expected_value in expected_result.items():
        assert missingness_dict[key] == expected_value


@pytest.fixture()
def ht_for_check_array_struct_missingness():
    """Fixture to set up a sample Hail Table with array<struct> fields and global index dictionary ."""
    
    # Define rows of the Table.
    data = [
        {"idx": 1, "freq": [{"AC": 10, "AF": 0.1, "AN": 20}, {"AC": 5, "AF": 0.05, "AN": None}]},
        {"idx": 2, "freq": [{"AC": 8, "AF": 0.08, "AN": 16}, {"AC": 6, "AF": None, "AN": None}]},
        {"idx": 3, "freq": [{"AC": 88, "AF": 0.18, "AN": 26}, {"AC": 65, "AF": None, "AN": None}]},
        {"idx": 4, "freq": [{"AC": 8, "AF": 0.08, "AN": 16}, None]},

    ]
    
    # Create Table.
    ht = hl.Table.parallelize(data, hl.tstruct(
        idx=hl.tint32,
        freq=hl.tarray(hl.tstruct(AC=hl.tint32, AF=hl.tfloat64, AN=hl.tint32))
    ))
    
    # Define global annotation for freq_index_dict.
    freq_index_dict = {"adj": 0, "raw": 1}
    ht = ht.annotate_globals(freq_index_dict=freq_index_dict)
    
    return ht


def test_unfurl_array_annotations(ht_for_check_array_struct_missingness):
    """Test the unfurl_array_annotations function for all rows."""
    ht = ht_for_check_array_struct_missingness
    indexed_array_annotations = {"freq": "freq_index_dict"}
    
    # Call the unfurl_array_annotations function.
    result = unfurl_array_annotations(ht, indexed_array_annotations)
    
    # Define names of the expected keys after unfurling.
    expected_keys = {"AC_adj", "AF_adj", "AN_adj", "AC_raw", "AF_raw", "AN_raw"}
    assert set(result.keys()) == expected_keys, "Unfurled keys do not match expected keys."
    
    # Annotate table with unfurled fields.
    ht = ht.annotate(**result)
    rows = ht.collect()

    # Define the expected values for each unfurled annotation.
    expected_values = [
        {"AC_adj": 10, "AF_adj": 0.1, "AN_adj": 20, "AC_raw": 5, "AF_raw": 0.05, "AN_raw": None},
        {"AC_adj": 8, "AF_adj": 0.08, "AN_adj": 16, "AC_raw": 6, "AF_raw": None, "AN_raw": None},
        {"AC_adj": 88, "AF_adj": 0.18, "AN_adj": 26, "AC_raw": 65, "AF_raw": None, "AN_raw": None},
        {"AC_adj": 8, "AF_adj": 0.08, "AN_adj": 16, "AC_raw": None, "AF_raw": None, "AN_raw": None},
    ]

    # Validate each expected value.
    for row, expected in zip(rows, expected_values):
        for key, expected_value in expected.items():
            assert row[key] == expected_value, (
                f"Mismatch in row {row['idx']} for key '{key}': "
                f"expected {expected_value}, got {row[key]}"
        )


def test_check_array_struct_missingness(ht_for_check_array_struct_missingness):
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
        "AN_adj": 0.0,
        "AC_raw": 0.25,
        "AF_raw": 0.75,
        "AN_raw": 1.00 
    }

    # Validate each field's missingness percentage.
    for field, expected_value in expected_missingness.items():
        assert missingness[field] == expected_value, (
            f"Mismatch in missingness for field '{field}': "
            f"expected {expected_value}, got {missingness[field]}"
        )



pytest.main(["-v"])
