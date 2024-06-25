"""Tests for the gnomAD VEP utility functions."""
import hail as hl

from gnomad.utils.vep import (
    add_most_severe_consequence_to_consequence,
    add_most_severe_csq_to_tc_within_vep_root,
    filter_to_most_severe_consequences,
    filter_vep_transcript_csqs_expr,
    get_most_severe_consequence_expr,
    get_most_severe_csq_from_multiple_csq_lists,
)


def test_get_most_severe_consequence_expr():
    """Test the get_most_severe_consequence_expr function."""
    # Create a mock ArrayExpression of consequences
    csq_expr = hl.literal(["splice_region_variant", "intron_variant"])

    # Test with default csq_order
    result = get_most_severe_consequence_expr(csq_expr)
    assert (
        hl.eval(result) == "splice_region_variant"
    ), "Expected 'splice_region_variant'"

    # Test with custom csq_order
    custom_csq_order = ["intron_variant", "splice_region_variant"]
    result = get_most_severe_consequence_expr(csq_expr, custom_csq_order)
    assert hl.eval(result) == "intron_variant", "Expected 'intron_variant'"

    # Test with csq_expr that contains a consequence not in csq_order
    csq_expr = hl.literal(["non_existent_consequence"])
    result = get_most_severe_consequence_expr(csq_expr)
    assert hl.eval(result) is None, "Expected None for non-existent consequence"


def test_add_most_severe_consequence_to_consequence():
    """Test the add_most_severe_consequence_to_consequence function."""
    # Create a mock StructExpression of transcript consequence
    tc_expr = hl.struct(consequence_terms=["missense_variant", "synonymous_variant"])

    # Test with default csq_order
    result = add_most_severe_consequence_to_consequence(tc_expr)
    assert hl.eval(result.most_severe_consequence) == hl.eval(
        get_most_severe_consequence_expr(tc_expr.consequence_terms)
    ), "Expected 'missense_variant'"

    # Test with custom csq_order
    custom_csq_order = ["synonymous_variant", "missense_variant"]
    result = add_most_severe_consequence_to_consequence(tc_expr, custom_csq_order)
    assert hl.eval(result.most_severe_consequence) == hl.eval(
        get_most_severe_consequence_expr(tc_expr.consequence_terms, custom_csq_order)
    ), "Expected 'synonymous_variant'"

    # Test with csq_expr that contains a consequence not in csq_order
    tc_expr = hl.struct(consequence_terms=["non_existent_consequence"])
    result = add_most_severe_consequence_to_consequence(tc_expr)
    assert hl.eval(result.most_severe_consequence) == hl.eval(
        get_most_severe_consequence_expr(tc_expr.consequence_terms)
    ), "Expected None for non-existent consequence"


def test_add_most_severe_csq_to_tc_within_vep_root():
    """Test the add_most_severe_csq_to_tc_within_vep_root function."""
    # Create a mock MatrixTable with vep.transcript_consequences
    mt = hl.utils.range_matrix_table(1, 1)
    mt = mt.annotate_rows(
        vep=hl.struct(
            transcript_consequences=[
                hl.struct(consequence_terms=["missense_variant", "synonymous_variant"])
            ]
        )
    )

    # Apply the function
    result_mt = add_most_severe_csq_to_tc_within_vep_root(mt)

    # Check that the most_severe_consequence field is added
    assert (
        "most_severe_consequence" in result_mt.vep.transcript_consequences[0].keys()
    ), "Expected 'most_severe_consequence' field"

    # Check that the most_severe_consequence is correct
    assert (
        result_mt.vep.transcript_consequences[0].most_severe_consequence.take(1)[0]
        == "missense_variant"
    )


def test_filter_to_most_severe_consequences():
    """Test the filter_to_most_severe_consequences function."""
    # Create a mock array of VEP consequences
    mock_csq_expr = hl.literal(
        [
            hl.struct(
                consequence_terms=["missense_variant", "synonymous_variant"],
                lof="HC",
                lof_flags="HC",
                protein_coding=1,
            ),
            hl.struct(
                consequence_terms=["splice_region_variant", "synonymous_variant"],
                lof="LC",
                lof_flags="HC",
                protein_coding=1,
            ),
            hl.struct(
                consequence_terms=["frameshift_variant", "synonymous_variant"],
                lof="HC",
                lof_flags="",
                protein_coding=1,
            ),
        ]
    )

    # Call the function with the mock data
    result = filter_to_most_severe_consequences(mock_csq_expr, prioritize_loftee=True)

    # Assert that the output is as expected
    assert (
        hl.eval(result.most_severe_consequence) == "frameshift_variant"
    ), "Expected 'frameshift_variant'"
    assert hl.eval(result.lof) == "HC", "Expected 'HC'"
    assert hl.eval(result.no_lof_flags)
    assert hl.eval(hl.len(result.consequences)) == 1, "Expected 1 consequence"
    assert hl.eval(result.consequences[0].consequence_terms) == [
        "frameshift_variant",
        "synonymous_variant",
    ], "Expected ['frameshift_variant', 'synonymous_variant']"


def test_filter_vep_transcript_csqs_expr_loftee():
    """Test the filter_vep_consequences_by_loftee function."""
    # Create a mock ArrayExpression of VEP consequences with LOFTEE annotations
    csq_expr = hl.literal(
        [
            hl.struct(lof="HC", lof_flags=hl.missing(hl.tstr)),
            hl.struct(lof="HC", lof_flags=""),
            hl.struct(lof="LC", lof_flags="flag1"),
            hl.struct(lof="OS", lof_flags=hl.missing(hl.tstr)),
            hl.struct(lof=hl.missing(hl.tstr), lof_flags="flag2"),
        ]
    )

    # Test with default parameters
    result = filter_vep_transcript_csqs_expr(csq_expr)
    assert hl.eval(hl.len(result)) == 5, "Expected 5 consequences"

    # Test with loftee_labels
    result = filter_vep_transcript_csqs_expr(csq_expr, loftee_labels=["HC"])
    assert hl.eval(hl.len(result)) == 2, "Expected 2 consequences"
    assert hl.eval(result[0].lof) == "HC", "Expected 'HC'"

    # Test with no_lof_flags
    result = filter_vep_transcript_csqs_expr(csq_expr, no_lof_flags=True)
    assert hl.eval(hl.len(result)) == 3, "Expected 3 consequences"
    assert hl.eval(
        hl.all(
            result.map(lambda csq: hl.is_missing(csq.lof_flags) | (csq.lof_flags == ""))
        )
    ), "Expected no LOFTEE flags"

    # Test with loftee_labels and no_lof_flags
    result = filter_vep_transcript_csqs_expr(
        csq_expr, loftee_labels=["HC"], no_lof_flags=True
    )
    assert hl.eval(hl.len(result)) == 2, "Expected 2 consequence"
    assert hl.eval(result[0].lof) == "HC", "Expected 'HC'"
    assert hl.eval(
        hl.all(
            result.map(lambda csq: hl.is_missing(csq.lof_flags) | (csq.lof_flags == ""))
        )
    ), "Expected no LOFTEE flags"


def test_get_most_severe_csq_from_multiple_csq_lists():
    """Test the get_most_severe_csq_from_multiple_csq_lists function."""
    # Create a mock VEP expression
    vep_expr = hl.struct(
        transcript_consequences=[
            hl.struct(
                biotype="protein_coding",
                lof="HC",
                lof_flags="flag1",
                consequence_terms="missense_variant",
            ),
            hl.struct(
                biotype="protein_coding",
                lof="HC",
                lof_flags=hl.missing(hl.tstr),
                consequence_terms="synonymous_variant",
            ),
            hl.struct(
                biotype="protein_coding",
                lof="HC",
                lof_flags=hl.missing(hl.tstr),
                consequence_terms="missense_variant",
            ),
        ],
        intergenic_consequences=[
            hl.struct(
                biotype="intergenic",
                lof=hl.missing(hl.tstr),
                lof_flags=hl.missing(hl.tstr),
                consequence_terms="intergenic_variant",
            )
        ],
    )

    # Test the function
    result = get_most_severe_csq_from_multiple_csq_lists(vep_expr)

    # Check the most_severe_consequence
    assert (
        hl.eval(result.most_severe_consequence) == "missense_variant"
    ), "Expected 'missense_variant'"

    # Check the protein_coding
    assert hl.eval(result.protein_coding), "Expected True"

    # Check the lof
    assert hl.eval(result.lof) == "HC", "Expected 'HC'"

    # Check the no_lof_flags
    assert hl.eval(result.no_lof_flags), "Expected True"

    # Check the transcript_consequences
    result = get_most_severe_csq_from_multiple_csq_lists(vep_expr, include_csqs=True)
    assert (
        hl.eval(hl.len(result.transcript_consequences)) == 2
    ), "Expected 2 transcript consequences"

    # Create a mock VEP expression
    vep_expr = hl.struct(
        transcript_consequences=[
            hl.struct(
                biotype="protein_coding",
                lof="HC",
                lof_flags="flag1",
                consequence_terms="missense_variant",
            )
        ]
    )

    # Test the function
    result = get_most_severe_csq_from_multiple_csq_lists(vep_expr)

    # Check the no_lof_flags
    assert hl.eval(result.no_lof_flags) == False, "Expected False"
