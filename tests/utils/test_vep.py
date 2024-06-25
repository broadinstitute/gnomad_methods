"""Tests for the gnomAD VEP utility functions."""
import hail as hl
from gnomad.utils.vep import (
    add_most_severe_consequence_to_consequence,
    get_most_severe_consequence_expr,
    add_most_severe_csq_to_tc_within_vep_root,
    filter_to_most_severe_consequences,
    CSQ_ORDER,
    filter_vep_consequences_by_loftee,
    get_most_severe_csq_from_multiple_csq_lists,
    process_consequences,
)


def test_get_most_severe_consequence_expr():
    """Test the get_most_severe_consequence_expr function."""
    # Create a mock ArrayExpression of consequences
    csq_expr = hl.literal(['splice_region_variant', 'intron_variant'])

    # Test with default csq_order
    result = get_most_severe_consequence_expr(csq_expr)
    assert result.value == 'splice_region_variant', "Expected 'splice_region_variant'"

    # Test with custom csq_order
    custom_csq_order = ['intron_variant', 'splice_region_variant']
    result = get_most_severe_consequence_expr(csq_expr, custom_csq_order)
    assert result.value == 'intron_variant', "Expected 'intron_variant'"

    # Test with csq_expr that contains a consequence not in csq_order
    csq_expr = hl.literal(['non_existent_consequence'])
    result = get_most_severe_consequence_expr(csq_expr)
    assert result.value is None, "Expected None for non-existent consequence"


def test_add_most_severe_consequence_to_consequence():
    """Test the add_most_severe_consequence_to_consequence function."""
    # Create a mock StructExpression of transcript consequence
    tc_expr = hl.struct(consequence_terms=['missense_variant', 'synonymous_variant'])

    # Test with default csq_order
    result = add_most_severe_consequence_to_consequence(tc_expr)
    assert result.most_severe_consequence == get_most_severe_consequence_expr(tc_expr.consequence_terms), "Expected 'missense_variant'"

    # Test with custom csq_order
    custom_csq_order = ['synonymous_variant', 'missense_variant']
    result = add_most_severe_consequence_to_consequence(tc_expr, custom_csq_order)
    assert result.most_severe_consequence == get_most_severe_consequence_expr(tc_expr.consequence_terms, custom_csq_order), "Expected 'synonymous_variant'"

    # Test with csq_expr that contains a consequence not in csq_order
    tc_expr = hl.struct(consequence_terms=['non_existent_consequence'])
    result = add_most_severe_consequence_to_consequence(tc_expr)
    assert result.most_severe_consequence == get_most_severe_consequence_expr(tc_expr.consequence_terms), "Expected None for non-existent consequence"


def test_add_most_severe_csq_to_tc_within_vep_root():
    """Test the add_most_severe_csq_to_tc_within_vep_root function."""
    # Create a mock MatrixTable with vep.transcript_consequences
    mt = hl.utils.range_matrix_table(1, 1)
    mt = mt.annotate_rows(vep=hl.struct(transcript_consequences=[hl.struct(consequence_terms=['missense_variant', 'synonymous_variant'])]))

    # Apply the function
    result_mt = add_most_severe_csq_to_tc_within_vep_root(mt)

    # Check that the most_severe_consequence field is added
    assert 'most_severe_consequence' in result_mt.vep.transcript_consequences[0]

    # Check that the most_severe_consequence is correct
    assert result_mt.vep.transcript_consequences[0].most_severe_consequence == 'missense_variant'


def test_filter_to_most_severe_consequences():
    """Test the filter_to_most_severe_consequences function."""
    # Create a mock ArrayExpression of VEP consequences
    csq_expr = hl.literal([
        hl.struct(most_severe_consequence='missense_variant'),
        hl.struct(most_severe_consequence='synonymous_variant'),
        hl.struct(most_severe_consequence='missense_variant')
    ])

    # Apply the function
    result = filter_to_most_severe_consequences(
        csq_expr, CSQ_ORDER, 'most_severe_consequence'
    )

    # Check that the result only contains the most severe consequences
    assert len(result) == 2, "Expected 2 consequences"
    assert all(csq.most_severe_consequence == 'missense_variant' for csq in
               result), "Expected 'missense_variant'"


def test_filter_vep_consequences_by_loftee():
    """Test the filter_vep_consequences_by_loftee function."""
    # Create a mock ArrayExpression of VEP consequences with LOFTEE annotations
    csq_expr = hl.literal([
        hl.struct(lof='HC', lof_flags=None),
        hl.struct(lof='HC', lof_flags=''),
        hl.struct(lof='LC', lof_flags='flag1'),
        hl.struct(lof='OS', lof_flags=None),
        hl.struct(lof=None, lof_flags='flag2')
    ])

    # Test with default parameters
    result = filter_vep_consequences_by_loftee(csq_expr)
    assert len(result) == 6, "Expected 6 consequences"

    # Test with loftee_labels
    result = filter_vep_consequences_by_loftee(csq_expr, loftee_labels=['HC'])
    assert len(result) == 2, "Expected 2 consequences"
    assert result[0].lof == 'HC', "Expected 'HC'"

    # Test with no_lof_flags
    result = filter_vep_consequences_by_loftee(csq_expr, no_lof_flags=True)
    assert len(result) == 3, "Expected 3 consequences"
    assert all(csq.lof_flags is None for csq in result), "Expected no LOFTEE flags"

    # Test with loftee_labels and no_lof_flags
    result = filter_vep_consequences_by_loftee(csq_expr, loftee_labels=['HC'],
                                               no_lof_flags=True)
    assert len(result) == 1, "Expected 1 consequence"
    assert result[0].lof == 'HC', "Expected 'HC'"
    assert all(csq.lof_flags is None for csq in result), "Expected no LOFTEE flags"

    # Test with keep=False
    result = filter_vep_consequences_by_loftee(csq_expr, loftee_labels=['HC'],
                                               keep=False)
    assert len(result) == 3, "Expected 3 consequences"
    assert all(csq.lof != 'HC' for csq in result), "Expected no 'HC' consequences"


def test_get_most_severe_csq_from_multiple_csq_lists():
    """Test the get_most_severe_csq_from_multiple_csq_lists function."""
    # Create a mock VEP expression
    vep_expr = hl.struct(
        transcript_consequences=[
            hl.struct(
                biotype='protein_coding',
                lof='HC',
                most_severe_consequence='missense_variant'
            ),
            hl.struct(
                biotype='protein_coding',
                lof='HC',
                most_severe_consequence='synonymous_variant'
            )
        ],
        intergenic_consequences=[
            hl.struct(
                biotype='intergenic',
                lof=None,
                most_severe_consequence='intergenic_variant'
            )
        ]
    )

    # Test the function
    result = get_most_severe_csq_from_multiple_csq_lists(vep_expr)

    # Check the most_severe_csq
    assert result.most_severe_csq == 'missense_variant', "Expected 'missense_variant'"

    # Check the protein_coding
    assert result.protein_coding == True, "Expected True"

    # Check the lof
    assert result.lof == 'HC', "Expected 'HC'"

    # Check the no_lof_flags
    assert result.no_lof_flags == False, "Expected False"

    # Check the transcript_consequences
    assert len(result.transcript_consequences) == 2, "Expected 2 transcript consequences"
    assert result.transcript_consequences[0].most_severe_consequence == 'missense_variant', "Expected 'missense_variant' for the first transcript consequence"
    assert result.transcript_consequences[1].most_severe_consequence == 'synonymous_variant', "Expected 'synonymous_variant' for the second transcript consequence"
