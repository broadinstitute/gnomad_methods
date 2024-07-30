"""Tests for the transcript_annotation module."""

import hail as hl
import pytest

from gnomad.utils.transcript_annotation import (
    clean_tissue_name_for_browser,
    create_tx_annotation_by_region,
)


class TestCleanTissueNameForBrowser:
    """Tests for the clean_tissue_name_for_browser function."""

    @pytest.mark.parametrize(
        "input_name,expected_output",
        [
            ("BasalGanglia", "basal_ganglia"),
            ("NucleusAccumbens", "nucleus_accumbens"),
            ("SpinalCord", "spinal_cord"),
            ("CervicalC", "cervical_c"),
            ("SubstantiaNigra", "substantia_nigra"),
            ("CulturedFibroblasts", "cultured_fibroblasts"),
            ("LowerLeg", "lower_leg"),
            ("TransformedLymphocytes", "transformed_lymphocytes"),
            ("AnteriorCingulateCortex", "anterior_cingulate_cortex"),
            ("BA24", "ba24"),
            ("BA9", "ba9"),
            ("EBV", "ebv"),
            ("Heart", "heart"),
            ("Liver", "liver"),
            ("braincortex", "braincortex"),
        ],
    )
    def test_clean_tissue_name_for_browser(
        self, input_name: str, expected_output: str
    ) -> None:
        """Test the clean_tissue_name_for_browser function."""
        result = clean_tissue_name_for_browser(input_name)
        assert result == expected_output


class TestCreateTxAnnotationByRegion:
    """Tests for the create_tx_annotation_by_region function."""

    @pytest.fixture
    def sample_hail_table(self) -> hl.Table:
        """Fixture to create a sample Hail Table."""
        return (
            hl.Table.parallelize(
                [
                    {
                        "locus": hl.Locus("chr1", 1),
                        "gene_id": "gene1",
                        "exp_prop_mean": 0.5,
                        "tissue1": 0.2,
                        "tissue2": 0.3,
                    },
                    {
                        "locus": hl.Locus("chr1", 2),
                        "gene_id": "gene1",
                        "exp_prop_mean": 0.5,
                        "tissue1": 0.2,
                        "tissue2": 0.3,
                    },
                    {
                        "locus": hl.Locus("chr1", 3),
                        "gene_id": "gene1",
                        "exp_prop_mean": 0.6,
                        "tissue1": 0.3,
                        "tissue2": 0.4,
                    },
                    {
                        "locus": hl.Locus("chr1", 4),
                        "gene_id": "gene2",
                        "exp_prop_mean": 0.7,
                        "tissue1": 0.5,
                        "tissue2": 0.6,
                    },
                    {
                        "locus": hl.Locus("chr1", 5),
                        "gene_id": "gene2",
                        "exp_prop_mean": 0.7,
                        "tissue1": 0.5,
                        "tissue2": 0.6,
                    },
                    {
                        "locus": hl.Locus("chr1", 6),
                        "gene_id": "gene2",
                        "exp_prop_mean": 0.8,
                        "tissue1": 0.6,
                        "tissue2": 0.7,
                    },
                ],
                hl.tstruct(
                    locus=hl.tlocus(),
                    gene_id=hl.tstr,
                    exp_prop_mean=hl.tfloat64,
                    tissue1=hl.tfloat64,
                    tissue2=hl.tfloat64,
                ),
            )
            .key_by("locus")
            .annotate_globals(tissues=["tissue1", "tissue2"])
        )

    def test_create_tx_annotation_by_region(self, sample_hail_table: hl.Table) -> None:
        """Test the create_tx_annotation_by_region function."""
        result_ht = create_tx_annotation_by_region(sample_hail_table)

        # Expected result
        expected_result = [
            hl.Struct(
                gene_id="gene1",
                regions=[
                    hl.Struct(
                        chrom="chr1",
                        start=1,
                        stop=2,
                        mean=0.5,
                        tissues=hl.Struct(tissue1=0.2, tissue2=0.3),
                    ),
                    hl.Struct(
                        chrom="chr1",
                        start=3,
                        stop=3,
                        mean=0.6,
                        tissues=hl.Struct(tissue1=0.3, tissue2=0.4),
                    ),
                ],
            ),
            hl.Struct(
                gene_id="gene2",
                regions=[
                    hl.Struct(
                        chrom="chr1",
                        start=4,
                        stop=5,
                        mean=0.7,
                        tissues=hl.Struct(tissue1=0.5, tissue2=0.6),
                    ),
                    hl.Struct(
                        chrom="chr1",
                        start=6,
                        stop=6,
                        mean=0.8,
                        tissues=hl.Struct(tissue1=0.6, tissue2=0.7),
                    ),
                ],
            ),
        ]

        # Collect results
        result = result_ht.collect()

        # Verify the result
        assert result == expected_result
