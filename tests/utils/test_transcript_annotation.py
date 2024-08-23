"""Tests for the transcript_annotation module."""

import hail as hl
import pytest

from gnomad.utils.transcript_annotation import (
    clean_tissue_name_for_browser,
    create_tx_annotation_by_region,
    tx_filter_variants_by_csqs,
)


class TestCleanTissueNameForBrowser:
    """Tests for the clean_tissue_name_for_browser function."""

    @pytest.mark.parametrize(
        "input_name,expected_output",
        [
            ("basalganglia", "basal_ganglia"),
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
            ("b_a24", "ba24"),
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
                    {
                        "locus": hl.Locus("chr1", 7),
                        "gene_id": "gene2",
                        "exp_prop_mean": 0.8,
                        "tissue1": 0.6,
                        "tissue2": 0.7,
                    },
                    {
                        "locus": hl.Locus("chr1", 8),
                        "gene_id": "gene2",
                        "exp_prop_mean": 0.8,
                        "tissue1": 0.6,
                        "tissue2": 0.7,
                    },
                    {
                        "locus": hl.Locus("chr1", 10),
                        "gene_id": "gene2",
                        "exp_prop_mean": 0.8,
                        "tissue1": 0.6,
                        "tissue2": 0.7,
                    },
                    {
                        "locus": hl.Locus("chr1", 11),
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
                        stop=8,
                        mean=0.8,
                        tissues=hl.Struct(tissue1=0.6, tissue2=0.7),
                    ),
                    hl.Struct(
                        chrom="chr1",
                        start=10,
                        stop=11,
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


@pytest.fixture
def mock_vep_annotated_ht():
    """Create a mock Hail Table with VEP annotations."""
    return hl.Table.parallelize(
        [
            {
                "locus": hl.Locus("1", 861393, reference_genome="GRCh37"),
                "alleles": ["G", "A"],
                "vep": {
                    "transcript_consequences": [
                        {
                            "gene_id": "ENSG00000187634",
                            "gene_symbol": "SAMD11",
                            "transcript_id": "ENST00000342066",
                            "consequence_terms": [
                                "splice_region_variant",
                                "synonymous_variant"
                            ],
                            "amino_acids": "V",
                            "biotype": "protein_coding",
                            "lof": None,
                            "lof_flags": None,
                            "canonical": 0,
                        },
                        {
                            "gene_id": "ENSG00000268179",
                            "gene_symbol": "AL645608.1",
                            "transcript_id": "ENST00000598827",
                            "consequence_terms": ["synonymous_variant"],
                            "amino_acids": "T",
                            "biotype": "protein_coding",
                            "lof": None,
                            "lof_flags": None,
                            "canonical": 1,
                        },
                    ],
                },
            },
            {
                "locus": hl.Locus("1", 871274, reference_genome="GRCh37"),
                "alleles": ["C", "A"],
                "vep": {
                    "transcript_consequences": [
                        {
                            "gene_id": "ENSG00000187634",
                            "gene_symbol": "SAMD11",
                            "transcript_id": "ENST00000420190",
                            "consequence_terms": ["splice_region_variant"],
                            "amino_acids": None,
                            "biotype": "protein_coding",
                            "lof": None,
                            "lof_flags": None,
                            "canonical": 0,
                        }
                    ]
                },
            },
            {
                "locus": hl.Locus("1", 871275, reference_genome="GRCh37"),
                "alleles": ["C", "A"],
                "vep": {
                    "transcript_consequences": [
                        {
                            "gene_id": "ENSG00000187634",
                            "gene_symbol": "SAMD11",
                            "transcript_id": "ENST00000420190",
                            "consequence_terms": [
                                "splice_region_variant",
                                "synonymous_variant"
                            ],
                            "amino_acids": "A",
                            "biotype": "protein_coding",
                            "lof": None,
                            "lof_flags": None,
                            "canonical": 1,
                        }
                    ]
                },
            },
            {
                "locus": hl.Locus("1", 1000, reference_genome="GRCh37"),
                "alleles": ["T", "G"],
                "vep": {
                    "transcript_consequences": [
                        {
                            "gene_id": "ENSG1",
                            "gene_symbol": "gene1",
                            "transcript_id": "ENST1",
                            "consequence_terms": ["stop_gained"],
                            "amino_acids": "Q/*",
                            "biotype": "protein_coding",
                            "lof": None,
                            "lof_flags": None,
                            "canonical": 1,
                        }
                    ]
                },
            },
            {
                "locus": hl.Locus("1", 2000, reference_genome="GRCh37"),
                "alleles": ["A", "T"],
                "vep": {
                    "transcript_consequences": [
                        {
                            "gene_id": "ENSG2",
                            "gene_symbol": "gene2",
                            "transcript_id": "ENST2",
                            "consequence_terms": ["missense_variant"],
                            "amino_acids": "K/R",
                            "biotype": "nonsense_mediated_decay",
                            "lof": None,
                            "lof_flags": None,
                            "canonical": 1,
                        }
                    ]
                },
            },
        ],
        hl.tstruct(
            locus=hl.tlocus(),
            alleles=hl.tarray(hl.tstr),
            vep=hl.tstruct(
                transcript_consequences=hl.tarray(
                    hl.tstruct(
                        gene_id=hl.tstr,
                        gene_symbol=hl.tstr,
                        transcript_id=hl.tstr,
                        consequence_terms=hl.tarray(hl.tstr),
                        amino_acids=hl.tstr,
                        biotype=hl.tstr,
                        lof=hl.tstr,
                        lof_flags=hl.tstr,
                        canonical=hl.tint,
                    )
                )
            ),
        ),
    )


class TestTxFilterVariantsByCsqs:
    """Tests for the tx_filter_variants_by_csqs function."""

    def test_filter_to_cds(self, mock_vep_annotated_ht):
        """Test filtering to CDS variants."""
        result_ht = tx_filter_variants_by_csqs(
            mock_vep_annotated_ht,
            filter_to_cds=True,
            ignore_splicing=False,
            filter_to_protein_coding=False,
        )
        assert result_ht.count() == 2

    def test_filter_to_genes(self, mock_vep_annotated_ht):
        """Test filtering to specific genes."""
        result_ht = tx_filter_variants_by_csqs(
            mock_vep_annotated_ht,
            filter_to_genes=["ENSG1", "ENSG2"],
            filter_to_cds=False,
            ignore_splicing=False,
            filter_to_protein_coding=False,
        )
        assert result_ht.count() == 2

    def test_ignore_splicing(self, mock_vep_annotated_ht):
        """Test ignoring splicing variants."""
        result_ht = tx_filter_variants_by_csqs(
            mock_vep_annotated_ht,
            filter_to_cds=False,
            ignore_splicing=True,
            filter_to_protein_coding=False,
        )
        assert result_ht.count() == 4

    def test_filter_to_protein_coding(self, mock_vep_annotated_ht):
        """Test filtering to protein coding variants."""
        result_ht = tx_filter_variants_by_csqs(
            mock_vep_annotated_ht,
            filter_to_cds=False,
            ignore_splicing=False,
            filter_to_protein_coding=True,
        )
        assert result_ht.count() == 4
