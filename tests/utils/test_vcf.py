"""Tests for the vcf utility module."""

import hail as hl

from gnomad.utils.vcf import build_vcf_export_reference

GRCH38_PRIMARY = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


class TestBuildVcfExportReference:
    """Test the build_vcf_export_reference function."""

    def test_default_contigs(self) -> None:
        """Test that the default keeps chr1-22, X, Y, and M for GRCh38."""
        ref = build_vcf_export_reference("test_default_contigs")
        assert ref.contigs == GRCH38_PRIMARY + ["chrM"]
        assert ref.mt_contigs == ["chrM"]
        assert ref.x_contigs == ["chrX"]
        assert ref.y_contigs == ["chrY"]

    def test_default_contigs_no_chrm(self) -> None:
        """Test that keep_chrM=False drops chrM."""
        ref = build_vcf_export_reference("test_no_chrm", keep_chrM=False)
        assert ref.contigs == GRCH38_PRIMARY

    def test_default_follows_build(self) -> None:
        """Test that the default keep_contigs follows the build, not GRCh38 naming."""
        ref = build_vcf_export_reference("test_grch37", build="GRCh37")
        assert ref.contigs == [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]

    def test_repeated_calls_do_not_accumulate_chrm(self) -> None:
        """Test that calling twice does not mutate a shared default and duplicate chrM."""
        build_vcf_export_reference("test_repeat_1")
        ref = build_vcf_export_reference("test_repeat_2")
        assert ref.contigs.count("chrM") == 1

    def test_explicit_keep_contigs_not_mutated(self) -> None:
        """Test that a caller-supplied keep_contigs list is not mutated."""
        keep = ["chr1", "chrX", "chrY"]
        ref = build_vcf_export_reference("test_explicit", keep_contigs=keep)
        assert keep == ["chr1", "chrX", "chrY"]
        assert ref.contigs == ["chr1", "chrX", "chrY", "chrM"]
        assert ref.lengths["chr1"] == hl.get_reference("GRCh38").lengths["chr1"]

    def test_keep_contigs_with_chrm_not_duplicated(self) -> None:
        """Test that chrM already in keep_contigs is not added again."""
        ref = build_vcf_export_reference(
            "test_chrm_dedup", keep_contigs=["chr1", "chrM", "chrX", "chrY"]
        )
        assert ref.contigs == ["chr1", "chrM", "chrX", "chrY"]
        assert ref.mt_contigs == ["chrM"]
