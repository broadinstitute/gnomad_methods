"""Tests for the GRCh38 gnomAD public resource definitions."""

import pytest

from gnomad.resources.grch38 import gnomad
from gnomad.resources.resource_utils import DataException, VersionedTableResource


def _assert_default_version_available(resource: VersionedTableResource) -> None:
    """Assert that a versioned resource's default version is one of its versions."""
    assert resource.default_version in resource.versions


class TestPublicRelease:
    """Tests for public_release."""

    @pytest.mark.parametrize("data_type", ["exomes", "genomes", "joint"])
    def test_default_version_available(self, data_type: str) -> None:
        """Test that the resource builds and its default is a known version."""
        _assert_default_version_available(gnomad.public_release(data_type))


class TestCoverage:
    """Tests for coverage."""

    @pytest.mark.parametrize("data_type", ["exomes", "genomes"])
    def test_default_version_available(self, data_type: str) -> None:
        """Test that the resource builds and its default is a known version."""
        _assert_default_version_available(gnomad.coverage(data_type))


class TestAllSitesAn:
    """Tests for all_sites_an."""

    @pytest.mark.parametrize("data_type", ["exomes", "genomes"])
    def test_default_version_available(self, data_type: str) -> None:
        """Test that the resource builds and its default is a known version."""
        _assert_default_version_available(gnomad.all_sites_an(data_type))


class TestConstraint:
    """Tests for constraint."""

    def test_default_version_available(self) -> None:
        """Test that the resource builds and its default is a known version."""
        _assert_default_version_available(gnomad.constraint())


class TestConstraintMutationRate:
    """Tests for constraint_mutation_rate."""

    def test_default_version_available(self) -> None:
        """
        Test that the resource builds and its default is a known version.

        Guards against the default drifting to a constraint release that has no
        matching mutation rate table.
        """
        _assert_default_version_available(gnomad.constraint_mutation_rate())


class TestBrowserVariant:
    """Tests for browser_variant."""

    def test_default_version_available(self) -> None:
        """Test that the resource builds and its default is a known version."""
        _assert_default_version_available(gnomad.browser_variant())


class TestSvSitesVcfPath:
    """Tests for sv_sites_vcf_path."""

    def test_default_version(self) -> None:
        """Test that the default version gives the released v4.1 SV sites VCF."""
        assert gnomad.sv_sites_vcf_path() == (
            "gs://gcp-public-data--gnomad/release/4.1/genome_sv/"
            "gnomad.v4.1.sv.sites.vcf.gz"
        )

    def test_unreleased_version_raises(self) -> None:
        """
        Test that v4.0 is rejected.

        v4.0 SVs were released only as per-chromosome VCFs, so the sites path
        pattern doesn't exist for it.
        """
        with pytest.raises(DataException):
            gnomad.sv_sites_vcf_path("4.0")
