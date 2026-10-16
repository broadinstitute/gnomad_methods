"""Tests for the GRCh38 gnomAD public resource definitions."""

import pytest

from gnomad.resources.grch38 import gnomad
from gnomad.resources.resource_utils import VersionedTableResource


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
