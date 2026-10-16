"""Tests for the GRCh38 gnomAD public resource definitions."""

from functools import partial

import pytest

from gnomad.resources.grch38 import gnomad

VERSIONED_RESOURCES = {
    "public_release_exomes": partial(gnomad.public_release, "exomes"),
    "public_release_genomes": partial(gnomad.public_release, "genomes"),
    "public_release_joint": partial(gnomad.public_release, "joint"),
    "coverage_exomes": partial(gnomad.coverage, "exomes"),
    "coverage_genomes": partial(gnomad.coverage, "genomes"),
    "all_sites_an_exomes": partial(gnomad.all_sites_an, "exomes"),
    "all_sites_an_genomes": partial(gnomad.all_sites_an, "genomes"),
    "constraint": gnomad.constraint,
    "constraint_mutation_rate": gnomad.constraint_mutation_rate,
    "browser_variant": gnomad.browser_variant,
}


@pytest.mark.parametrize(
    "get_resource", VERSIONED_RESOURCES.values(), ids=VERSIONED_RESOURCES.keys()
)
def test_versioned_resource_default_version_is_available(get_resource):
    """Test that each versioned resource can be built and its default is a known version."""
    resource = get_resource()
    assert resource.default_version in resource.versions
