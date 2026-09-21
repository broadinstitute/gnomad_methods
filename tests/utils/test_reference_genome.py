"""Tests for the reference_genome utility module."""

from typing import List

import pytest

from gnomad.utils.reference_genome import get_primary_contigs

GRCH37_AUTOSOMES = [str(i) for i in range(1, 23)]
GRCH38_AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


class TestGetPrimaryContigs:
    """Test the get_primary_contigs function."""

    @pytest.mark.parametrize(
        "build, include_sex, keep_chrM, expected",
        [
            ("GRCh38", True, False, GRCH38_AUTOSOMES + ["chrX", "chrY"]),
            ("GRCh38", False, False, GRCH38_AUTOSOMES),
            ("GRCh38", True, True, GRCH38_AUTOSOMES + ["chrX", "chrY", "chrM"]),
            ("GRCh38", False, True, GRCH38_AUTOSOMES + ["chrM"]),
            ("GRCh37", True, False, GRCH37_AUTOSOMES + ["X", "Y"]),
            ("GRCh37", False, False, GRCH37_AUTOSOMES),
            ("GRCh37", True, True, GRCH37_AUTOSOMES + ["X", "Y", "MT"]),
        ],
    )
    def test_get_primary_contigs(
        self, build: str, include_sex: bool, keep_chrM: bool, expected: List[str]
    ) -> None:
        """Test that the expected contigs are returned in reference order."""
        assert get_primary_contigs(build, include_sex, keep_chrM) == expected

    def test_default_build_is_grch38(self) -> None:
        """Test that the default build is GRCh38 with sex contigs and no MT."""
        assert get_primary_contigs() == GRCH38_AUTOSOMES + ["chrX", "chrY"]

    def test_no_alt_contigs(self) -> None:
        """Test that alt/decoy contigs are never included."""
        assert all(
            "_" not in c for c in get_primary_contigs("GRCh38", include_chrM=True)
        )

    def test_unsupported_build_raises(self) -> None:
        """Test that non-human builds raise instead of returning positional contigs."""
        with pytest.raises(NotImplementedError, match="GRCm38"):
            get_primary_contigs("GRCm38")
