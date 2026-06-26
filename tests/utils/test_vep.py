"""Tests for the VEP utility module."""

import io
import json
from unittest.mock import MagicMock, patch

import hail as hl
import pytest

from gnomad.utils.vep import (
    get_loftee_end_trunc_filter_expr,
    get_vep_help,
    update_loftee_end_trunc_filter,
    vep_or_lookup_vep,
)


def _mock_open_returning(content):
    """Build a mock for `hfs.open` whose context manager yields `content`.

    `content` is wrapped in a `StringIO` so both `json.load(f)` and `f.read()`
    work against it.
    """
    ctx = MagicMock()
    ctx.__enter__.return_value = io.StringIO(content)
    ctx.__exit__.return_value = False
    return MagicMock(return_value=ctx)


def _fake_vep_context(default_version, versions):
    """Build a fake VEP context resource standing in for `get_vep_context`.

    `versions` maps a version string to the context Hail Table that its
    `.ht()` should return; `default_version` is exposed as `.default_version`.
    """
    context = MagicMock()
    context.default_version = default_version
    context.versions = {
        version: MagicMock(ht=MagicMock(return_value=ht))
        for version, ht in versions.items()
    }
    return context


class TestGetLofteeEndTruncFilterExpr:
    """Test the get_loftee_end_trunc_filter_expr function."""

    @pytest.fixture
    def sample_csq_structs(self):
        """Fixture to create sample consequence structs with different LOFTEE annotations."""
        return [
            # Case 1: GERP_DIST < 0, 50_BP_RULE != PASS -> should be True for default.
            hl.Struct(lof_info="GERP_DIST:-2.5,50_BP_RULE:FAIL,OTHER:value"),
            # Case 2: GERP_DIST >= 0, 50_BP_RULE != PASS -> should be False for default.
            hl.Struct(lof_info="GERP_DIST:1.5,50_BP_RULE:FAIL,OTHER:value"),
            # Case 3: GERP_DIST < 0, 50_BP_RULE = PASS -> should be False for default.
            hl.Struct(lof_info="GERP_DIST:-1.0,50_BP_RULE:PASS,OTHER:value"),
            # Case 4: GERP_DIST >= 0, 50_BP_RULE = PASS -> should be False for default.
            hl.Struct(lof_info="GERP_DIST:0.5,50_BP_RULE:PASS,OTHER:value"),
            # Case 5: GERP_DIST >= 0, 50_BP_RULE != PASS -> should be False for default.
            hl.Struct(lof_info="GERP_DIST:0.5,50_BP_RULE:FAIL,OTHER:value"),
            # Case 6: Missing GERP_DIST (defaults to 0), 50_BP_RULE != PASS -> should
            # be False for default.
            hl.Struct(lof_info="50_BP_RULE:FAIL,OTHER:value"),
            # Case 7: GERP_DIST < 0, missing 50_BP_RULE (defaults to empty) -> should
            # be True for default.
            hl.Struct(lof_info="GERP_DIST:-1.5,OTHER:value"),
            # Case 8: Empty lof_info -> should be False for default.
            hl.Struct(lof_info=""),
        ]

    def test_default_cutoff(self, sample_csq_structs):
        """Test the function with default cutoff of 0.0."""
        ht = hl.Table.parallelize(
            [{"csq": csq} for csq in sample_csq_structs],
            hl.tstruct(csq=hl.tstruct(lof_info=hl.tstr)),
        )

        # Apply the function.
        ht = ht.annotate(end_trunc=get_loftee_end_trunc_filter_expr(ht.csq))

        # Collect results
        results = ht.collect()

        # Expected results for default cutoff (0.0).
        expected = [True, False, False, False, False, False, True, False]

        assert [r.end_trunc for r in results] == expected

    def test_custom_cutoff_positive(self, sample_csq_structs):
        """Test the function with a positive cutoff."""
        ht = hl.Table.parallelize(
            [{"csq": csq} for csq in sample_csq_structs],
            hl.tstruct(csq=hl.tstruct(lof_info=hl.tstr)),
        )

        # Apply the function with cutoff 1.0.
        ht = ht.annotate(
            end_trunc=get_loftee_end_trunc_filter_expr(ht.csq, gerp_dist_cutoff=1.0)
        )

        # Collect results
        results = ht.collect()

        # Expected results for cutoff 1.0.
        expected = [True, False, False, False, True, True, True, True]

        assert [r.end_trunc for r in results] == expected

    def test_custom_cutoff_negative(self, sample_csq_structs):
        """Test the function with a negative cutoff."""
        ht = hl.Table.parallelize(
            [{"csq": csq} for csq in sample_csq_structs],
            hl.tstruct(csq=hl.tstruct(lof_info=hl.tstr)),
        )

        # Apply the function with cutoff -1.0.
        ht = ht.annotate(
            end_trunc=get_loftee_end_trunc_filter_expr(ht.csq, gerp_dist_cutoff=-1.0)
        )

        # Collect results
        results = ht.collect()

        # Expected results for cutoff -1.0.
        expected = [True, False, False, False, False, False, True, False]

        assert [r.end_trunc for r in results] == expected


class TestUpdateLofteeEndTruncFilter:
    """Test the update_loftee_end_trunc_filter function."""

    @pytest.fixture
    def sample_csq_structs_with_filters(self):
        """Fixture to create sample consequence structs with lof_filter and lof annotations."""
        return [
            # Case 1: Should add END_TRUNC filter.
            hl.Struct(
                lof_info="GERP_DIST:-2.5,50_BP_RULE:FAIL",
                lof_filter="SINGLE_EXON",
                lof="HC",
            ),
            # Case 2: Should not add END_TRUNC filter.
            hl.Struct(
                lof_info="GERP_DIST:1.5,50_BP_RULE:PASS",
                lof_filter="SINGLE_EXON",
                lof="HC",
            ),
            # Case 3: Should remove existing END_TRUNC filter.
            hl.Struct(
                lof_info="GERP_DIST:1.0,50_BP_RULE:PASS",
                lof_filter="SINGLE_EXON,END_TRUNC",
                lof="LC",
            ),
            # Case 4: Should add END_TRUNC.
            hl.Struct(
                lof_info="GERP_DIST:-1.0,50_BP_RULE:FAIL", lof_filter="", lof="HC"
            ),
            # Case 5: Missing lof_filter.
            hl.Struct(
                lof_info="GERP_DIST:-1.5,50_BP_RULE:FAIL", lof_filter=None, lof="HC"
            ),
        ]

    def test_update_single_struct(self, sample_csq_structs_with_filters):
        """Test updating a single consequence struct."""
        ht = hl.Table.parallelize(
            [{"csq": csq} for csq in sample_csq_structs_with_filters],
            hl.tstruct(
                csq=hl.tstruct(lof_info=hl.tstr, lof_filter=hl.tstr, lof=hl.tstr)
            ),
        )

        # Apply the function.
        ht = ht.annotate(updated_csq=update_loftee_end_trunc_filter(ht.csq))

        # Collect results.
        results = ht.collect()

        # Check results.
        assert results[0].updated_csq.lof_filter == "END_TRUNC,SINGLE_EXON"
        assert results[0].updated_csq.lof == "LC"

        # Still LC because filter is not empty.
        assert results[1].updated_csq.lof_filter == "SINGLE_EXON"
        assert results[1].updated_csq.lof == "LC"

        # Still LC because filter is not empty.
        assert results[2].updated_csq.lof_filter == "SINGLE_EXON"
        assert results[2].updated_csq.lof == "LC"

        assert results[3].updated_csq.lof_filter == "END_TRUNC"
        assert results[3].updated_csq.lof == "LC"

        assert results[4].updated_csq.lof_filter == "END_TRUNC"
        assert results[4].updated_csq.lof == "LC"

    def test_update_array_of_structs(self, sample_csq_structs_with_filters):
        """Test updating an array of consequence structs."""
        # Create a table with arrays of consequences.
        ht = hl.Table.parallelize(
            [
                {"csqs": sample_csq_structs_with_filters[:2]},
                {"csqs": sample_csq_structs_with_filters[2:]},
            ],
            hl.tstruct(
                csqs=hl.tarray(
                    hl.tstruct(lof_info=hl.tstr, lof_filter=hl.tstr, lof=hl.tstr)
                )
            ),
        )

        # Apply the function.
        ht = ht.annotate(updated_csqs=update_loftee_end_trunc_filter(ht.csqs))

        # Collect results.
        results = ht.collect()

        # Check first array.
        first_array = results[0].updated_csqs
        assert first_array[0].lof_filter == "END_TRUNC,SINGLE_EXON"
        assert first_array[0].lof == "LC"
        assert first_array[1].lof_filter == "SINGLE_EXON"
        assert first_array[1].lof == "LC"

        # Check second array.
        second_array = results[1].updated_csqs
        assert second_array[0].lof_filter == "SINGLE_EXON"
        assert second_array[0].lof == "LC"
        assert second_array[1].lof_filter == "END_TRUNC"
        assert second_array[1].lof == "LC"
        assert second_array[2].lof_filter == "END_TRUNC"
        assert second_array[2].lof == "LC"

    def test_missing_lof_annotation(self):
        """Test updating when lof annotation is missing."""
        csq_with_missing_lof = hl.Struct(
            lof_info="GERP_DIST:-2.5,50_BP_RULE:FAIL",
            lof_filter="SINGLE_EXON",
            lof=None,
        )

        ht = hl.Table.parallelize(
            [{"csq": csq_with_missing_lof}],
            hl.tstruct(
                csq=hl.tstruct(lof_info=hl.tstr, lof_filter=hl.tstr, lof=hl.tstr)
            ),
        )

        # Apply the function.
        ht = ht.annotate(updated_csq=update_loftee_end_trunc_filter(ht.csq))

        # Collect results.
        results = ht.collect()

        # This case shouldn't happen. If lof_filter is defined, lof should be defined
        # too. However, we should handle it gracefully by adding END_TRUNC, but
        # maintaining the lof missingness status.
        assert results[0].updated_csq.lof_filter == "END_TRUNC,SINGLE_EXON"
        assert results[0].updated_csq.lof is None

    def test_empty_filter_handling(self):
        """Test handling of empty and None filters."""
        test_cases = [
            hl.Struct(
                lof_info="GERP_DIST:-2.5,50_BP_RULE:FAIL", lof_filter="", lof="HC"
            ),
            hl.Struct(
                lof_info="GERP_DIST:-2.5,50_BP_RULE:FAIL", lof_filter=None, lof="HC"
            ),
        ]

        ht = hl.Table.parallelize(
            [{"csq": csq} for csq in test_cases],
            hl.tstruct(
                csq=hl.tstruct(lof_info=hl.tstr, lof_filter=hl.tstr, lof=hl.tstr)
            ),
        )

        # Apply the function.
        ht = ht.annotate(updated_csq=update_loftee_end_trunc_filter(ht.csq))

        # Collect results.
        results = ht.collect()

        assert results[0].updated_csq.lof_filter == "END_TRUNC"
        assert results[0].updated_csq.lof == "LC"

        assert results[1].updated_csq.lof_filter == "END_TRUNC"
        assert results[1].updated_csq.lof == "LC"


class TestGetVepHelp:
    """Test the get_vep_help function."""

    def test_uses_explicit_path_and_returns_help(self):
        """Test that the given config is read and `vep --help` output returned."""
        config = json.dumps({"command": ["/path/to/vep", "--offline"]})
        with patch(
            "gnomad.utils.vep.hfs.open", _mock_open_returning(config)
        ) as mock_open:
            with patch(
                "gnomad.utils.vep.subprocess.check_output",
                return_value=b"ensembl-vep 110",
            ) as mock_check_output:
                result = get_vep_help("gs://bucket/vep.json")

        assert result == "ensembl-vep 110"
        mock_open.assert_called_once_with("gs://bucket/vep.json")
        # Only the command itself is invoked, not its config flags.
        mock_check_output.assert_called_once_with(["/path/to/vep"])

    def test_falls_back_to_env_var(self, monkeypatch):
        """Test that VEP_CONFIG_URI is used when no path is supplied."""
        monkeypatch.setenv("VEP_CONFIG_URI", "gs://env/vep.json")
        config = json.dumps({"command": ["/path/to/vep"]})
        with patch(
            "gnomad.utils.vep.hfs.open", _mock_open_returning(config)
        ) as mock_open:
            with patch("gnomad.utils.vep.subprocess.check_output", return_value=b"vep"):
                get_vep_help()

        mock_open.assert_called_once_with("gs://env/vep.json")


class TestVepOrLookupVep:
    """Test the vep_or_lookup_vep function."""

    def test_invalid_reference_raises(self):
        """Test that an unrecognized reference build raises a ValueError."""
        with patch("gnomad.utils.vep.get_vep_help", return_value="help"):
            with patch("gnomad.utils.vep.hfs.open", _mock_open_returning("config")):
                with pytest.raises(ValueError, match="Expected one of"):
                    vep_or_lookup_vep(
                        None, reference="hg99", vep_config_path="gs://bucket/vep.json"
                    )

    def test_falls_back_to_full_vep_when_version_unavailable(self):
        """Test that a missing context version VEPs all variants via `hl.vep`.

        When the requested VEP version has no context Table, the function should
        warn and return ``hl.vep(ht, vep_config_path)`` directly. ``hl.vep`` is
        mocked so no VEP process is launched.
        """
        ht = hl.utils.range_table(1)
        sentinel = hl.utils.range_table(2)
        # default_version is set but `versions` is empty, so it is never found.
        context = _fake_vep_context("105", {})
        with patch("gnomad.utils.vep.get_vep_help", return_value="help"):
            with patch("gnomad.utils.vep.hfs.open", _mock_open_returning("config")):
                with patch("gnomad.utils.vep.get_vep_context", return_value=context):
                    with patch(
                        "gnomad.utils.vep.hl.vep", return_value=sentinel
                    ) as mock_vep:
                        result = vep_or_lookup_vep(
                            ht,
                            reference="GRCh38",
                            vep_config_path="gs://bucket/vep.json",
                        )

        # The full-VEP fallback result is returned unchanged.
        assert result is sentinel
        mock_vep.assert_called_once_with(ht, "gs://bucket/vep.json")

    def test_raises_on_vep_help_mismatch(self):
        """Test that a VEP version mismatch between config and context HT raises."""
        context_ht = hl.utils.range_table(1).annotate_globals(
            vep_help="CONTEXT_HELP", vep_config="CONTEXT_CONFIG"
        )
        context = _fake_vep_context("105", {"105": context_ht})
        with patch("gnomad.utils.vep.get_vep_help", return_value="DIFFERENT_HELP"):
            with patch(
                "gnomad.utils.vep.hfs.open", _mock_open_returning("CONTEXT_CONFIG")
            ):
                with patch("gnomad.utils.vep.get_vep_context", return_value=context):
                    with pytest.raises(AssertionError, match="does not match"):
                        vep_or_lookup_vep(
                            hl.utils.range_table(1),
                            reference="GRCh38",
                            vep_config_path="gs://bucket/vep.json",
                        )

    def test_raises_on_vep_config_mismatch(self):
        """Test that a VEP config mismatch between config file and context HT raises."""
        context_ht = hl.utils.range_table(1).annotate_globals(
            vep_help="CONTEXT_HELP", vep_config="CONTEXT_CONFIG"
        )
        context = _fake_vep_context("105", {"105": context_ht})
        # Help matches so the first assert passes; the config read does not.
        with patch("gnomad.utils.vep.get_vep_help", return_value="CONTEXT_HELP"):
            with patch(
                "gnomad.utils.vep.hfs.open", _mock_open_returning("DIFFERENT_CONFIG")
            ):
                with patch("gnomad.utils.vep.get_vep_context", return_value=context):
                    with pytest.raises(AssertionError, match="configuration does"):
                        vep_or_lookup_vep(
                            hl.utils.range_table(1),
                            reference="GRCh38",
                            vep_config_path="gs://bucket/vep.json",
                        )

    def test_uses_context_ht_when_all_variants_present(self):
        """Test the lookup path: variants found in the context HT are not re-VEPed.

        All input keys exist in the context Table, so the re-VEP subset is empty
        and the returned Table carries the context's annotations plus the VEP
        version/help/config globals. ``hl.vep`` is mocked to an identity so no
        VEP process is launched, and it is only applied to the (empty) subset.
        """
        context_ht = hl.utils.range_table(3).annotate(vep=hl.struct(found=True))
        context_ht = context_ht.annotate_globals(
            vep_help="CONTEXT_HELP", vep_config="CONTEXT_CONFIG"
        )
        context = _fake_vep_context("105", {"105": context_ht})
        ht = hl.utils.range_table(3)
        with patch("gnomad.utils.vep.get_vep_help", return_value="CONTEXT_HELP"):
            with patch(
                "gnomad.utils.vep.hfs.open", _mock_open_returning("CONTEXT_CONFIG")
            ):
                with patch("gnomad.utils.vep.get_vep_context", return_value=context):
                    with patch(
                        "gnomad.utils.vep.hl.vep", side_effect=lambda t, p: t
                    ) as mock_vep:
                        result = vep_or_lookup_vep(
                            ht,
                            reference="GRCh38",
                            vep_config_path="gs://bucket/vep.json",
                        )

        # Every variant was annotated from the context HT.
        assert result.count() == 3
        assert result.aggregate(hl.agg.all(result.vep.found))
        # Globals reflect the resolved version and the matched help/config.
        assert hl.eval(result.vep_version) == "v105"
        assert hl.eval(result.vep_help) == "CONTEXT_HELP"
        assert hl.eval(result.vep_config) == "CONTEXT_CONFIG"
        # hl.vep was only ever applied to the (empty) re-VEP subset, never skipped.
        mock_vep.assert_called_once()

    def test_unions_context_lookup_with_revep_of_missing_variants(self):
        """Test the hybrid path: missing variants are re-VEPed and unioned back.

        The context Table covers only a subset of the input keys, so the
        variants absent from it must be VEPed via ``hl.vep`` and unioned with
        the looked-up variants. ``hl.vep`` is mocked to an identity so no VEP
        process is launched; this exercises the non-empty re-VEP/union branch
        that ``test_uses_context_ht_when_all_variants_present`` deliberately
        does not.
        """
        # Context HT holds keys {0, 1}; input has keys {0, 1, 2}, so key 2 is
        # absent from the context and must be re-VEPed.
        context_ht = hl.utils.range_table(2).annotate(vep=hl.struct(found=True))
        context_ht = context_ht.annotate_globals(
            vep_help="CONTEXT_HELP", vep_config="CONTEXT_CONFIG"
        )
        context = _fake_vep_context("105", {"105": context_ht})
        ht = hl.utils.range_table(3)
        with patch("gnomad.utils.vep.get_vep_help", return_value="CONTEXT_HELP"):
            with patch(
                "gnomad.utils.vep.hfs.open", _mock_open_returning("CONTEXT_CONFIG")
            ):
                with patch("gnomad.utils.vep.get_vep_context", return_value=context):
                    with patch(
                        "gnomad.utils.vep.hl.vep", side_effect=lambda t, p: t
                    ) as mock_vep:
                        result = vep_or_lookup_vep(
                            ht,
                            reference="GRCh38",
                            vep_config_path="gs://bucket/vep.json",
                        )

        # Looked-up and re-VEPed variants are unioned back into one Table.
        assert result.count() == 3
        # Only the 2 context-matched variants carry the context annotation; the
        # re-VEPed key 2 keeps a missing vep (identity-mocked hl.vep adds none).
        assert result.aggregate(hl.agg.count_where(result.vep.found)) == 2
        # hl.vep ran exactly once, on the single missing variant, with the config.
        mock_vep.assert_called_once()
        revep_subset, passed_path = mock_vep.call_args[0]
        assert revep_subset.count() == 1
        assert passed_path == "gs://bucket/vep.json"
        # Globals survive the union and reflect the resolved version/help/config.
        assert hl.eval(result.vep_version) == "v105"
        assert hl.eval(result.vep_help) == "CONTEXT_HELP"
        assert hl.eval(result.vep_config) == "CONTEXT_CONFIG"
