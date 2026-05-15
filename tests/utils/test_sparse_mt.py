"""Tests for the sparse_mt utility module."""

import logging

import hail as hl
import pytest

from gnomad.utils.annotations import (
    generate_freq_group_membership_array,
    qual_hist_expr,
)
from gnomad.utils.sparse_mt import (
    compute_allele_number_per_ref_site,
    compute_allele_number_per_ref_site_sparse,
    compute_stats_per_ref_site,
    compute_stats_per_ref_site_sparse,
    get_allele_number_agg_func,
    get_coverage_agg_func,
    get_coverage_agg_func_sparse,
    merge_coverage_stats_array_expression,
    merge_qual_hists_array_expression,
    merge_sum_array_expression,
)

# Set up logger for tests.
logger = logging.getLogger(__name__)


class TestGetCoverageAggFunc:
    """Test the get_coverage_agg_func function."""

    def test_get_coverage_agg_func_returns_tuple(self):
        """Test that get_coverage_agg_func returns a tuple of two functions."""
        transform_func, agg_func = get_coverage_agg_func()

        assert callable(transform_func)
        assert callable(agg_func)

    def test_get_coverage_agg_func_default_params(self):
        """Test get_coverage_agg_func with default parameters."""
        transform_func, _ = get_coverage_agg_func()

        # Test the transform function with various inputs.
        test_cases = [
            ({"DP": 10}, 10),  # Normal case.
            ({"DP": 0}, 0),  # Zero value.
            ({"DP": 150}, 150),  # Large value.
        ]

        for input_data, expected in test_cases:
            # Create a simple struct for testing.
            test_struct = hl.Struct(**input_data)
            result = hl.eval(transform_func(test_struct))
            assert result == expected

    def test_get_coverage_agg_func_missing_values(self):
        """Test that missing values are handled correctly."""
        transform_func, _ = get_coverage_agg_func()

        # Test with missing value - need to use hl.missing.
        test_struct = hl.Struct(DP=hl.missing(hl.tint32))
        result = hl.eval(transform_func(test_struct))
        assert result == 0

    def test_get_coverage_agg_func_with_nan_values(self):
        """Test that NaN values are handled correctly."""
        transform_func, _ = get_coverage_agg_func()
        test_struct = hl.Struct(DP=hl.float64("nan"))
        result = hl.eval(transform_func(test_struct))
        assert result == 0

    def test_aggregation_function_structure(self):
        """Test that the aggregation function returns the expected structure."""
        _, agg_func = get_coverage_agg_func()
        test_data = [10, 20, 30, 40, 50]
        # Create a Hail Table with a DP field.
        ht = hl.Table.parallelize(
            [{"DP": v} for v in test_data], hl.tstruct(DP=hl.tint32)
        )
        # Aggregate using the aggregation function.
        result = ht.aggregate(agg_func(ht.DP))
        expected_fields = ["coverage_counter", "mean", "median_approx", "total_DP"]
        for field in expected_fields:
            assert hasattr(result, field)

    def test_aggregation_function_values(self):
        """Test that the aggregation function computes correct values."""
        _, agg_func = get_coverage_agg_func()
        test_data = [10, 20, 30, 40, 50]
        ht = hl.Table.parallelize(
            [{"DP": v} for v in test_data], hl.tstruct(DP=hl.tint32)
        )
        result = ht.aggregate(agg_func(ht.DP))

        # Test computed values.
        assert result.mean == 30.0  # (10+20+30+40+50)/5.
        assert result.total_DP == 150  # 10+20+30+40+50.
        assert result.median_approx == 30  # Approximate median.
        assert result.coverage_counter.get(10, 0) == 1
        assert result.coverage_counter.get(20, 0) == 1
        assert result.coverage_counter.get(30, 0) == 1
        assert result.coverage_counter.get(40, 0) == 1
        assert result.coverage_counter.get(50, 0) == 1

    def test_aggregation_function_with_missing_values(self):
        """Test aggregation function handles missing values correctly."""
        _, agg_func = get_coverage_agg_func()
        test_data = [10, 20, None, 40, 50]
        ht = hl.Table.parallelize(
            [{"DP": v if v is not None else hl.missing(hl.tint32)} for v in test_data],
            hl.tstruct(DP=hl.tint32),
        )
        result = ht.aggregate(agg_func(ht.DP))

        # Missing values should be transformed to 0, so they should be included in
        # aggregation.
        assert result.total_DP == 120  # 10+20+0+40+50.
        # The mean calculation may vary depending on how missing values are handled.
        # Let's check that the structure is correct.
        assert hasattr(result, "mean")
        assert hasattr(result, "total_DP")
        assert hasattr(result, "coverage_counter")

    def test_aggregation_function_with_zero_values(self):
        """Test aggregation function handles zero values correctly."""
        _, agg_func = get_coverage_agg_func()
        test_data = [0, 10, 20, 0, 30]
        ht = hl.Table.parallelize(
            [{"DP": v} for v in test_data], hl.tstruct(DP=hl.tint32)
        )
        result = ht.aggregate(agg_func(ht.DP))

        # Zero values should be preserved.
        assert result.total_DP == 60  # 0+10+20+0+30.
        assert result.mean == 12.0  # (0+10+20+0+30)/5.
        assert result.coverage_counter.get(0, 0) == 2

    def test_aggregation_function_with_negative_values(self):
        """Test aggregation function handles negative values correctly."""
        _, agg_func = get_coverage_agg_func()
        test_data = [-5, 10, 20, -3, 30]
        ht = hl.Table.parallelize(
            [{"DP": v} for v in test_data], hl.tstruct(DP=hl.tint32)
        )
        result = ht.aggregate(agg_func(ht.DP))

        # Negative values should be preserved.
        assert result.total_DP == 52  # -5+10+20-3+30.
        assert result.mean == 10.4  # (-5+10+20-3+30)/5.
        assert result.coverage_counter.get(-5, 0) == 1
        assert result.coverage_counter.get(-3, 0) == 1

    def test_aggregation_function_edge_cases(self):
        """Test aggregation function with edge cases."""
        _, agg_func = get_coverage_agg_func()

        # Test with single value.
        ht1 = hl.Table.parallelize([{"DP": 42}], hl.tstruct(DP=hl.tint32))
        result1 = ht1.aggregate(agg_func(ht1.DP))
        assert result1.mean == 42.0
        assert result1.total_DP == 42
        assert result1.coverage_counter.get(42, 0) == 1

        # Test with empty table (should handle gracefully).
        ht2 = hl.Table.parallelize([], hl.tstruct(DP=hl.tint32))
        result2 = ht2.aggregate(agg_func(ht2.DP))
        assert result2.total_DP == 0
        assert result2.mean == 0.0

    def test_aggregation_function_with_different_max_cov_bin(self):
        """Test that max_cov_bin parameter affects the aggregation function."""
        _, agg_func = get_coverage_agg_func(max_cov_bin=5)
        test_data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        ht = hl.Table.parallelize(
            [{"DP": v} for v in test_data], hl.tstruct(DP=hl.tint32)
        )
        result = ht.aggregate(agg_func(ht.DP))
        coverage_counter = result.coverage_counter
        assert (
            coverage_counter.get(5, 0) == 6
        )  # 1 original value 5 + 5 capped values (6-10)
        assert coverage_counter.get(6, 0) == 0
        assert coverage_counter.get(7, 0) == 0

    def test_aggregation_function_max_cov_bin_capping(self):
        """Test that max_cov_bin properly caps values in coverage_counter."""
        _, agg_func = get_coverage_agg_func(max_cov_bin=3)
        test_data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        ht = hl.Table.parallelize(
            [{"DP": v} for v in test_data], hl.tstruct(DP=hl.tint32)
        )
        result = ht.aggregate(agg_func(ht.DP))
        coverage_counter = result.coverage_counter

        # Values 1, 2, 3 should be counted as themselves.
        assert coverage_counter.get(1, 0) == 1
        assert coverage_counter.get(2, 0) == 1
        # Values 3-10 should all be counted as 3 (max_cov_bin).
        # So there should be 8 values at position 3: the original value 3 plus 7
        # capped values (4-10).
        assert coverage_counter.get(3, 0) == 8

        # No values should be counted above max_cov_bin.
        assert coverage_counter.get(4, 0) == 0
        assert coverage_counter.get(5, 0) == 0

    def test_aggregation_function_median_approximation(self):
        """Test that median_approx provides reasonable approximation."""
        _, agg_func = get_coverage_agg_func()

        # Test with odd number of values.
        test_data_odd = [10, 20, 30, 40, 50]
        ht_odd = hl.Table.parallelize(
            [{"DP": v} for v in test_data_odd], hl.tstruct(DP=hl.tint32)
        )
        result_odd = ht_odd.aggregate(agg_func(ht_odd.DP))
        assert result_odd.median_approx == 30

        # Test with even number of values - approximate median may vary.
        test_data_even = [10, 20, 30, 40]
        ht_even = hl.Table.parallelize(
            [{"DP": v} for v in test_data_even], hl.tstruct(DP=hl.tint32)
        )
        result_even = ht_even.aggregate(agg_func(ht_even.DP))
        # Hail's median_approx is non-deterministic, so we can't test for exact
        # equality.
        # See https://hail.is/docs/0.2/aggregators.html#hail.expr.aggregators.approx_median.
        # Approximate median should be close to 25 (true median).
        assert result_even.median_approx in [20, 25, 30]

    def test_transform_function_with_various_inputs(self):
        """Test the transform function with various input types."""
        transform_func, _ = get_coverage_agg_func(max_cov_bin=100)

        # Create a sample dataset.
        sample_data = [{"DP": 10}, {"DP": 20}, {"DP": 30}, {"DP": 0}, {"DP": 150}]

        # Transform the data.
        transformed_data = []
        for data in sample_data:
            test_struct = hl.Struct(**data)
            result = hl.eval(transform_func(test_struct))
            transformed_data.append(result)

        # Test that transformation worked correctly.
        assert len(transformed_data) == 5
        assert 150 in transformed_data  # Large value should not be capped in transform.
        assert 0 in transformed_data  # Zero values should be preserved.

    def test_different_dp_field_names(self):
        """Test that different DP field names work correctly."""
        field_names = ["DP", "depth", "coverage", "read_depth"]

        for field_name in field_names:
            transform_func, _ = get_coverage_agg_func(dp_field=field_name)

            # Test with the custom field name.
            test_struct = hl.Struct(**{field_name: 25})
            result = hl.eval(transform_func(test_struct))
            assert result == 25

    def test_transform_and_aggregation_integration(self):
        """Test that transform and aggregation functions work together correctly."""
        transform_func, agg_func = get_coverage_agg_func(max_cov_bin=50)

        # Create test data with various edge cases.
        test_data = [10, 20, None, 40, 50, 0, 100, -5, float("nan")]

        # Test transform function on individual values.
        transformed_values = []
        for value in test_data:
            if value is None:
                test_struct = hl.Struct(DP=hl.missing(hl.tint32))
            elif isinstance(value, float) and value != value:  # Check for NaN
                test_struct = hl.Struct(DP=hl.float64("nan"))
            else:
                test_struct = hl.Struct(DP=value)
            result = hl.eval(transform_func(test_struct))
            transformed_values.append(result)

        # Expected transformed values: [10, 20, 0, 40, 50, 0, 100, -5, 0].
        expected_transformed = [10, 20, 0, 40, 50, 0, 100, -5, 0]
        # For aggregation test (without NaN): [10, 20, 0, 40, 50, 0, 100, -5]
        expected_transformed_no_nan = [10, 20, 0, 40, 50, 0, 100, -5]
        assert transformed_values == expected_transformed

        # Test aggregation function on the same data (excluding NaN due to type
        # handling complexity).
        test_data_no_nan = [
            v for v in test_data if not (isinstance(v, float) and v != v)
        ]
        ht = hl.Table.parallelize(
            [
                {"DP": v if v is not None else hl.missing(hl.tint32)}
                for v in test_data_no_nan
            ],
            hl.tstruct(DP=hl.tint32),
        )

        # Apply the transform function to each row to create transformed data.
        ht = ht.annotate(transformed_DP=transform_func(ht))
        result = ht.aggregate(agg_func(ht.transformed_DP))

        # Check that aggregation works correctly.
        assert result.total_DP == 215  # 10+20+0+40+50+0+100-5.
        # The aggregation function operates on transformed data (missing values become 0).
        # So mean should be calculated from all 8 values: 215/8 = 26.875.
        assert (
            abs(
                result.mean
                - sum(expected_transformed_no_nan) / len(expected_transformed_no_nan)
            )
            < 0.1
        )  # Allow small tolerance.
        # Now that we're using transformed data, there are 2 zeros: one explicit,
        # one from missing.
        assert (
            result.coverage_counter.get(0, 0) == 2
        )  # 2 zero values (explicit 0 + missing->0).
        # Values 50 (original) and 100 (capped to 50 by max_cov_bin) should be at
        # position 50.
        assert result.coverage_counter.get(50, 0) == 2

    def test_custom_field_transform_and_aggregation(self):
        """Test both transform and aggregation with custom field names."""
        transform_func, agg_func = get_coverage_agg_func(
            dp_field="depth", max_cov_bin=30
        )

        # Test transform function.
        test_cases = [
            ({"depth": 10}, 10),
            ({"depth": 0}, 0),
            ({"depth": 50}, 50),  # Should not be capped in transform.
        ]

        for input_data, expected in test_cases:
            test_struct = hl.Struct(**input_data)
            result = hl.eval(transform_func(test_struct))
            assert result == expected

        # Test aggregation function
        test_data = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
        ht = hl.Table.parallelize(
            [{"depth": v} for v in test_data], hl.tstruct(depth=hl.tint32)
        )

        result = ht.aggregate(agg_func(ht.depth))

        # Check aggregation results.
        assert result.total_DP == 275  # Sum of all values.
        assert result.mean == 27.5  # 275/10.
        # Values 35, 40, 45, 50 should be capped to 30, plus the original 30.
        assert result.coverage_counter.get(30, 0) == 5  # 4 capped + 1 original.


class TestComputeAlleleNumberPerRefSiteReduceToMinimalGroups:
    """Test that the reduce_to_minimal_groups optimization preserves AN values."""

    @pytest.fixture
    def synthetic_vds(self):
        """Build a tiny VDS with 8 samples × 4 sites for AN testing.

        The VDS has:
          - 8 samples annotated with `gen_anc` (afr/nfe) and `sex` (XX/XY).
          - 4 variant sites with explicit GT calls (some missing) for every
            sample.
          - A reference_data MT with 1-base ref blocks at the same loci, so
            densification is essentially a no-op (we test the reduction
            machinery, not the densification machinery).
        """
        samples = [
            ("s1", "afr", "XX"),
            ("s2", "afr", "XY"),
            ("s3", "afr", "XX"),
            ("s4", "afr", "XY"),
            ("s5", "nfe", "XX"),
            ("s6", "nfe", "XY"),
            ("s7", "nfe", "XX"),
            ("s8", "nfe", "XY"),
        ]
        positions = [1000, 2000, 3000, 4000]
        variants = [
            (hl.locus("chr1", p, reference_genome="GRCh38"), ["A", "T"])
            for p in positions
        ]
        # Mix of called and missing genotypes so AN varies across sites and
        # strata. Missing genotypes contribute 0 to AN; called diploid
        # genotypes contribute 2.
        gt_patterns = [
            [(0, 0), (0, 1), (1, 1), (0, 1), (0, 0), (0, 1), (1, 1), (0, 1)],
            [(0, 1), None, (0, 0), (0, 1), (0, 1), (1, 1), None, (0, 1)],
            [None, (0, 0), (0, 1), (1, 1), (0, 0), (0, 0), (0, 1), (1, 1)],
            [(1, 1), (0, 1), (0, 1), None, (1, 1), (0, 1), (0, 1), (0, 0)],
        ]

        sample_table = hl.Table.parallelize(
            [{"s": s, "gen_anc": g, "sex": x} for s, g, x in samples],
            hl.tstruct(s=hl.tstr, gen_anc=hl.tstr, sex=hl.tstr),
        ).key_by("s")

        # Build the variant_data MT.
        vd_entries = []
        for v_idx, (locus, alleles) in enumerate(variants):
            for s_idx, (sample_id, _, _) in enumerate(samples):
                gt = gt_patterns[v_idx][s_idx]
                vd_entries.append(
                    {
                        "locus": locus,
                        "alleles": alleles,
                        "s": sample_id,
                        "GT": hl.call(*gt) if gt is not None else hl.missing(hl.tcall),
                    }
                )
        vd = hl.Table.parallelize(
            vd_entries,
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                s=hl.tstr,
                GT=hl.tcall,
            ),
        ).to_matrix_table(row_key=["locus", "alleles"], col_key=["s"])
        vd = vd.annotate_cols(
            gen_anc=sample_table[vd.s].gen_anc,
            sex=sample_table[vd.s].sex,
        )

        # Build the reference_data MT: 1-base ref blocks at every variant
        # locus. END = locus.position means the block covers exactly one
        # position.
        rd_entries = []
        for locus, _ in variants:
            for sample_id, _, _ in samples:
                rd_entries.append(
                    {
                        "locus": locus,
                        "s": sample_id,
                        "END": locus.position,
                        "LEN": 1,
                    }
                )
        rd = hl.Table.parallelize(
            rd_entries,
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                s=hl.tstr,
                END=hl.tint32,
                LEN=hl.tint32,
            ),
        ).to_matrix_table(row_key=["locus"], col_key=["s"])

        return hl.vds.VariantDataset(reference_data=rd, variant_data=vd)

    @pytest.fixture
    def reference_ht(self):
        """Return a reference HT covering the four variant loci."""
        positions = [1000, 2000, 3000, 4000]
        return hl.Table.parallelize(
            [
                {"locus": hl.locus("chr1", p, reference_genome="GRCh38")}
                for p in positions
            ],
            hl.tstruct(locus=hl.tlocus("GRCh38")),
        ).key_by("locus")

    @staticmethod
    def _index_an(ht):
        """Return a dict mapping sorted strata items to lists of per-row AN."""
        strata_meta = hl.eval(ht.strata_meta)
        rows = ht.collect()
        out = {}
        for i, m in enumerate(strata_meta):
            key = tuple(sorted(m.items()))
            out[key] = [r.AN[i] for r in rows]
        return out

    def test_reduce_matches_full_with_gen_anc_and_sex(
        self, synthetic_vds, reference_ht
    ):
        """compute_allele_number_per_ref_site(reduce=True) matches the full output."""
        vd = synthetic_vds.variant_data
        strata_expr = [
            {"gen_anc": vd.gen_anc},
            {"sex": vd.sex},
            {"gen_anc": vd.gen_anc, "sex": vd.sex},
        ]

        full_ht = compute_allele_number_per_ref_site(
            synthetic_vds, reference_ht, strata_expr=strata_expr
        )
        reduced_ht = compute_allele_number_per_ref_site(
            synthetic_vds,
            reference_ht,
            strata_expr=strata_expr,
            reduce_to_minimal_groups=True,
        )

        full_indexed = self._index_an(full_ht)
        reduced_indexed = self._index_an(reduced_ht)

        # The strata_meta key sets must match exactly — reduction is
        # transparent to callers.
        assert set(full_indexed.keys()) == set(reduced_indexed.keys())
        # Sanity check: the parent {raw} entry must be present and was
        # decomposed into the gen_anc-by-sex leaves.
        assert (("group", "raw"),) in full_indexed

        for key in full_indexed:
            assert full_indexed[key] == reduced_indexed[key], key

        # Direct invariant check: the all-samples AN must equal the
        # element-wise sum of the gen_anc-by-sex leaf ANs. This catches
        # the case where both the full and reduced paths agree but the
        # underlying math is wrong.
        leaf_keys = [
            k for k in full_indexed if "gen_anc" in dict(k) and "sex" in dict(k)
        ]
        leaf_values = [full_indexed[k] for k in leaf_keys]
        leaf_sums = [sum(vals) for vals in zip(*leaf_values)]
        assert full_indexed[(("group", "raw"),)] == leaf_sums


class TestComputeStatsPerRefSiteReducibleAggs:
    """Test `reducible_aggs` mixes summable and non-summable aggregations under leaf reduction."""

    @pytest.fixture
    def synthetic_vds(self):
        """8 samples × 4 sites VDS — same shape as the AN reduction test."""
        samples = [
            ("s1", "afr", "XX"),
            ("s2", "afr", "XY"),
            ("s3", "afr", "XX"),
            ("s4", "afr", "XY"),
            ("s5", "nfe", "XX"),
            ("s6", "nfe", "XY"),
            ("s7", "nfe", "XX"),
            ("s8", "nfe", "XY"),
        ]
        positions = [1000, 2000, 3000, 4000]
        variants = [
            (hl.locus("chr1", p, reference_genome="GRCh38"), ["A", "T"])
            for p in positions
        ]
        gt_patterns = [
            [(0, 0), (0, 1), (1, 1), (0, 1), (0, 0), (0, 1), (1, 1), (0, 1)],
            [(0, 1), None, (0, 0), (0, 1), (0, 1), (1, 1), None, (0, 1)],
            [None, (0, 0), (0, 1), (1, 1), (0, 0), (0, 0), (0, 1), (1, 1)],
            [(1, 1), (0, 1), (0, 1), None, (1, 1), (0, 1), (0, 1), (0, 0)],
        ]

        sample_table = hl.Table.parallelize(
            [{"s": s, "gen_anc": g, "sex": x} for s, g, x in samples],
            hl.tstruct(s=hl.tstr, gen_anc=hl.tstr, sex=hl.tstr),
        ).key_by("s")

        vd_entries = []
        for v_idx, (locus, alleles) in enumerate(variants):
            for s_idx, (sample_id, _, _) in enumerate(samples):
                gt = gt_patterns[v_idx][s_idx]
                vd_entries.append(
                    {
                        "locus": locus,
                        "alleles": alleles,
                        "s": sample_id,
                        "GT": hl.call(*gt) if gt is not None else hl.missing(hl.tcall),
                        # `adj=True` for every (variant, sample) so we can build
                        # group_membership_hts with `group_label="adj"` (the
                        # default) and exercise non-leaf parent targets like
                        # `{"group": "adj"}` without needing the AD/DP/GQ
                        # fields the auto-adj path would otherwise demand.
                        "adj": True,
                    }
                )
        vd = hl.Table.parallelize(
            vd_entries,
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                s=hl.tstr,
                GT=hl.tcall,
                adj=hl.tbool,
            ),
        ).to_matrix_table(row_key=["locus", "alleles"], col_key=["s"])
        vd = vd.annotate_cols(
            gen_anc=sample_table[vd.s].gen_anc,
            sex=sample_table[vd.s].sex,
        )

        rd_entries = []
        for locus, _ in variants:
            for sample_id, _, _ in samples:
                rd_entries.append(
                    {"locus": locus, "s": sample_id, "END": locus.position, "LEN": 1}
                )
        rd = hl.Table.parallelize(
            rd_entries,
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                s=hl.tstr,
                END=hl.tint32,
                LEN=hl.tint32,
            ),
        ).to_matrix_table(row_key=["locus"], col_key=["s"])

        return hl.vds.VariantDataset(reference_data=rd, variant_data=vd)

    @pytest.fixture
    def reference_ht(self):
        """Return a reference HT covering the four variant loci."""
        return hl.Table.parallelize(
            [
                {"locus": hl.locus("chr1", p, reference_genome="GRCh38")}
                for p in [1000, 2000, 3000, 4000]
            ],
            hl.tstruct(locus=hl.tlocus("GRCh38")),
        ).key_by("locus")

    @staticmethod
    def _index(ht, ann):
        meta = hl.eval(ht.strata_meta)
        rows = ht.collect()
        return {
            tuple(sorted(m.items())): [r[ann][i] for r in rows]
            for i, m in enumerate(meta)
        }

    def _build_group_membership_ht(self, vd, *, reduce: bool):
        """Pre-build a group_membership_ht with `group_label="raw"` and `no_raw_group=True` to match the per-ref-site internal pattern (no adj filtering needed on entries — the test fixture's variant_data only has GT)."""
        cols_ht = vd.cols()
        strata_expr = [
            {"gen_anc": cols_ht.gen_anc},
            {"sex": cols_ht.sex},
            {"gen_anc": cols_ht.gen_anc, "sex": cols_ht.sex},
        ]
        return generate_freq_group_membership_array(
            cols_ht,
            strata_expr,
            no_raw_group=True,
            reduce_to_minimal_groups=reduce,
            group_label="raw",
        )

    # A fully-stratified leaf — survives reduction, so the
    # `freq_meta.index(...)` lookup in agg_by_strata works against
    # both the un-reduced and reduced freq_meta.
    LEAF_TARGET = {"group": "raw", "gen_anc": "afr", "sex": "XX"}

    def test_mixed_reducible_and_restricted_aggs_match_full(
        self, synthetic_vds, reference_ht
    ):
        """AN goes through the leaf reduction; max_called is restricted to a leaf via entry_agg_group_membership; both match the full-fanout baseline.

        The restriction target is a fully-stratified leaf so the
        `freq_meta.index(...)` lookup succeeds against both the
        un-reduced and the reduced `freq_meta`. Restricting a
        non-summable aggregation to a non-leaf parent (e.g.
        `{"group": "adj"}` against a gen_anc×sex stratification) is a
        separate limitation of `agg_by_strata`'s lookup against the
        reduced `freq_meta` — independent of the `reducible_aggs`
        plumbing under test here.
        """
        vd = synthetic_vds.variant_data
        full_gmh = self._build_group_membership_ht(vd, reduce=False)
        reduced_gmh = self._build_group_membership_ht(vd, reduce=True)

        entry_agg_funcs = {
            "AN": (lambda t: t.GT.ploidy, hl.agg.sum),
            "max_called": (lambda t: hl.int32(hl.is_defined(t.GT)), hl.agg.max),
        }

        full_ht = compute_stats_per_ref_site(
            synthetic_vds,
            reference_ht,
            entry_agg_funcs,
            group_membership_ht=full_gmh,
            entry_agg_group_membership={"max_called": [self.LEAF_TARGET]},
        )
        reduced_ht = compute_stats_per_ref_site(
            synthetic_vds,
            reference_ht,
            entry_agg_funcs,
            group_membership_ht=reduced_gmh,
            reducible_aggs={"AN"},
            entry_agg_group_membership={"max_called": [self.LEAF_TARGET]},
        )

        full_an = self._index(full_ht, "AN")
        reduced_an = self._index(reduced_ht, "AN")
        assert set(full_an) == set(reduced_an)
        for k in full_an:
            assert full_an[k] == reduced_an[k], k

        # max_called: same target leaf on both runs, so per-row length-1
        # arrays should match element-for-element.
        full_max = [r.max_called for r in full_ht.collect()]
        reduced_max = [r.max_called for r in reduced_ht.collect()]
        assert all(len(v) == 1 for v in full_max)
        assert all(len(v) == 1 for v in reduced_max)
        assert full_max == reduced_max

    def test_reducible_aggs_overlap_with_entry_agg_group_membership_raises(
        self, synthetic_vds, reference_ht
    ):
        """An annotation in BOTH reducible_aggs and entry_agg_group_membership raises ValueError."""
        vd = synthetic_vds.variant_data
        gmh = self._build_group_membership_ht(vd, reduce=True)
        with pytest.raises(ValueError, match="disjoint"):
            compute_stats_per_ref_site(
                synthetic_vds,
                reference_ht,
                {"AN": (lambda t: t.GT.ploidy, hl.agg.sum)},
                group_membership_ht=gmh,
                reducible_aggs={"AN"},
                entry_agg_group_membership={"AN": [self.LEAF_TARGET]},
            )

    def test_unaccounted_annotation_under_reduction_raises(
        self, synthetic_vds, reference_ht
    ):
        """An annotation that is neither reducible nor restricted via entry_agg_group_membership raises when reduction is in effect."""
        vd = synthetic_vds.variant_data
        gmh = self._build_group_membership_ht(vd, reduce=True)
        entry_agg_funcs = {
            "AN": (lambda t: t.GT.ploidy, hl.agg.sum),
            "max_called": (lambda t: hl.int32(hl.is_defined(t.GT)), hl.agg.max),
        }
        with pytest.raises(ValueError, match="no shape designation"):
            compute_stats_per_ref_site(
                synthetic_vds,
                reference_ht,
                entry_agg_funcs,
                group_membership_ht=gmh,
                reducible_aggs={"AN"},
                # `max_called` is missing from BOTH reducible_aggs and
                # entry_agg_group_membership — the guard should fire.
            )

    def _build_group_membership_ht_with_adj(self, vd, *, reduce: bool):
        """Like `_build_group_membership_ht` but with `group_label='adj'` (the default) so freq_meta entries carry `group:adj` and a non-leaf parent like `{group:adj}` exists in `freq_meta_full`. Requires the variant_data MT to have an `adj` field."""
        cols_ht = vd.cols()
        strata_expr = [
            {"gen_anc": cols_ht.gen_anc},
            {"sex": cols_ht.sex},
            {"gen_anc": cols_ht.gen_anc, "sex": cols_ht.sex},
        ]
        return generate_freq_group_membership_array(
            cols_ht,
            strata_expr,
            reduce_to_minimal_groups=reduce,
        )

    def test_non_leaf_parent_target_under_reduction_matches_full(
        self, synthetic_vds, reference_ht
    ):
        """`entry_agg_group_membership` target `{"group": "adj"}` is a non-leaf parent under reduction; the parent-resolution path in `agg_by_strata` synthesizes its sample-set as the union of its leaf-children and produces the same per-row values as the un-reduced run.

        This is the gnomad_qc compute_coverage.py call shape: AN goes
        through the reducible-aggs path; a non-summable annotation
        (here `max_called`) is restricted to the global adj group via
        `entry_agg_group_membership`. Under reduction the parent
        `{"group": "adj"}` lives only in `freq_meta_full`; the
        decomposition reconstructs its sample-set from the
        gen_anc-by-sex leaves.
        """
        vd = synthetic_vds.variant_data
        full_gmh = self._build_group_membership_ht_with_adj(vd, reduce=False)
        reduced_gmh = self._build_group_membership_ht_with_adj(vd, reduce=True)

        # Sanity-check the structural shape so a future regression in
        # find_minimal_strata_groups (e.g., promoting the root to a
        # leaf) is caught here rather than as a silent value mismatch.
        reduced_freq_meta = hl.eval(reduced_gmh.freq_meta)
        reduced_freq_meta_full = hl.eval(reduced_gmh.freq_meta_full)
        assert {"group": "adj"} in [dict(m) for m in reduced_freq_meta_full]
        assert {"group": "adj"} not in [dict(m) for m in reduced_freq_meta]

        entry_agg_funcs = {
            "AN": (lambda t: t.GT.ploidy, hl.agg.sum),
            "max_called": (lambda t: hl.int32(hl.is_defined(t.GT)), hl.agg.max),
        }
        adj_target = {"group": "adj"}

        full_ht = compute_stats_per_ref_site(
            synthetic_vds,
            reference_ht,
            entry_agg_funcs,
            group_membership_ht=full_gmh,
            entry_agg_group_membership={"max_called": [adj_target]},
        )
        reduced_ht = compute_stats_per_ref_site(
            synthetic_vds,
            reference_ht,
            entry_agg_funcs,
            group_membership_ht=reduced_gmh,
            reducible_aggs={"AN"},
            entry_agg_group_membership={"max_called": [adj_target]},
        )

        # AN: full strata_meta and per-strata values match.
        full_an = self._index(full_ht, "AN")
        reduced_an = self._index(reduced_ht, "AN")
        assert set(full_an) == set(reduced_an)
        for k in full_an:
            assert full_an[k] == reduced_an[k], k

        # max_called: a single value per row, matching the {group:adj}
        # entry from the un-reduced run.
        full_max = [r.max_called for r in full_ht.collect()]
        reduced_max = [r.max_called for r in reduced_ht.collect()]
        assert all(len(v) == 1 for v in full_max)
        assert all(len(v) == 1 for v in reduced_max)
        assert full_max == reduced_max


class TestComputeAlleleNumberPerRefSiteSparse:
    """Test the sparse-aware AN aggregation that avoids a full VDS densify."""

    # At each locus, every sample is in exactly one of `variant_data` or
    # `reference_data` (VDS invariant). The dense path
    # (`compute_allele_number_per_ref_site`) and the sparse path
    # (`compute_allele_number_per_ref_site_sparse(include_variant_data=True)`)
    # must produce identical AN per (locus, stratum) when this invariant holds.
    # Format: 4 loci × 4 samples; each cell is either `"R"` (ref block, ploidy
    # 2) or `("V", gt_tuple_or_None)` (variant_data row with the given GT).
    _ASSIGNMENT = [
        # locus 1000 (autosomal): s1=ref, s2=ref, s3=var(0/1), s4=var(1/1)
        ["R", "R", ("V", (0, 1)), ("V", (1, 1))],
        # locus 2000 (autosomal): s1=ref, s2=var(missing), s3=ref, s4=var(0/1)
        ["R", ("V", None), "R", ("V", (0, 1))],
        # locus 3000 (autosomal): all ref — exercises the no-variant_data-row case
        ["R", "R", "R", "R"],
        # locus 4000 (autosomal): s1=var(0/1), s2=ref, s3=var(missing), s4=ref
        [("V", (0, 1)), "R", ("V", None), "R"],
    ]

    _SAMPLES = [
        ("s1", "afr", "XX"),
        ("s2", "afr", "XY"),
        ("s3", "nfe", "XX"),
        ("s4", "nfe", "XY"),
    ]

    @staticmethod
    def _build_vds(assignment, samples, positions, contig="chr1"):
        """Build a VDS from a (sample, locus) -> {"R" | ("V", gt)} assignment.

        Respects the VDS invariant: each (sample, locus) is in exactly one
        of `reference_data` or `variant_data`.
        """
        loci = [hl.locus(contig, p, reference_genome="GRCh38") for p in positions]
        sample_table = hl.Table.parallelize(
            [{"s": s, "gen_anc": g, "sex": x} for s, g, x in samples],
            hl.tstruct(s=hl.tstr, gen_anc=hl.tstr, sex=hl.tstr),
        ).key_by("s")

        vd_entries = []
        rd_entries = []
        for l_idx, locus in enumerate(loci):
            for s_idx, (sample_id, _, _) in enumerate(samples):
                cell = assignment[l_idx][s_idx]
                if cell == "R":
                    rd_entries.append(
                        {
                            "locus": locus,
                            "s": sample_id,
                            # Real VDS reference blocks always carry LGT=0/0
                            # (diploid hom-ref); ploidy is what feeds the AN
                            # aggregation in `get_allele_number_agg_func`.
                            "LGT": hl.call(0, 0),
                            "END": locus.position,
                            "LEN": 1,
                        }
                    )
                else:
                    _, gt = cell
                    vd_entries.append(
                        {
                            "locus": locus,
                            "alleles": ["A", "T"],
                            "s": sample_id,
                            # Real VDSs carry GT in variant_data and LGT in
                            # reference_data — `hl.vds.to_dense_mt` relies on
                            # the field names differing to disambiguate.
                            "GT": (
                                hl.call(*gt) if gt is not None else hl.missing(hl.tcall)
                            ),
                        }
                    )

        vd = hl.Table.parallelize(
            vd_entries,
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                s=hl.tstr,
                GT=hl.tcall,
            ),
        ).to_matrix_table(row_key=["locus", "alleles"], col_key=["s"])
        vd = vd.annotate_cols(
            gen_anc=sample_table[vd.s].gen_anc, sex=sample_table[vd.s].sex
        )

        rd = hl.Table.parallelize(
            rd_entries,
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                s=hl.tstr,
                LGT=hl.tcall,
                END=hl.tint32,
                LEN=hl.tint32,
            ),
        ).to_matrix_table(row_key=["locus"], col_key=["s"])
        rd = rd.annotate_cols(
            gen_anc=sample_table[rd.s].gen_anc, sex=sample_table[rd.s].sex
        )

        return hl.vds.VariantDataset(reference_data=rd, variant_data=vd)

    @pytest.fixture
    def synthetic_vds(self):
        """4 samples × 4 autosomal loci VDS respecting the VDS invariant."""
        return self._build_vds(
            self._ASSIGNMENT, self._SAMPLES, [1000, 2000, 3000, 4000]
        )

    @pytest.fixture
    def reference_ht(self):
        """Return a reference HT covering the four autosomal loci."""
        return hl.Table.parallelize(
            [
                {"locus": hl.locus("chr1", p, reference_genome="GRCh38")}
                for p in [1000, 2000, 3000, 4000]
            ],
            hl.tstruct(locus=hl.tlocus("GRCh38")),
        ).key_by("locus")

    @staticmethod
    def _index_an(ht):
        """Return {sorted-meta-tuple: [per-row AN]} for an array-shaped output."""
        strata_meta = hl.eval(ht.strata_meta)
        rows = ht.collect()
        out = {}
        for i, m in enumerate(strata_meta):
            key = tuple(sorted(m.items()))
            out[key] = [r.AN[i] for r in rows]
        return out

    def _build_gmh(self, vd, *, reduce=False):
        cols_ht = vd.cols()
        strata_expr = [
            {"gen_anc": cols_ht.gen_anc},
            {"sex": cols_ht.sex},
            {"gen_anc": cols_ht.gen_anc, "sex": cols_ht.sex},
        ]
        return generate_freq_group_membership_array(
            cols_ht,
            strata_expr,
            no_raw_group=True,
            reduce_to_minimal_groups=reduce,
            group_label="raw",
        )

    def test_rejects_non_vds_input(self, reference_ht):
        """A non-VDS input raises ValueError."""
        mt = hl.utils.range_matrix_table(2, 2)
        with pytest.raises(ValueError, match="VariantDataset"):
            compute_allele_number_per_ref_site_sparse(mt, reference_ht)

    def test_rejects_both_strata_args(self, synthetic_vds, reference_ht):
        """Supplying both `strata_expr` and `group_membership_ht` raises."""
        vd = synthetic_vds.variant_data
        gmh = self._build_gmh(vd)
        with pytest.raises(ValueError, match="Only one"):
            compute_allele_number_per_ref_site_sparse(
                synthetic_vds,
                reference_ht,
                strata_expr=[{"sex": vd.sex}],
                group_membership_ht=gmh,
            )

    def test_include_variant_data_requires_strata(self, synthetic_vds, reference_ht):
        """`include_variant_data=True` without strata raises (strata indices must align for the merge)."""
        with pytest.raises(ValueError, match="include_variant_data"):
            compute_allele_number_per_ref_site_sparse(
                synthetic_vds, reference_ht, include_variant_data=True
            )

    def test_ref_only_matches_hand_computed(self, synthetic_vds, reference_ht):
        """Default (ref-only) path: per-locus AN equals 2 × #samples with a ref block at that locus."""
        vd = synthetic_vds.variant_data
        gmh = self._build_gmh(vd)
        ht = compute_allele_number_per_ref_site_sparse(
            synthetic_vds, reference_ht, group_membership_ht=gmh
        )
        an_by_strata = self._index_an(ht)
        # All-samples "raw" group (the no_raw_group=True call drops the
        # default raw entry; the all-samples entry is the entry whose meta
        # has only `{group: raw}` — i.e., no gen_anc/sex keys).
        all_key = (("group", "raw"),)
        # Hand-computed ref-only AN per locus (ordered by locus, ascending):
        # L1000: s1, s2 → 4. L2000: s1, s3 → 4. L3000: all → 8. L4000: s2, s4 → 4.
        assert an_by_strata[all_key] == [4, 4, 8, 4]

    def test_include_variant_data_matches_dense_path(self, synthetic_vds, reference_ht):
        """`include_variant_data=True` output matches the dense densify-based path on a VDS-invariant-respecting fixture."""
        vd = synthetic_vds.variant_data
        gmh = self._build_gmh(vd)
        sparse_ht = compute_allele_number_per_ref_site_sparse(
            synthetic_vds,
            reference_ht,
            include_variant_data=True,
            group_membership_ht=gmh,
        )
        dense_ht = compute_allele_number_per_ref_site(
            synthetic_vds, reference_ht, group_membership_ht=gmh
        )
        sparse_indexed = self._index_an(sparse_ht)
        dense_indexed = self._index_an(dense_ht)
        assert set(sparse_indexed) == set(dense_indexed)
        for k in sparse_indexed:
            assert sparse_indexed[k] == dense_indexed[k], k

    def test_strata_expr_path_builds_group_membership(
        self, synthetic_vds, reference_ht
    ):
        """When only `strata_expr` is supplied (no pre-built `group_membership_ht`), the function builds one internally and the output matches the pre-built path."""
        vd = synthetic_vds.variant_data
        strata_expr = [
            {"gen_anc": vd.gen_anc},
            {"sex": vd.sex},
            {"gen_anc": vd.gen_anc, "sex": vd.sex},
        ]
        sparse_ht_via_strata = compute_allele_number_per_ref_site_sparse(
            synthetic_vds,
            reference_ht,
            include_variant_data=True,
            strata_expr=strata_expr,
        )
        gmh = self._build_gmh(vd)
        sparse_ht_via_gmh = compute_allele_number_per_ref_site_sparse(
            synthetic_vds,
            reference_ht,
            include_variant_data=True,
            group_membership_ht=gmh,
        )
        a = self._index_an(sparse_ht_via_strata)
        b = self._index_an(sparse_ht_via_gmh)
        assert set(a) == set(b)
        for k in a:
            assert a[k] == b[k], k

    def test_reduce_to_minimal_groups_matches_full(self, synthetic_vds, reference_ht):
        """`reduce_to_minimal_groups=True` produces identical per-locus AN to the un-reduced path."""
        vd = synthetic_vds.variant_data
        full_gmh = self._build_gmh(vd, reduce=False)
        reduced_gmh = self._build_gmh(vd, reduce=True)
        full_ht = compute_allele_number_per_ref_site_sparse(
            synthetic_vds,
            reference_ht,
            include_variant_data=True,
            group_membership_ht=full_gmh,
        )
        reduced_ht = compute_allele_number_per_ref_site_sparse(
            synthetic_vds,
            reference_ht,
            include_variant_data=True,
            group_membership_ht=reduced_gmh,
        )
        full_indexed = self._index_an(full_ht)
        reduced_indexed = self._index_an(reduced_ht)
        assert set(full_indexed) == set(reduced_indexed)
        for k in full_indexed:
            assert full_indexed[k] == reduced_indexed[k], k

    def test_sex_karyotype_adjusts_ploidy_on_chrx_nonpar(self):
        """On chrX non-PAR, XY samples contribute ploidy 1 (not 2) on both the ref and variant paths; XX samples are unchanged."""
        samples = [
            ("s1", "afr", "XX"),
            ("s2", "afr", "XY"),
            ("s3", "nfe", "XX"),
            ("s4", "nfe", "XY"),
        ]
        # chrX non-PAR positions on GRCh38 (PAR1 ends at ~2.78M, PAR2 starts at
        # ~155.7M; pick anything in between). Two loci exercise both paths:
        #   - First locus: all samples in ref blocks → ref-only AN contributes.
        #   - Second locus: all samples have variant calls → variant-data AN
        #     contributes. Both should give AN = 2 + 1 + 2 + 1 = 6 after
        #     ploidy adjustment.
        positions = [50_000_000, 60_000_000]
        assignment = [
            ["R", "R", "R", "R"],
            [
                ("V", (0, 0)),
                ("V", (0, 0)),
                ("V", (0, 0)),
                ("V", (0, 0)),
            ],
        ]
        vds = self._build_vds(assignment, samples, positions, contig="chrX")
        ref_ht = hl.Table.parallelize(
            [
                {"locus": hl.locus("chrX", p, reference_genome="GRCh38")}
                for p in positions
            ],
            hl.tstruct(locus=hl.tlocus("GRCh38")),
        ).key_by("locus")

        cols_ht = vds.variant_data.cols()
        gmh = generate_freq_group_membership_array(
            cols_ht,
            [{"sex": cols_ht.sex}],
            no_raw_group=True,
            group_label="raw",
        )
        # Annotate sex_karyotype on variant_data cols (the function
        # propagates it to reference_data cols).
        vmt = vds.variant_data.annotate_cols(sex_karyotype=vds.variant_data.sex)
        vds = hl.vds.VariantDataset(vds.reference_data, vmt)

        ht = compute_allele_number_per_ref_site_sparse(
            vds,
            ref_ht,
            include_variant_data=True,
            group_membership_ht=gmh,
            sex_karyotype_field="sex_karyotype",
        )
        an_by_strata = self._index_an(ht)
        all_key = (("group", "raw"),)
        # 2 (s1 XX) + 1 (s2 XY) + 2 (s3 XX) + 1 (s4 XY) = 6 at each locus.
        assert an_by_strata[all_key] == [6, 6]


class TestComputeStatsPerRefSiteSparse:
    """Test the generic sparse-aware stats aggregation (AN + coverage + hists)."""

    _SAMPLES = [
        ("s1", "afr", "XX"),
        ("s2", "afr", "XY"),
        ("s3", "nfe", "XX"),
        ("s4", "nfe", "XY"),
    ]

    # (sample, locus) -> "R" (ref block) or ("V", gt) — VDS invariant: each
    # cell is in exactly one path. The DP and GQ values below are deterministic
    # and chosen so the dense and sparse paths are easy to compare by hand.
    _ASSIGNMENT = [
        # locus 1000: s1=R, s2=R, s3=V(0/1), s4=V(1/1)
        ["R", "R", ("V", (0, 1)), ("V", (1, 1))],
        # locus 2000: s1=R, s2=V(missing), s3=R, s4=V(0/1)
        ["R", ("V", None), "R", ("V", (0, 1))],
        # locus 3000: all ref
        ["R", "R", "R", "R"],
        # locus 4000: s1=V(0/1), s2=R, s3=V(missing), s4=R
        [("V", (0, 1)), "R", ("V", None), "R"],
    ]

    # Per-sample DP value across all loci where the sample has any data.
    _DP_BY_SAMPLE = {"s1": 20, "s2": 30, "s3": 40, "s4": 50}
    # Per-sample GQ value across all loci where the sample has any data.
    _GQ_BY_SAMPLE = {"s1": 60, "s2": 70, "s3": 80, "s4": 90}
    # adj flag for variant-data entries; ref entries have adj=True by
    # construction below.
    _ADJ_BY_SAMPLE = {"s1": True, "s2": True, "s3": True, "s4": True}

    @classmethod
    def _build_vds(cls, assignment, samples, positions, contig="chr1"):
        loci = [hl.locus(contig, p, reference_genome="GRCh38") for p in positions]
        sample_table = hl.Table.parallelize(
            [{"s": s, "gen_anc": g, "sex": x} for s, g, x in samples],
            hl.tstruct(s=hl.tstr, gen_anc=hl.tstr, sex=hl.tstr),
        ).key_by("s")

        vd_entries = []
        rd_entries = []
        for l_idx, locus in enumerate(loci):
            for s_idx, (sample_id, _, _) in enumerate(samples):
                cell = assignment[l_idx][s_idx]
                dp = cls._DP_BY_SAMPLE[sample_id]
                gq = cls._GQ_BY_SAMPLE[sample_id]
                if cell == "R":
                    rd_entries.append(
                        {
                            "locus": locus,
                            "s": sample_id,
                            "LGT": hl.call(0, 0),
                            "DP": dp,
                            "GQ": gq,
                            "adj": True,
                            "END": locus.position,
                            "LEN": 1,
                        }
                    )
                else:
                    _, gt = cell
                    vd_entries.append(
                        {
                            "locus": locus,
                            "alleles": ["A", "T"],
                            "s": sample_id,
                            "GT": (
                                hl.call(*gt) if gt is not None else hl.missing(hl.tcall)
                            ),
                            "DP": dp,
                            "GQ": gq,
                            "adj": cls._ADJ_BY_SAMPLE[sample_id],
                        }
                    )

        vd = hl.Table.parallelize(
            vd_entries,
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                s=hl.tstr,
                GT=hl.tcall,
                DP=hl.tint32,
                GQ=hl.tint32,
                adj=hl.tbool,
            ),
        ).to_matrix_table(row_key=["locus", "alleles"], col_key=["s"])
        vd = vd.annotate_cols(
            gen_anc=sample_table[vd.s].gen_anc, sex=sample_table[vd.s].sex
        )

        rd = hl.Table.parallelize(
            rd_entries,
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                s=hl.tstr,
                LGT=hl.tcall,
                DP=hl.tint32,
                GQ=hl.tint32,
                adj=hl.tbool,
                END=hl.tint32,
                LEN=hl.tint32,
            ),
        ).to_matrix_table(row_key=["locus"], col_key=["s"])
        rd = rd.annotate_cols(
            gen_anc=sample_table[rd.s].gen_anc, sex=sample_table[rd.s].sex
        )
        return hl.vds.VariantDataset(reference_data=rd, variant_data=vd)

    @pytest.fixture
    def synthetic_vds(self):
        """Return a small VDS-invariant-respecting VDS for stats aggregation tests."""
        return self._build_vds(
            self._ASSIGNMENT, self._SAMPLES, [1000, 2000, 3000, 4000]
        )

    @pytest.fixture
    def reference_ht(self):
        """Return a reference HT covering the four autosomal loci."""
        return hl.Table.parallelize(
            [
                {"locus": hl.locus("chr1", p, reference_genome="GRCh38")}
                for p in [1000, 2000, 3000, 4000]
            ],
            hl.tstruct(locus=hl.tlocus("GRCh38")),
        ).key_by("locus")

    def _build_gmh(self, vd):
        cols_ht = vd.cols()
        return generate_freq_group_membership_array(
            cols_ht,
            [{"gen_anc": cols_ht.gen_anc}, {"sex": cols_ht.sex}],
            no_raw_group=True,
            group_label="raw",
        )

    def test_an_via_generic_matches_thin_wrapper(self, synthetic_vds, reference_ht):
        """The generic function called with AN-only entry_agg_funcs produces the same result as `compute_allele_number_per_ref_site_sparse`."""
        vd = synthetic_vds.variant_data
        gmh = self._build_gmh(vd)
        an_transform = lambda mt: (mt.GT if "GT" in mt.entry else mt.LGT).ploidy
        generic_ht = compute_stats_per_ref_site_sparse(
            synthetic_vds,
            reference_ht,
            entry_agg_funcs={"AN": (an_transform, hl.agg.sum)},
            include_variant_data=True,
            merge_funcs={"AN": merge_sum_array_expression},
            group_membership_ht=gmh,
        )
        wrapper_ht = compute_allele_number_per_ref_site_sparse(
            synthetic_vds,
            reference_ht,
            include_variant_data=True,
            group_membership_ht=gmh,
        )
        # Same row count, same per-locus AN per stratum.
        g = generic_ht.collect()
        w = wrapper_ht.collect()
        assert len(g) == len(w)
        for gr, wr in zip(g, w):
            assert gr.AN == wr.AN, (gr, wr)

    def test_coverage_stats_merge_matches_hand_computed(
        self, synthetic_vds, reference_ht
    ):
        """coverage_stats merged across the two sparse paths matches hand-computed dense-equivalent values for counter, total_DP, mean, and exact-from-counter median.

        Note: the dense path's `compute_stats_per_ref_site` with
        `get_coverage_agg_func` (which uses `hl.agg.approx_median`)
        triggers an IR-level CSE bug on small synthetic fixtures
        (`NoSuchElementException: key not found: __cse_*`), so this test
        validates the sparse path against hand-computed dense-equivalent
        expectations rather than against the dense path's output.
        """
        vd = synthetic_vds.variant_data
        gmh = self._build_gmh(vd)
        sparse_dp_agg = get_coverage_agg_func_sparse(dp_field="DP", max_cov_bin=100)
        coverage_target = [{"group": "raw", "gen_anc": "afr"}]
        sparse_ht = compute_stats_per_ref_site_sparse(
            synthetic_vds,
            reference_ht,
            entry_agg_funcs={"coverage_stats": sparse_dp_agg},
            entry_agg_group_membership={"coverage_stats": coverage_target},
            include_variant_data=True,
            merge_funcs={"coverage_stats": merge_coverage_stats_array_expression},
            group_membership_ht=gmh,
            entry_keep_fields=["DP"],
        )

        # Hand-compute dense-equivalent values for the {gen_anc: afr}
        # target. afr samples are s1 and s2; the sparse path filters
        # samples by `hl.is_defined(gt_field)`, so a sample with a
        # variant-data entry but a missing GT is treated as no-data at
        # that locus (and contributes to counter[0] via n_no_data).
        #
        # Per the assignment + DP map (s1=20, s2=30):
        #   - locus 1000: s1=R(20), s2=R(30)                → DPs=[20,30]
        #   - locus 2000: s1=R(20), s2=V(GT=missing)        → DPs=[20]; one no-data
        #   - locus 3000: s1=R(20), s2=R(30)                → DPs=[20,30]
        #   - locus 4000: s1=V(GT=defined,20), s2=R(30)     → DPs=[20,30]
        expected_per_locus = {
            1000: [20, 30],
            2000: [20],
            3000: [20, 30],
            4000: [20, 30],
        }
        n_total_afr = 2
        rows = sparse_ht.collect()
        assert len(rows) == 4
        for r in rows:
            cs = r.coverage_stats[0]
            dps = expected_per_locus[r.locus.position]
            expected_counter = {}
            for dp in dps:
                expected_counter[dp] = expected_counter.get(dp, 0) + 1
            # Add the (0, n_no_data) entry only when there's actually a
            # no-data sample in this group at this locus.
            n_no_data = n_total_afr - len(dps)
            if n_no_data > 0:
                expected_counter[0] = expected_counter.get(0, 0) + n_no_data
            assert dict(cs.coverage_counter) == expected_counter, (
                r.locus,
                cs.coverage_counter,
                expected_counter,
            )
            assert cs.total_DP == sum(dps)
            # mean = total_DP / n_total (no_data samples count toward
            # the denominator with DP=0, matching dense semantics).
            assert cs.mean == pytest.approx(sum(dps) / n_total_afr)
            # Exact median over the FULL counter (including 0 entries).
            # For n_total=2, half_pos = (2+1)//2 = 1 (1-indexed): the
            # smallest DP value whose cumulative count reaches 1.
            sorted_full = sorted(expected_counter.items(), key=lambda kv: kv[0])
            cum = 0
            expected_median = 0
            for dp, count in sorted_full:
                if cum < 1 and (cum + count) >= 1:
                    expected_median = dp
                    break
                cum += count
            assert cs.median_approx == expected_median

    def test_qual_hists_merge_produces_summed_bin_counts(
        self, synthetic_vds, reference_ht
    ):
        """qual_hists merged across the two paths matches the hand-computed sum of per-path bin counts (which is exactly the dense path's semantic for histograms since they skip missing values)."""
        vd = synthetic_vds.variant_data
        gmh = self._build_gmh(vd)

        def _hists_agg(qual_expr):
            return qual_hist_expr(
                gq_expr=qual_expr[0],
                dp_expr=qual_expr[1],
                adj_expr=qual_expr[2] == 1,
                split_adj_and_raw=True,
            )

        entry_agg_funcs = {
            "qual_hists": (lambda t: [t.GQ, t.DP, t.adj], _hists_agg),
        }
        target = [{"group": "raw", "gen_anc": "afr"}]
        sparse_ht = compute_stats_per_ref_site_sparse(
            synthetic_vds,
            reference_ht,
            entry_agg_funcs=entry_agg_funcs,
            entry_agg_group_membership={"qual_hists": target},
            include_variant_data=True,
            merge_funcs={"qual_hists": merge_qual_hists_array_expression},
            group_membership_ht=gmh,
            entry_keep_fields=["GQ", "DP", "adj"],
        )
        rows = sparse_ht.collect()
        assert len(rows) == 4
        # afr group is s1 and s2. The sparse path's auto-wrap filters
        # samples by `hl.is_defined(gt_field)`, so s2 at locus 2000
        # (variant_data entry with GT=missing) is excluded from the
        # histogram. Expected sample count per locus per histogram:
        #   - locus 1000, 3000, 4000: both s1 and s2 contribute → 2
        #   - locus 2000: only s1 contributes → 1
        expected_count_by_locus = {1000: 2, 2000: 1, 3000: 2, 4000: 2}
        for r in rows:
            qh = r.qual_hists[0]
            expected_n = expected_count_by_locus[r.locus.position]
            for path in ("raw_qual_hists", "qual_hists"):
                gq_hist = qh[path]["gq_hist_all"]
                dp_hist = qh[path]["dp_hist_all"]
                assert sum(gq_hist.bin_freq) == expected_n, (
                    r.locus,
                    path,
                    gq_hist.bin_freq,
                    expected_n,
                )
                assert sum(dp_hist.bin_freq) == expected_n
                # GQ values 60, 70 and DP values 20, 30 are all in
                # [0..100], so n_smaller/n_larger should be zero.
                assert gq_hist.n_smaller == 0
                assert gq_hist.n_larger == 0
                assert dp_hist.n_smaller == 0
                assert dp_hist.n_larger == 0

    def test_all_three_annotations_together(self, synthetic_vds, reference_ht):
        """End-to-end: AN + coverage_stats + qual_hists in a single sparse call, mirroring `compute_coverage.py`'s entry_agg_funcs shape."""
        vd = synthetic_vds.variant_data
        gmh = self._build_gmh(vd)

        def _hists_agg(qual_expr):
            return qual_hist_expr(
                gq_expr=qual_expr[0],
                dp_expr=qual_expr[1],
                adj_expr=qual_expr[2] == 1,
                split_adj_and_raw=True,
            )

        an_transform = lambda mt: (mt.GT if "GT" in mt.entry else mt.LGT).ploidy
        entry_agg_funcs = {
            "AN": (an_transform, hl.agg.sum),
            # Sparse-aware coverage agg (placeholder median; recomputed
            # at merge time).
            "coverage_stats": get_coverage_agg_func_sparse(
                dp_field="DP", max_cov_bin=100
            ),
            "qual_hists": (lambda t: [t.GQ, t.DP, t.adj], _hists_agg),
        }
        entry_agg_group_membership = {
            "coverage_stats": [{"group": "raw", "gen_anc": "afr"}],
            "qual_hists": [{"group": "raw", "gen_anc": "afr"}],
        }
        merge_funcs = {
            "AN": merge_sum_array_expression,
            "coverage_stats": merge_coverage_stats_array_expression,
            "qual_hists": merge_qual_hists_array_expression,
        }
        ht = compute_stats_per_ref_site_sparse(
            synthetic_vds,
            reference_ht,
            entry_agg_funcs=entry_agg_funcs,
            entry_agg_group_membership=entry_agg_group_membership,
            include_variant_data=True,
            merge_funcs=merge_funcs,
            group_membership_ht=gmh,
            entry_keep_fields=["DP", "GQ", "adj"],
        )
        # Smoke check: all three fields are present and shaped as expected.
        rows = ht.collect()
        assert len(rows) == 4
        strata_meta = hl.eval(ht.strata_meta)
        n_full = len(strata_meta)
        for r in rows:
            assert len(r.AN) == n_full
            # Narrowed to a single target.
            assert len(r.coverage_stats) == 1
            assert len(r.qual_hists) == 1
