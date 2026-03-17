"""Tests for the constraint utility module."""

import hail as hl
import pytest

from gnomad.utils.constraint import (
    build_constraint_consequence_groups,
    calculate_gerp_cutoffs,
    calculate_raw_z_score,
    compute_oe_upper_percentile_thresholds,
    get_counts_agg_expr,
    oe_confidence_interval,
    rank_and_assign_bins,
    single_variant_count_expr,
)


class TestOeConfidenceInterval:
    """Test the oe_confidence_interval function."""

    def test_gamma_returns_lower_and_upper(self):
        """Test that gamma method returns a struct with lower < upper."""
        ht = hl.Table.parallelize(
            [{"obs": 10, "exp": 20.0}],
            hl.tstruct(obs=hl.tint64, exp=hl.tfloat64),
        )
        ht = ht.annotate(ci=oe_confidence_interval(ht.obs, ht.exp, method="gamma"))
        result = ht.collect()[0]

        assert result.ci.lower < result.ci.upper

    def test_poisson_returns_lower_and_upper(self):
        """Test that poisson method returns a struct with lower < upper."""
        ht = hl.Table.parallelize(
            [{"obs": 10, "exp": 20.0}],
            hl.tstruct(obs=hl.tint64, exp=hl.tfloat64),
        )
        ht = ht.annotate(ci=oe_confidence_interval(ht.obs, ht.exp, method="poisson"))
        result = ht.collect()[0]

        assert result.ci.lower < result.ci.upper

    def test_poisson_obs_zero_lower_is_zero(self):
        """Test that poisson method returns lower=0 when obs=0."""
        ht = hl.Table.parallelize(
            [{"obs": 0, "exp": 10.0}],
            hl.tstruct(obs=hl.tint64, exp=hl.tfloat64),
        )
        ht = ht.annotate(ci=oe_confidence_interval(ht.obs, ht.exp, method="poisson"))
        result = ht.collect()[0]

        assert result.ci.lower == 0

    def test_invalid_method_raises(self):
        """Test that an invalid method raises ValueError."""
        ht = hl.Table.parallelize(
            [{"obs": 5, "exp": 10.0}],
            hl.tstruct(obs=hl.tint64, exp=hl.tfloat64),
        )
        with pytest.raises(ValueError, match="Unknown CI method"):
            ht.annotate(ci=oe_confidence_interval(ht.obs, ht.exp, method="invalid"))

    def test_gamma_and_poisson_give_similar_results(self):
        """Test that gamma and poisson methods give roughly similar CIs."""
        ht = hl.Table.parallelize(
            [{"obs": 15, "exp": 20.0}],
            hl.tstruct(obs=hl.tint64, exp=hl.tfloat64),
        )
        ht = ht.annotate(
            ci_gamma=oe_confidence_interval(ht.obs, ht.exp, method="gamma"),
            ci_poisson=oe_confidence_interval(ht.obs, ht.exp, method="poisson"),
        )
        result = ht.collect()[0]

        # Both should be in a similar ballpark (within 0.2 of each other).
        assert abs(result.ci_gamma.lower - result.ci_poisson.lower) < 0.2
        assert abs(result.ci_gamma.upper - result.ci_poisson.upper) < 0.2


class TestCalculateRawZScore:
    """Test the calculate_raw_z_score function."""

    def test_obs_less_than_exp_positive_z(self):
        """Test that fewer observed than expected gives a positive z-score."""
        ht = hl.Table.parallelize(
            [{"obs": 5, "exp": 20.0}],
            hl.tstruct(obs=hl.tint64, exp=hl.tfloat64),
        )
        ht = ht.annotate(z=calculate_raw_z_score(ht.obs, ht.exp))
        result = ht.collect()[0]

        assert result.z > 0

    def test_obs_greater_than_exp_negative_z(self):
        """Test that more observed than expected gives a negative z-score."""
        ht = hl.Table.parallelize(
            [{"obs": 30, "exp": 10.0}],
            hl.tstruct(obs=hl.tint64, exp=hl.tfloat64),
        )
        ht = ht.annotate(z=calculate_raw_z_score(ht.obs, ht.exp))
        result = ht.collect()[0]

        assert result.z < 0

    def test_obs_equals_exp_zero_z(self):
        """Test that obs == exp gives z-score of 0."""
        ht = hl.Table.parallelize(
            [{"obs": 10, "exp": 10.0}],
            hl.tstruct(obs=hl.tint64, exp=hl.tfloat64),
        )
        ht = ht.annotate(z=calculate_raw_z_score(ht.obs, ht.exp))
        result = ht.collect()[0]

        assert result.z == 0.0


class TestCalculateGerpCutoffs:
    """Test the calculate_gerp_cutoffs function."""

    def test_cutoffs_within_expected_range(self):
        """Test that GERP cutoffs fall within the data range."""
        ht = hl.Table.parallelize(
            [{"gerp": float(i)} for i in range(100)],
            hl.tstruct(gerp=hl.tfloat64),
        )
        lower, upper = calculate_gerp_cutoffs(ht)

        # 5th percentile should be around 5, 95th around 95.
        assert 0 <= lower <= 10
        assert 90 <= upper <= 99

    def test_custom_percentiles(self):
        """Test with custom percentile thresholds."""
        ht = hl.Table.parallelize(
            [{"gerp": float(i)} for i in range(100)],
            hl.tstruct(gerp=hl.tfloat64),
        )
        lower, upper = calculate_gerp_cutoffs(
            ht, lower_percentile=0.25, upper_percentile=0.75
        )

        assert 20 <= lower <= 30
        assert 70 <= upper <= 80

    def test_custom_gerp_expr(self):
        """Test with a custom GERP expression."""
        ht = hl.Table.parallelize(
            [{"my_gerp": float(i)} for i in range(100)],
            hl.tstruct(my_gerp=hl.tfloat64),
        )
        lower, upper = calculate_gerp_cutoffs(ht, gerp_expr=ht.my_gerp)

        assert 0 <= lower <= 10
        assert 90 <= upper <= 99


class TestBuildConstraintConsequenceGroups:
    """Test the build_constraint_consequence_groups function."""

    @pytest.fixture
    def sample_table(self):
        """Fixture to create a table with consequence and modifier fields."""
        return hl.Table.parallelize(
            [
                {"csq": "synonymous_variant", "modifier": "NA"},
                {"csq": "missense_variant", "modifier": "NA"},
                {"csq": "stop_gained", "modifier": "HC"},
                {"csq": "splice_donor_variant", "modifier": "LC"},
                {"csq": "missense_variant", "modifier": "NA"},
            ],
            hl.tstruct(csq=hl.tstr, modifier=hl.tstr),
        )

    def test_meta_contains_expected_keys(self, sample_table):
        """Test that returned meta dicts contain expected keys."""
        _, meta = build_constraint_consequence_groups(
            sample_table.csq, sample_table.modifier
        )

        meta_keys = [set(m.keys()) for m in meta]
        assert {"csq_set"} in meta_keys
        assert {"lof"} in meta_keys

    def test_default_groups_count(self, sample_table):
        """Test that the default number of groups is 5 (syn, mis, classic, hc_lc, hc)."""
        filters, meta = build_constraint_consequence_groups(
            sample_table.csq, sample_table.modifier
        )

        assert len(meta) == 5
        assert len(filters) == 5

    def test_syn_filter_selects_correct_rows(self, sample_table):
        """Test that the synonymous filter selects only synonymous variants."""
        filters, meta = build_constraint_consequence_groups(
            sample_table.csq, sample_table.modifier
        )

        # Find the syn filter.
        syn_idx = next(i for i, m in enumerate(meta) if m.get("csq_set") == "syn")
        ht = sample_table.annotate(is_syn=filters[syn_idx])
        results = ht.collect()

        assert results[0].is_syn is True  # synonymous_variant
        assert results[1].is_syn is False  # missense_variant
        assert results[2].is_syn is False  # stop_gained

    def test_hc_lof_filter_selects_correct_rows(self, sample_table):
        """Test that the HC LoF filter selects only HC modifier rows."""
        filters, meta = build_constraint_consequence_groups(
            sample_table.csq, sample_table.modifier
        )

        hc_idx = next(i for i, m in enumerate(meta) if m.get("lof") == "hc")
        ht = sample_table.annotate(is_hc=filters[hc_idx])
        results = ht.collect()

        assert results[2].is_hc is True  # stop_gained, HC
        assert results[3].is_hc is False  # splice_donor_variant, LC

    def test_meta_values(self, sample_table):
        """Test that meta contains the expected value combinations."""
        _, meta = build_constraint_consequence_groups(
            sample_table.csq, sample_table.modifier
        )

        meta_values = [tuple(sorted(m.items())) for m in meta]
        assert (("csq_set", "syn"),) in meta_values
        assert (("csq_set", "mis"),) in meta_values
        assert (("lof", "classic"),) in meta_values
        assert (("lof", "hc_lc"),) in meta_values
        assert (("lof", "hc"),) in meta_values


class TestRankAndAssignBins:
    """Test the rank_and_assign_bins function."""

    @pytest.fixture
    def sample_table(self):
        """Fixture to create a keyed table with float values."""
        return hl.Table.parallelize(
            [{"id": i, "val": float(i)} for i in range(10)],
            hl.tstruct(id=hl.tint32, val=hl.tfloat64),
            key="id",
        )

    def test_rank_ordering_ascending(self, sample_table):
        """Test that ranks are assigned in ascending order of value."""
        result = sample_table.annotate(bins=rank_and_assign_bins(sample_table.val))
        rows = result.order_by("id").collect()

        # Values 0..9 should get ranks 0..9.
        for i, r in enumerate(rows):
            assert r.bins.rank == i

    def test_default_bin_fields(self, sample_table):
        """Test that default bin granularities produce expected fields."""
        result = sample_table.annotate(bins=rank_and_assign_bins(sample_table.val))
        row = result.collect()[0]

        assert hasattr(row.bins, "bin_percentile")
        assert hasattr(row.bins, "bin_decile")
        assert hasattr(row.bins, "bin_sextile")

    def test_custom_bin_granularities(self, sample_table):
        """Test with custom bin granularities."""
        result = sample_table.annotate(
            bins=rank_and_assign_bins(
                sample_table.val, bin_granularities={"quintile": 5}
            )
        )
        rows = result.order_by("id").collect()

        # 10 rows, quintile bins: 0,0, 1,1, 2,2, 3,3, 4,4
        bins = [r.bins.bin_quintile for r in rows]
        assert bins == [0, 0, 1, 1, 2, 2, 3, 3, 4, 4]

    def test_decile_bins_range(self, sample_table):
        """Test that decile bins are in the range [0, 9]."""
        result = sample_table.annotate(bins=rank_and_assign_bins(sample_table.val))
        rows = result.collect()
        deciles = [r.bins.bin_decile for r in rows]

        assert min(deciles) == 0
        assert max(deciles) == 9

    def test_single_row(self):
        """Test ranking a single-row table."""
        ht = hl.Table.parallelize(
            [{"id": 1, "val": 5.0}],
            hl.tstruct(id=hl.tint32, val=hl.tfloat64),
            key="id",
        )
        result = ht.annotate(bins=rank_and_assign_bins(ht.val))
        row = result.collect()[0]

        assert row.bins.rank == 0
        assert row.bins.bin_decile == 0
        assert row.bins.bin_percentile == 0

    def test_tied_values(self):
        """Test that tied values all receive the same bin."""
        ht = hl.Table.parallelize(
            [{"id": i, "val": 1.0} for i in range(10)],
            hl.tstruct(id=hl.tint32, val=hl.tfloat64),
            key="id",
        )
        result = ht.annotate(
            bins=rank_and_assign_bins(ht.val, bin_granularities={"half": 2})
        )
        rows = result.collect()
        bins = {r.bins.bin_half for r in rows}

        # All same value — ranks are arbitrary but bins should be contiguous.
        assert bins.issubset({0, 1})

    def test_descending_input_still_ranks_ascending(self):
        """Test that input order does not affect ascending rank assignment."""
        ht = hl.Table.parallelize(
            [{"id": i, "val": float(9 - i)} for i in range(10)],
            hl.tstruct(id=hl.tint32, val=hl.tfloat64),
            key="id",
        )
        result = ht.annotate(bins=rank_and_assign_bins(ht.val))
        rows = result.order_by("id").collect()

        # id=0 has val=9 (largest) -> rank 9; id=9 has val=0 (smallest) -> rank 0.
        assert rows[0].bins.rank == 9
        assert rows[9].bins.rank == 0

    def test_negative_values(self):
        """Test ranking with negative values."""
        ht = hl.Table.parallelize(
            [{"id": i, "val": float(i - 5)} for i in range(10)],
            hl.tstruct(id=hl.tint32, val=hl.tfloat64),
            key="id",
        )
        result = ht.annotate(bins=rank_and_assign_bins(ht.val))
        rows = result.order_by("id").collect()

        # val ranges from -5 to 4; id=0 (val=-5) should get rank 0.
        assert rows[0].bins.rank == 0
        assert rows[9].bins.rank == 9


class TestComputeOeUpperPercentileThresholds:
    """Test the compute_oe_upper_percentile_thresholds function."""

    @pytest.fixture
    def sample_table(self):
        """Fixture to create a table with known metric values."""
        return hl.Table.parallelize(
            [{"metric": float(i) / 100.0} for i in range(101)],
            hl.tstruct(metric=hl.tfloat64),
        )

    def test_returns_dict_with_correct_keys(self, sample_table):
        """Test that the returned dict has keys matching requested percentiles."""
        result = compute_oe_upper_percentile_thresholds(
            sample_table, sample_table.metric, percentiles=(10, 50, 90)
        )

        assert set(result.keys()) == {10, 50, 90}

    def test_thresholds_in_ascending_order(self, sample_table):
        """Test that thresholds are in ascending order for ascending percentiles."""
        result = compute_oe_upper_percentile_thresholds(
            sample_table, sample_table.metric, percentiles=(10, 25, 50, 75)
        )

        values = [result[p] for p in (10, 25, 50, 75)]
        assert values == sorted(values)

    def test_thresholds_within_data_range(self, sample_table):
        """Test that thresholds fall within the data range."""
        result = compute_oe_upper_percentile_thresholds(
            sample_table, sample_table.metric, percentiles=(10, 50, 90)
        )

        for v in result.values():
            assert 0.0 <= v <= 1.0

    def test_outlier_filtering(self):
        """Test that outlier rows are excluded from threshold computation."""
        ht = hl.Table.parallelize(
            [{"metric": float(i) / 10.0, "is_outlier": i > 8} for i in range(10)],
            hl.tstruct(metric=hl.tfloat64, is_outlier=hl.tbool),
        )

        result = compute_oe_upper_percentile_thresholds(
            ht, ht.metric, outlier_expr=ht.is_outlier, percentiles=(50,)
        )

        # With outliers excluded, median should be around 0.4 (values 0..8 / 10).
        assert result[50] < 0.5

    def test_transcript_filter(self):
        """Test that transcript_filter_expr restricts which rows are used."""
        ht = hl.Table.parallelize(
            [
                {"metric": 0.1, "include": True},
                {"metric": 0.2, "include": True},
                {"metric": 0.9, "include": False},
                {"metric": 1.0, "include": False},
            ],
            hl.tstruct(metric=hl.tfloat64, include=hl.tbool),
        )

        result = compute_oe_upper_percentile_thresholds(
            ht,
            ht.metric,
            transcript_filter_expr=ht.include,
            percentiles=(50,),
        )

        # Only 0.1 and 0.2 included, so median should be ~0.15.
        assert result[50] < 0.3

    def test_default_percentiles(self):
        """Test that default percentiles (1,5,10,15,25,50,75) are used."""
        ht = hl.Table.parallelize(
            [{"metric": float(i) / 100.0} for i in range(101)],
            hl.tstruct(metric=hl.tfloat64),
        )

        result = compute_oe_upper_percentile_thresholds(ht, ht.metric)
        assert set(result.keys()) == {1, 5, 10, 15, 25, 50, 75}

    def test_missing_metric_excluded(self):
        """Test that rows with missing metric values are excluded."""
        ht = hl.Table.parallelize(
            [
                {"metric": 0.1},
                {"metric": 0.2},
                {"metric": None},
                {"metric": None},
            ],
            hl.tstruct(metric=hl.tfloat64),
        )

        result = compute_oe_upper_percentile_thresholds(
            ht, ht.metric, percentiles=(50,)
        )

        # Only 0.1 and 0.2 are non-missing, median ~ 0.15.
        assert result[50] < 0.3

    def test_no_filters_includes_all_rows(self):
        """Test that omitting outlier and transcript filters includes all rows."""
        ht = hl.Table.parallelize(
            [{"metric": float(i)} for i in range(10)],
            hl.tstruct(metric=hl.tfloat64),
        )

        result_all = compute_oe_upper_percentile_thresholds(
            ht, ht.metric, percentiles=(50,)
        )
        result_filtered = compute_oe_upper_percentile_thresholds(
            ht,
            ht.metric,
            outlier_expr=hl.literal(False),
            transcript_filter_expr=hl.literal(True),
            percentiles=(50,),
        )

        # Should be the same when filters are effectively no-ops.
        assert abs(result_all[50] - result_filtered[50]) < 0.01

    def test_combined_outlier_and_transcript_filter(self):
        """Test that outlier and transcript filters are combined correctly."""
        ht = hl.Table.parallelize(
            [
                {"metric": 0.1, "include": True, "is_outlier": False},
                {"metric": 0.5, "include": True, "is_outlier": True},
                {"metric": 0.9, "include": False, "is_outlier": False},
            ],
            hl.tstruct(metric=hl.tfloat64, include=hl.tbool, is_outlier=hl.tbool),
        )

        result = compute_oe_upper_percentile_thresholds(
            ht,
            ht.metric,
            outlier_expr=ht.is_outlier,
            transcript_filter_expr=ht.include,
            percentiles=(50,),
        )

        # Only first row passes both filters, so median is ~0.1.
        assert abs(result[50] - 0.1) < 0.05


class TestRankVsThresholdBinning:
    """Compare rank_and_assign_bins and compute_oe_upper_percentile_thresholds.

    rank_and_assign_bins assigns bins by exact row rank (deterministic for
    distinct values, arbitrary tie-breaking for duplicates).
    compute_oe_upper_percentile_thresholds computes approximate quantile
    thresholds — binning by value comparison puts all tied rows in the same bin.

    These should agree when all values are distinct but can diverge when many
    values are tied.
    """

    def test_distinct_values_bins_agree(self):
        """Test that both methods agree when all values are distinct."""
        n = 100
        ht = hl.Table.parallelize(
            [{"id": i, "val": float(i) / n} for i in range(n)],
            hl.tstruct(id=hl.tint32, val=hl.tfloat64),
            key="id",
        )

        # Rank-based decile bins.
        ranked = ht.annotate(
            bins=rank_and_assign_bins(ht.val, bin_granularities={"decile": 10})
        )
        rank_rows = ranked.order_by("id").collect()
        rank_bins = {r.id: r.bins.bin_decile for r in rank_rows}

        # Threshold-based decile bins (boundaries at 10th, 20th, ..., 90th).
        thresholds = compute_oe_upper_percentile_thresholds(
            ht, ht.val, percentiles=tuple(range(10, 100, 10))
        )
        threshold_list = [thresholds[p] for p in sorted(thresholds)]

        rows = ht.order_by("id").collect()
        for r in rows:
            rank_bin = rank_bins[r.id]
            # Threshold bin: count how many thresholds the value exceeds.
            thresh_bin = sum(1 for t in threshold_list if r.val >= t)
            assert (
                rank_bin == thresh_bin
            ), f"id={r.id}, val={r.val}: rank_bin={rank_bin}, thresh_bin={thresh_bin}"

    def test_tied_values_threshold_bins_consistent(self):
        """Test that threshold-based bins put all tied values in the same bin.

        rank_and_assign_bins may split tied values across bins because it
        assigns unique ranks to each row regardless of ties.
        """
        # 50 rows at val=0.1 and 50 rows at val=0.9.
        ht = hl.Table.parallelize(
            [{"id": i, "val": 0.1 if i < 50 else 0.9} for i in range(100)],
            hl.tstruct(id=hl.tint32, val=hl.tfloat64),
            key="id",
        )

        # Threshold-based: all 0.1 rows must be in the same bin, all 0.9 in another.
        thresholds = compute_oe_upper_percentile_thresholds(
            ht, ht.val, percentiles=(50,)
        )
        rows = ht.collect()
        low_bins = set()
        high_bins = set()
        for r in rows:
            bin_val = sum(1 for t in thresholds.values() if r.val >= t)
            if r.val < 0.5:
                low_bins.add(bin_val)
            else:
                high_bins.add(bin_val)

        # All 0.1 rows get the same bin, all 0.9 rows get the same bin.
        assert len(low_bins) == 1
        assert len(high_bins) == 1
        assert low_bins != high_bins

    def test_tied_values_rank_bins_may_split(self):
        """Test that rank_and_assign_bins splits tied values across bins.

        With 100 identical values and 10 decile bins, rank-based binning will
        arbitrarily assign ~10 rows per bin even though the values are the same.
        """
        ht = hl.Table.parallelize(
            [{"id": i, "val": 1.0} for i in range(100)],
            hl.tstruct(id=hl.tint32, val=hl.tfloat64),
            key="id",
        )

        ranked = ht.annotate(
            bins=rank_and_assign_bins(ht.val, bin_granularities={"decile": 10})
        )
        rows = ranked.collect()
        deciles = {r.bins.bin_decile for r in rows}

        # All values are identical, but rank bins still spread across deciles.
        assert len(deciles) > 1

    def test_mostly_tied_with_outlier_diverges(self):
        """Test divergence when most values are tied with one outlier.

        Threshold-based: outlier alone in top bin, all tied values in bin 0.
        Rank-based: tied values spread across multiple bins, outlier in top bin.
        """
        ht = hl.Table.parallelize(
            [{"id": i, "val": 0.5 if i < 99 else 10.0} for i in range(100)],
            hl.tstruct(id=hl.tint32, val=hl.tfloat64),
            key="id",
        )

        # Threshold at the 50th percentile.
        thresholds = compute_oe_upper_percentile_thresholds(
            ht, ht.val, percentiles=(50,)
        )
        # The 50th pct threshold should be ~0.5 (the tied value), so all 0.5
        # rows fall at or above it → threshold bin 1, and the outlier also bin 1.
        # But with rank-based bins, about half the 0.5 rows are in bin 0 and
        # half in bin 1.
        ranked = ht.annotate(
            bins=rank_and_assign_bins(ht.val, bin_granularities={"half": 2})
        )
        rank_rows = ranked.collect()

        # Rank-based: the 0.5-valued rows are split across bin 0 and bin 1.
        tied_rank_bins = {r.bins.bin_half for r in rank_rows if r.val == 0.5}
        assert len(tied_rank_bins) == 2, "Rank-based should split tied values"

        # Threshold-based: all 0.5-valued rows are in the same bin.
        tied_thresh_bins = set()
        for r in rank_rows:
            if r.val == 0.5:
                tied_thresh_bins.add(sum(1 for t in thresholds.values() if r.val >= t))
        assert len(tied_thresh_bins) == 1, "Threshold-based should not split ties"


class TestSingleVariantCountExpr:
    """Test the single_variant_count_expr function."""

    def test_ac_positive_counts_as_one(self):
        """Test that a variant with AC > 0 counts as 1."""
        ht = hl.Table.parallelize(
            [{"freq": hl.Struct(AC=5, AF=0.01)}],
            hl.tstruct(freq=hl.tstruct(AC=hl.tint32, AF=hl.tfloat64)),
        )
        ht = ht.annotate(count=single_variant_count_expr(freq_expr=ht.freq))
        result = ht.collect()[0]

        assert result.count == 1

    def test_ac_zero_counts_as_zero(self):
        """Test that a variant with AC == 0 counts as 0."""
        ht = hl.Table.parallelize(
            [{"freq": hl.Struct(AC=0, AF=0.0)}],
            hl.tstruct(freq=hl.tstruct(AC=hl.tint32, AF=hl.tfloat64)),
        )
        ht = ht.annotate(count=single_variant_count_expr(freq_expr=ht.freq))
        result = ht.collect()[0]

        assert result.count == 0

    def test_singleton_counts_only_ac_one(self):
        """Test that singleton mode only counts variants with AC == 1."""
        ht = hl.Table.parallelize(
            [
                {"freq": hl.Struct(AC=1, AF=0.001)},
                {"freq": hl.Struct(AC=5, AF=0.01)},
            ],
            hl.tstruct(freq=hl.tstruct(AC=hl.tint32, AF=hl.tfloat64)),
        )
        ht = ht.annotate(
            count=single_variant_count_expr(freq_expr=ht.freq, singleton=True)
        )
        results = ht.collect()

        assert results[0].count == 1
        assert results[1].count == 0

    def test_max_af_filters_by_frequency(self):
        """Test that max_af filters variants by allele frequency."""
        ht = hl.Table.parallelize(
            [
                {"freq": hl.Struct(AC=5, AF=0.001)},
                {"freq": hl.Struct(AC=5, AF=0.1)},
            ],
            hl.tstruct(freq=hl.tstruct(AC=hl.tint32, AF=hl.tfloat64)),
        )
        ht = ht.annotate(
            count=single_variant_count_expr(freq_expr=ht.freq, max_af=0.01)
        )
        results = ht.collect()

        assert results[0].count == 1
        assert results[1].count == 0

    def test_no_freq_counts_as_one(self):
        """Test that when no freq_expr is provided and no filtering, count is 1."""
        ht = hl.Table.parallelize(
            [{"x": 1}],
            hl.tstruct(x=hl.tint32),
        )
        ht = ht.annotate(count=single_variant_count_expr(ht=ht))
        result = ht.collect()[0]

        assert result.count == 1

    def test_raises_when_no_ht_or_freq(self):
        """Test that ValueError is raised when neither ht nor freq_expr is given."""
        with pytest.raises(ValueError, match="Either `ht` or `freq_expr`"):
            single_variant_count_expr()


class TestGetCountsAggExpr:
    """Test the get_counts_agg_expr function."""

    @pytest.fixture
    def sample_table(self):
        """Fixture to create a table with frequency data."""
        return hl.Table.parallelize(
            [
                {"freq": hl.Struct(AC=1, AF=0.001)},
                {"freq": hl.Struct(AC=5, AF=0.01)},
                {"freq": hl.Struct(AC=0, AF=0.0)},
                {"freq": hl.Struct(AC=3, AF=0.05)},
            ],
            hl.tstruct(freq=hl.tstruct(AC=hl.tint32, AF=hl.tfloat64)),
        )

    def test_variant_count_no_filter(self, sample_table):
        """Test variant count with no filtering (AC > 0)."""
        result = sample_table.aggregate(
            get_counts_agg_expr(freq_expr=sample_table.freq)
        )

        # 3 variants have AC > 0.
        assert result.variant_count == 3

    def test_variant_count_with_max_af(self, sample_table):
        """Test variant count with max_af filter."""
        result = sample_table.aggregate(
            get_counts_agg_expr(freq_expr=sample_table.freq, max_af=0.01)
        )

        # AC=1/AF=0.001 and AC=5/AF=0.01 pass the filter.
        assert result.variant_count == 2

    def test_singleton_count(self, sample_table):
        """Test that singleton count is returned when requested."""
        result = sample_table.aggregate(
            get_counts_agg_expr(freq_expr=sample_table.freq, count_singletons=True)
        )

        assert result.singleton_count == 1
        assert result.variant_count == 3

    def test_no_singleton_key_by_default(self, sample_table):
        """Test that singleton_count is not present when not requested."""
        result = sample_table.aggregate(
            get_counts_agg_expr(freq_expr=sample_table.freq)
        )

        assert not hasattr(result, "singleton_count")

    def test_raises_when_no_ht_or_freq(self):
        """Test that ValueError is raised when neither ht nor freq_expr is given."""
        with pytest.raises(ValueError, match="Either `ht` or `freq_expr`"):
            get_counts_agg_expr(freq_expr=None, ht=None)
