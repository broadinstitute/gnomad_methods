"""Tests for the constraint utility module."""

import hail as hl
import pytest

from gnomad.utils.constraint import (
    build_constraint_consequence_groups,
    calculate_gerp_cutoffs,
    calculate_raw_z_score,
    compute_oe_upper_percentile_thresholds,
    count_observed_and_possible_by_group,
    counts_agg_expr,
    get_constraint_grouping_expr,
    oe_confidence_interval,
    rank_and_assign_bins,
    variant_observed_expr,
    weighted_sum_agg_expr,
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
    """Test the variant_observed_expr function."""

    def test_ac_positive_counts_as_one(self):
        """Test that a variant with AC > 0 counts as 1."""
        ht = hl.Table.parallelize(
            [{"freq": hl.Struct(AC=5, AF=0.01)}],
            hl.tstruct(freq=hl.tstruct(AC=hl.tint32, AF=hl.tfloat64)),
        )
        ht = ht.annotate(count=variant_observed_expr(freq_expr=ht.freq))
        result = ht.collect()[0]

        assert result.count == 1

    def test_ac_zero_counts_as_zero(self):
        """Test that a variant with AC == 0 counts as 0."""
        ht = hl.Table.parallelize(
            [{"freq": hl.Struct(AC=0, AF=0.0)}],
            hl.tstruct(freq=hl.tstruct(AC=hl.tint32, AF=hl.tfloat64)),
        )
        ht = ht.annotate(count=variant_observed_expr(freq_expr=ht.freq))
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
        ht = ht.annotate(count=variant_observed_expr(freq_expr=ht.freq, singleton=True))
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
        ht = ht.annotate(count=variant_observed_expr(freq_expr=ht.freq, max_af=0.01))
        results = ht.collect()

        assert results[0].count == 1
        assert results[1].count == 0

    def test_no_freq_counts_as_one(self):
        """Test that when no freq_expr is provided and no filtering, count is 1."""
        ht = hl.Table.parallelize(
            [{"x": 1}],
            hl.tstruct(x=hl.tint32),
        )
        ht = ht.annotate(count=variant_observed_expr(ht=ht))
        result = ht.collect()[0]

        assert result.count == 1

    def test_raises_when_no_ht_or_freq(self):
        """Test that ValueError is raised when neither ht nor freq_expr is given."""
        with pytest.raises(ValueError, match="Either `ht` or `freq_expr`"):
            variant_observed_expr()

    def test_max_af_zero_filters_all(self):
        """Test that max_af=0.0 filters out all variants (AF cannot be <= 0 with AC > 0)."""
        ht = hl.Table.parallelize(
            [
                {"freq": hl.Struct(AC=1, AF=0.001)},
                {"freq": hl.Struct(AC=5, AF=0.01)},
            ],
            hl.tstruct(freq=hl.tstruct(AC=hl.tint32, AF=hl.tfloat64)),
        )
        ht = ht.annotate(count=variant_observed_expr(freq_expr=ht.freq, max_af=0.0))
        results = ht.collect()

        assert results[0].count == 0
        assert results[1].count == 0

    def test_ht_fallback_to_freq_field(self):
        """Test that passing ht without freq_expr falls back to ht.freq."""
        ht = hl.Table.parallelize(
            [
                {"freq": hl.Struct(AC=1, AF=0.001)},
                {"freq": hl.Struct(AC=0, AF=0.0)},
            ],
            hl.tstruct(freq=hl.tstruct(AC=hl.tint32, AF=hl.tfloat64)),
        )
        ht = ht.annotate(count=variant_observed_expr(ht=ht, max_af=0.01))
        results = ht.collect()

        assert results[0].count == 1
        assert results[1].count == 0


class TestGetCountsAggExpr:
    """Test the counts_agg_expr function."""

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
        result = sample_table.aggregate(counts_agg_expr(freq_expr=sample_table.freq))

        # 3 variants have AC > 0.
        assert result.variant_count == 3

    def test_variant_count_with_max_af(self, sample_table):
        """Test variant count with max_af filter."""
        result = sample_table.aggregate(
            counts_agg_expr(freq_expr=sample_table.freq, max_af=0.01)
        )

        # AC=1/AF=0.001 and AC=5/AF=0.01 pass the filter.
        assert result.variant_count == 2

    def test_singleton_count(self, sample_table):
        """Test that singleton count is returned when requested."""
        result = sample_table.aggregate(
            counts_agg_expr(freq_expr=sample_table.freq, count_singletons=True)
        )

        assert result.singleton_count == 1
        assert result.variant_count == 3

    def test_no_singleton_key_by_default(self, sample_table):
        """Test that singleton_count is not present when not requested."""
        result = sample_table.aggregate(counts_agg_expr(freq_expr=sample_table.freq))

        assert not hasattr(result, "singleton_count")

    def test_raises_when_no_ht_or_freq(self):
        """Test that ValueError is raised when neither ht nor freq_expr is given."""
        with pytest.raises(ValueError, match="Either `ht` or `freq_expr`"):
            counts_agg_expr(freq_expr=None, ht=None)

    def test_ht_only_no_freq_expr(self, sample_table):
        """Test that passing ht without freq_expr falls back to ht.freq."""
        result = sample_table.aggregate(counts_agg_expr(ht=sample_table))

        # Same as test_variant_count_no_filter: 3 variants have AC > 0.
        assert result.variant_count == 3

    def test_max_af_zero_counts_none(self, sample_table):
        """Test that max_af=0.0 counts no variants."""
        result = sample_table.aggregate(
            counts_agg_expr(freq_expr=sample_table.freq, max_af=0.0)
        )

        assert result.variant_count == 0


class TestWeightedAggSumExpr:
    """Test the weighted_sum_agg_expr function."""

    def test_scalar_scalar(self):
        """Test weighted sum with two scalar expressions."""
        ht = hl.Table.parallelize(
            [{"val": 2.0, "weight": 3.0}, {"val": 4.0, "weight": 5.0}],
            hl.tstruct(val=hl.tfloat64, weight=hl.tfloat64),
        )
        result = ht.aggregate(weighted_sum_agg_expr(ht.val, ht.weight))

        # 2*3 + 4*5 = 26
        assert result == 26.0

    def test_array_array(self):
        """Test weighted sum with two array expressions (pairwise multiply)."""
        ht = hl.Table.parallelize(
            [
                {"val": [1.0, 2.0], "weight": [3.0, 4.0]},
                {"val": [5.0, 6.0], "weight": [7.0, 8.0]},
            ],
            hl.tstruct(val=hl.tarray(hl.tfloat64), weight=hl.tarray(hl.tfloat64)),
        )
        result = ht.aggregate(weighted_sum_agg_expr(ht.val, ht.weight))

        # element 0: 1*3 + 5*7 = 38, element 1: 2*4 + 6*8 = 56
        assert result == [38.0, 56.0]

    def test_scalar_array_mixed(self):
        """Test weighted sum with scalar expr and array weight (broadcast)."""
        ht = hl.Table.parallelize(
            [
                {"val": 2.0, "weight": [1.0, 10.0]},
                {"val": 3.0, "weight": [1.0, 10.0]},
            ],
            hl.tstruct(val=hl.tfloat64, weight=hl.tarray(hl.tfloat64)),
        )
        result = ht.aggregate(weighted_sum_agg_expr(ht.val, ht.weight))

        # element 0: 2*1 + 3*1 = 5, element 1: 2*10 + 3*10 = 50
        assert result == [5.0, 50.0]

    def test_array_scalar_mixed(self):
        """Test weighted sum with array expr and scalar weight (broadcast)."""
        ht = hl.Table.parallelize(
            [
                {"val": [1.0, 2.0], "weight": 3.0},
                {"val": [4.0, 5.0], "weight": 6.0},
            ],
            hl.tstruct(val=hl.tarray(hl.tfloat64), weight=hl.tfloat64),
        )
        result = ht.aggregate(weighted_sum_agg_expr(ht.val, ht.weight))

        # element 0: 1*3 + 4*6 = 27, element 1: 2*3 + 5*6 = 36
        assert result == [27.0, 36.0]

    def test_single_row(self):
        """Test weighted sum with a single row."""
        ht = hl.Table.parallelize(
            [{"val": 5.0, "weight": 2.0}],
            hl.tstruct(val=hl.tfloat64, weight=hl.tfloat64),
        )
        result = ht.aggregate(weighted_sum_agg_expr(ht.val, ht.weight))

        assert result == 10.0


class TestCountObservedAndPossibleByGroup:
    """Test the count_observed_and_possible_by_group function."""

    @pytest.fixture
    def sample_table(self):
        """Fixture to create a table with context, ref, alt, and count fields."""
        return hl.Table.parallelize(
            [
                {
                    "context": "ACG",
                    "ref": "C",
                    "alt": "T",
                    "methylation_level": 0,
                    "obs": [1, 0],
                    "poss": 1,
                },
                {
                    "context": "ACG",
                    "ref": "C",
                    "alt": "T",
                    "methylation_level": 0,
                    "obs": [1, 1],
                    "poss": 1,
                },
                {
                    "context": "TCG",
                    "ref": "C",
                    "alt": "A",
                    "methylation_level": 1,
                    "obs": [0, 0],
                    "poss": 1,
                },
            ],
            hl.tstruct(
                context=hl.tstr,
                ref=hl.tstr,
                alt=hl.tstr,
                methylation_level=hl.tint32,
                obs=hl.tarray(hl.tint32),
                poss=hl.tint32,
            ),
        )

    def test_basic_grouping(self, sample_table):
        """Test that rows are grouped by context, ref, alt, methylation_level."""
        result = count_observed_and_possible_by_group(
            sample_table,
            possible_expr=sample_table.poss,
            observed_expr=sample_table.obs,
        )
        rows = result.collect()

        # Two groups: (ACG, C, T, 0) and (TCG, C, A, 1).
        assert len(rows) == 2

    def test_observed_summed(self, sample_table):
        """Test that observed arrays are summed within groups."""
        result = count_observed_and_possible_by_group(
            sample_table,
            possible_expr=sample_table.poss,
            observed_expr=sample_table.obs,
        )
        rows = {r.context: r for r in result.collect()}

        # ACG group: [1,0] + [1,1] = [2,1]
        assert rows["ACG"].observed_variants == [2, 1]
        # TCG group: [0,0]
        assert rows["TCG"].observed_variants == [0, 0]

    def test_possible_summed(self, sample_table):
        """Test that possible counts are summed within groups."""
        result = count_observed_and_possible_by_group(
            sample_table,
            possible_expr=sample_table.poss,
            observed_expr=sample_table.obs,
        )
        rows = {r.context: r for r in result.collect()}

        assert rows["ACG"].possible_variants == 2
        assert rows["TCG"].possible_variants == 1

    def test_no_additional_grouping(self, sample_table):
        """Test grouping without methylation_level."""
        result = count_observed_and_possible_by_group(
            sample_table,
            possible_expr=sample_table.poss,
            observed_expr=sample_table.obs,
            additional_grouping=(),
        )
        rows = result.collect()

        # Without methylation_level, still two groups because context/ref/alt differ.
        assert len(rows) == 2

    def test_weight_exprs_dict(self, sample_table):
        """Test that weight_exprs produces a weighted sum field."""
        sample_table = sample_table.annotate(mu=0.5)
        result = count_observed_and_possible_by_group(
            sample_table,
            possible_expr=sample_table.poss,
            observed_expr=sample_table.obs,
            weight_exprs={"weighted_poss": sample_table.mu},
        )
        rows = {r.context: r for r in result.collect()}

        # ACG: poss=1*0.5 + poss=1*0.5 = 1.0
        assert abs(rows["ACG"].weighted_poss - 1.0) < 1e-6

    def test_additional_agg_sum_exprs(self, sample_table):
        """Test that additional_agg_sum_exprs sums extra fields."""
        sample_table = sample_table.annotate(extra=1)
        result = count_observed_and_possible_by_group(
            sample_table,
            possible_expr=sample_table.poss,
            observed_expr=sample_table.obs,
            additional_agg_sum_exprs={"extra": sample_table.extra},
        )
        rows = {r.context: r for r in result.collect()}

        assert rows["ACG"].extra == 2
        assert rows["TCG"].extra == 1


class TestGetConstraintGroupingExpr:
    """Test the get_constraint_grouping_expr function."""

    @pytest.fixture
    def vep_table(self):
        """Fixture to create a table with VEP annotation struct."""
        return hl.Table.parallelize(
            [
                {
                    "vep": hl.Struct(
                        most_severe_consequence="missense_variant",
                        lof=None,
                        polyphen_prediction="probably_damaging",
                        gene_symbol="BRCA1",
                        gene_id="ENSG00000012048",
                        transcript_id="ENST00000357654",
                        canonical=1,
                        mane_select="NM_007294.4",
                    ),
                },
                {
                    "vep": hl.Struct(
                        most_severe_consequence="stop_gained",
                        lof="HC",
                        polyphen_prediction=None,
                        gene_symbol="TP53",
                        gene_id="ENSG00000141510",
                        transcript_id="ENST00000269305",
                        canonical=0,
                        mane_select=None,
                    ),
                },
            ],
            hl.tstruct(
                vep=hl.tstruct(
                    most_severe_consequence=hl.tstr,
                    lof=hl.tstr,
                    polyphen_prediction=hl.tstr,
                    gene_symbol=hl.tstr,
                    gene_id=hl.tstr,
                    transcript_id=hl.tstr,
                    canonical=hl.tint32,
                    mane_select=hl.tstr,
                ),
            ),
        )

    def test_default_fields(self, vep_table):
        """Test that default groupings include annotation, modifier, gene, gene_id, transcript, canonical."""
        groupings = get_constraint_grouping_expr(vep_table.vep)

        assert set(groupings.keys()) == {
            "annotation",
            "modifier",
            "gene",
            "gene_id",
            "transcript",
            "canonical",
        }

    def test_modifier_uses_lof_over_polyphen(self, vep_table):
        """Test that modifier prefers lof when present, falls back to polyphen."""
        vep_table = vep_table.annotate(**get_constraint_grouping_expr(vep_table.vep))
        rows = vep_table.collect()

        # Row 0: lof is None, polyphen is "probably_damaging"
        assert rows[0].modifier == "probably_damaging"
        # Row 1: lof is "HC"
        assert rows[1].modifier == "HC"

    def test_modifier_falls_back_to_none_string(self):
        """Test that modifier is 'None' when both lof and polyphen are missing."""
        ht = hl.Table.parallelize(
            [
                {
                    "vep": hl.Struct(
                        most_severe_consequence="synonymous_variant",
                        lof=None,
                        polyphen_prediction=None,
                        gene_symbol="GENE1",
                        gene_id="ENSG00000000001",
                        transcript_id="ENST00000000001",
                        canonical=1,
                        mane_select=None,
                    ),
                },
            ],
            hl.tstruct(
                vep=hl.tstruct(
                    most_severe_consequence=hl.tstr,
                    lof=hl.tstr,
                    polyphen_prediction=hl.tstr,
                    gene_symbol=hl.tstr,
                    gene_id=hl.tstr,
                    transcript_id=hl.tstr,
                    canonical=hl.tint32,
                    mane_select=hl.tstr,
                ),
            ),
        )
        ht = ht.annotate(**get_constraint_grouping_expr(ht.vep))
        result = ht.collect()[0]

        assert result.modifier == "None"

    def test_canonical_is_boolean(self, vep_table):
        """Test that canonical is converted to a boolean."""
        vep_table = vep_table.annotate(**get_constraint_grouping_expr(vep_table.vep))
        rows = vep_table.collect()

        assert rows[0].canonical is True  # canonical=1
        assert rows[1].canonical is False  # canonical=0

    def test_include_mane_select(self, vep_table):
        """Test that mane_select is included when requested."""
        groupings = get_constraint_grouping_expr(
            vep_table.vep, include_mane_select_group=True
        )

        assert "mane_select" in groupings

        vep_table = vep_table.annotate(**groupings)
        rows = vep_table.collect()

        assert rows[0].mane_select is True  # has mane_select value
        assert rows[1].mane_select is False  # mane_select is None

    def test_exclude_transcript_and_canonical(self, vep_table):
        """Test that transcript and canonical can be excluded."""
        groupings = get_constraint_grouping_expr(
            vep_table.vep,
            include_transcript_group=False,
            include_canonical_group=False,
        )

        assert "transcript" not in groupings
        assert "canonical" not in groupings
        assert set(groupings.keys()) == {"annotation", "modifier", "gene", "gene_id"}

    def test_coverage_expr_included(self, vep_table):
        """Test that coverage is included when coverage_expr is provided."""
        vep_table = vep_table.annotate(cov=30)
        groupings = get_constraint_grouping_expr(
            vep_table.vep, coverage_expr=vep_table.cov
        )

        assert "coverage" in groupings

    def test_polyphen_missing_from_struct(self):
        """Test that missing polyphen_prediction field is handled gracefully."""
        ht = hl.Table.parallelize(
            [
                {
                    "vep": hl.Struct(
                        most_severe_consequence="stop_gained",
                        lof=None,
                        gene_symbol="GENE1",
                        gene_id="ENSG00000000001",
                        transcript_id="ENST00000000001",
                        canonical=1,
                        mane_select=None,
                    ),
                },
            ],
            hl.tstruct(
                vep=hl.tstruct(
                    most_severe_consequence=hl.tstr,
                    lof=hl.tstr,
                    gene_symbol=hl.tstr,
                    gene_id=hl.tstr,
                    transcript_id=hl.tstr,
                    canonical=hl.tint32,
                    mane_select=hl.tstr,
                ),
            ),
        )
        ht = ht.annotate(**get_constraint_grouping_expr(ht.vep))
        result = ht.collect()[0]

        # lof is None, polyphen is missing from struct → modifier should be "None"
        assert result.modifier == "None"
