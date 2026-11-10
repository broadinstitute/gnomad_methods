"""Tests for the regional mutation rate utility module."""

import logging

import hail as hl
import numpy as np
import pandas as pd
import pytest

from gnomad.utils.constraint import (
    choose_regional_mutation_rate_penalty,
    compute_regional_mutation_rate_ht,
    fit_regional_mutation_rate_model,
    fit_regional_mutation_rate_models,
    normalize_regional_mutation_rate_expr,
    predict_regional_mutation_rate_probability,
    regional_mutation_rate_expr,
    regional_mutation_rate_models_to_ht,
    regional_mutation_rate_normalization_factors_agg_expr,
    tune_regional_mutation_rate_penalty,
)

INTERVAL_TYPE = hl.tinterval(hl.tlocus("GRCh38"))
"""Type of the window intervals in the test feature and rate Tables."""


def _expit(x: np.ndarray) -> np.ndarray:
    """Return the logistic function of ``x``."""
    return 1.0 / (1.0 + np.exp(-x))


def _window(start: int, end: int) -> hl.Interval:
    """Return the chr1 interval ``[start, end)``."""
    return hl.Interval(
        hl.Locus("chr1", start, "GRCh38"), hl.Locus("chr1", end, "GRCh38")
    )


def _training_df(n: int = 4000, seed: int = 0) -> pd.DataFrame:
    """Return training rows whose labels follow a logistic model of ``a`` and ``b``."""
    rng = np.random.default_rng(seed)
    a = rng.normal(size=n)
    b = 3 * rng.normal(size=n) + 5
    y = rng.random(n) < _expit(-2 + 0.8 * a - 0.5 * (b - 5) / 3)
    return pd.DataFrame({"y": y.astype("int64"), "a": a, "b": b})


def _features_ht(rows: list) -> hl.Table:
    """Return a feature Table keyed by window with float features ``a`` and ``b``."""
    return hl.Table.parallelize(
        [{"interval": _window(s, e), "a": a, "b": b} for s, e, a, b in rows],
        hl.tstruct(interval=INTERVAL_TYPE, a=hl.tfloat64, b=hl.tfloat64),
        key="interval",
    )


MODELS = [
    hl.Struct(
        context="AAA",
        features=["a"],
        mean=[0.0],
        sd=[1.0],
        coefs=[1.0],
        intercept=0.0,
        ave=0.5,
    ),
    hl.Struct(
        context="ACG",
        features=["b"],
        mean=[0.0],
        sd=[2.0],
        coefs=[-1.0],
        intercept=0.0,
        ave=0.5,
    ),
]
"""Two one-feature models with hand-computable rates."""


class TestFitRegionalMutationRateModel:
    """Test the fit_regional_mutation_rate_model function."""

    def test_minimizes_the_penalized_log_loss(self) -> None:
        """
        Test that the fit minimizes mean log loss + penalty / 2 * sum(coef ** 2).

        The gradient of that objective, with the intercept unpenalized and the
        features standardized with the stored ``mean`` and ``sd``, is 0 at the fit.
        This checks the per-row penalty passed to sklearn and that the stored
        standardization is the one the model was fit with.
        """
        df = _training_df()
        penalty = 0.3
        model = fit_regional_mutation_rate_model(df, ["a", "b"], penalty)
        z = (df[["a", "b"]].to_numpy() - model["mean"]) / model["sd"]
        coefs = np.array(model["coefs"])
        residual = _expit(model["intercept"] + z @ coefs) - df["y"].to_numpy()
        assert residual.mean() == pytest.approx(0, abs=1e-4)
        assert z.T @ residual / len(df) + penalty * coefs == pytest.approx(
            [0, 0], abs=1e-4
        )
        assert model["sd"] == pytest.approx(list(df[["a", "b"]].std(ddof=0)))
        assert model["converged"]

    def test_ave(self) -> None:
        """Test that ave is the predicted probability with every feature at its mean."""
        model = fit_regional_mutation_rate_model(_training_df(), ["a", "b"], 0.3)
        assert model["ave"] == pytest.approx(_expit(model["intercept"]))

    def test_label_col(self) -> None:
        """Test that the labels can come from another column."""
        df = _training_df()
        df["label"] = df.pop("y")
        model = fit_regional_mutation_rate_model(df, ["a", "b"], 0.3, label_col="label")
        assert model["n_mutated"] == int(df["label"].sum())

    @pytest.mark.parametrize("value", [0.0, 0.4123])
    def test_constant_feature_is_left_out(self, value: float) -> None:
        """
        Test that a feature with a single value is left out of the model.

        Z-scoring a constant feature gives NaN. A constant like 0.4123 has a pandas
        sd of about 1e-17 rather than 0, so a check on the sd would keep it, and the
        rate at prediction would become 0 or 1 / ave.
        """
        df = _training_df().assign(c=value)
        model = fit_regional_mutation_rate_model(df, ["a", "c", "b"], 0.3)
        assert model["features"] == ["a", "b"]
        assert len(model["mean"]) == len(model["sd"]) == len(model["coefs"]) == 2
        assert np.isfinite(model["coefs"]).all()

    def test_rows_with_missing_features_are_dropped(self) -> None:
        """Test that rows missing any feature aren't counted or used."""
        df = _training_df()
        df.loc[:9, "a"] = np.nan
        model = fit_regional_mutation_rate_model(df, ["a", "b"], 0.3)
        assert model["n_mutated"] + model["n_background"] == len(df) - 10
        assert model["n_mutated"] == int(df["y"][10:].sum())

    @pytest.mark.parametrize("case", ["one_label", "no_varying_feature"])
    def test_no_model(self, case: str) -> None:
        """Test that there's no model without both labels or a varying feature."""
        df = _training_df()
        if case == "one_label":
            df = df[df["y"] == 0]
        else:
            df = df.assign(a=1.0, b=2.0)
        assert fit_regional_mutation_rate_model(df, ["a", "b"], 0.3) is None


class TestFitRegionalMutationRateModels:
    """Test the fit_regional_mutation_rate_models function."""

    def test_models_per_context(self) -> None:
        """
        Test that each context gets a model of its own rows and features.

        A context with no rows gets no model, and a feature constant in one context
        is left out of that context's model only.
        """
        df = pd.concat(
            [
                _training_df(seed=1).assign(ctx="AAA"),
                _training_df(seed=2).assign(ctx="ACG", b=0.0),
            ]
        )
        models = fit_regional_mutation_rate_models(
            df,
            {"AAA": ["a", "b"], "ACG": ["a", "b"], "CCC": ["a"]},
            0.3,
            context_col="ctx",
        )
        assert [m["context"] for m in models] == ["AAA", "ACG"]
        assert models[0]["features"] == ["a", "b"]
        assert models[1]["features"] == ["a"]
        expected = fit_regional_mutation_rate_model(
            df[df["ctx"] == "AAA"], ["a", "b"], 0.3
        )
        assert models[0]["coefs"] == pytest.approx(expected["coefs"])


class TestPredictRegionalMutationRateProbability:
    """Test the predict_regional_mutation_rate_probability function."""

    def test_prediction(self) -> None:
        """Test the prediction against the logistic model computed by hand."""
        model = {
            "features": ["b", "a"],
            "mean": [1.0, -1.0],
            "sd": [2.0, 0.5],
            "coefs": [0.3, -0.7],
            "intercept": -1.5,
        }
        df = pd.DataFrame({"a": [0.0, 1.0], "b": [3.0, -1.0]})
        eta = -1.5 + 0.3 * (df["b"] - 1.0) / 2.0 - 0.7 * (df["a"] + 1.0) / 0.5
        assert predict_regional_mutation_rate_probability(model, df) == pytest.approx(
            list(_expit(eta))
        )


class TestTuneRegionalMutationRatePenalty:
    """Test the tune_regional_mutation_rate_penalty function."""

    @staticmethod
    def _df() -> pd.DataFrame:
        """Return training rows of one context with a real effect, in two folds."""
        df = _training_df(n=8000, seed=2)
        df["context"] = "AAA"
        df["fold"] = ["odd" if i % 2 else "even" for i in range(len(df))]
        return df

    def test_held_out_log_loss(self) -> None:
        """
        Test that each fold is scored by a model fit on the other folds.

        The summed log loss of the odd fold matches that of a model fit on the even
        rows, computed by hand.
        """
        df = self._df()
        tuning_df = tune_regional_mutation_rate_penalty(df, {"AAA": ["a", "b"]}, [0.3])
        assert sorted(tuning_df["fold"]) == ["even", "odd"]
        odd = tuning_df.set_index("fold").loc["odd"]
        is_odd = df["fold"] == "odd"
        assert odd["n_test"] == is_odd.sum()

        model = fit_regional_mutation_rate_model(df[~is_odd], ["a", "b"], 0.3)
        p = predict_regional_mutation_rate_probability(model, df[is_odd])
        y = df["y"][is_odd].to_numpy()
        expected = -np.sum(y * np.log(p) + (1 - y) * np.log(1 - p))
        assert odd["log_loss"] == pytest.approx(expected, rel=1e-6)

    def test_rows_without_a_fold_are_only_fit(self) -> None:
        """Test that rows with a missing fold are used for fitting but not scored."""
        df = self._df()
        df.loc[df.index[:1000], "fold"] = None
        tuning_df = tune_regional_mutation_rate_penalty(df, {"AAA": ["a", "b"]}, [0.3])
        assert tuning_df["n_test"].sum() == len(df) - 1000

        model = fit_regional_mutation_rate_model(
            df[df["fold"] != "odd"], ["a", "b"], 0.3
        )
        test_df = df[df["fold"] == "odd"]
        p = predict_regional_mutation_rate_probability(model, test_df)
        y = test_df["y"].to_numpy()
        expected = -np.sum(y * np.log(p) + (1 - y) * np.log(1 - p))
        odd = tuning_df.set_index("fold").loc["odd"]
        assert odd["log_loss"] == pytest.approx(expected, rel=1e-6)

    def test_one_fold(self) -> None:
        """Test that nothing is scored when every row is in one fold."""
        df = self._df().assign(fold="even")
        tuning_df = tune_regional_mutation_rate_penalty(
            df, {"AAA": ["a", "b"]}, [0.01, 1000.0]
        )
        assert tuning_df.empty
        assert list(tuning_df.columns) == [
            "penalty",
            "context",
            "fold",
            "n_test",
            "log_loss",
        ]

    def test_too_strong_penalty_loses(self) -> None:
        """Test that a penalty that flattens a real effect has a higher held-out loss."""
        tuning_df = tune_regional_mutation_rate_penalty(
            self._df(), {"AAA": ["a", "b"]}, [0.01, 1000.0]
        )
        loss = tuning_df.groupby("penalty")["log_loss"].sum()
        assert loss[0.01] < loss[1000.0]


class TestChooseRegionalMutationRatePenalty:
    """Test the choose_regional_mutation_rate_penalty function."""

    TUNING_DF = pd.DataFrame(
        {
            "penalty": [0.1, 0.1, 0.3, 0.3, 1.0, 1.0],
            "context": ["AAA", "ACG"] * 3,
            "fold": ["odd"] * 6,
            "n_test": [10] * 6,
            "log_loss": [5.0, 3.0, 4.0, 3.5, 4.0, 4.0],
        }
    )
    """Summed losses of 8.0, 7.5, and 8.0 for penalties 0.1, 0.3, and 1.0."""

    def test_lowest_summed_loss(self, caplog: pytest.LogCaptureFixture) -> None:
        """Test that the penalty with the lowest loss summed over contexts is chosen."""
        with caplog.at_level(logging.WARNING):
            assert choose_regional_mutation_rate_penalty(self.TUNING_DF, 1.0) == 0.3
        assert not caplog.records

    def test_end_of_grid_warns(self, caplog: pytest.LogCaptureFixture) -> None:
        """Test that a chosen penalty at the end of those tried logs a warning."""
        tuning_df = self.TUNING_DF[self.TUNING_DF["penalty"] < 1.0]
        with caplog.at_level(logging.WARNING):
            assert choose_regional_mutation_rate_penalty(tuning_df, 1.0) == 0.3
        assert "at the end of those tried" in caplog.text

    def test_empty(self, caplog: pytest.LogCaptureFixture) -> None:
        """Test that the default penalty is returned, with a warning, when empty."""
        with caplog.at_level(logging.WARNING):
            assert (
                choose_regional_mutation_rate_penalty(self.TUNING_DF.iloc[:0], 0.3)
                == 0.3
            )
        assert "No held-out rows" in caplog.text


class TestRegionalMutationRateModelsToHt:
    """Test the regional_mutation_rate_models_to_ht function."""

    def test_round_trip(self) -> None:
        """Test that each model becomes a row keyed by context with its fields."""
        models = fit_regional_mutation_rate_models(
            _training_df().assign(context="AAA"), {"AAA": ["a", "b"]}, 0.3
        )
        ht = regional_mutation_rate_models_to_ht(models)
        assert list(ht.key) == ["context"]
        row = ht.collect()[0]
        assert row.context == "AAA"
        assert row.features == ["a", "b"]
        assert row.coefs == pytest.approx(models[0]["coefs"])
        assert row.n_mutated == models[0]["n_mutated"]


class TestRegionalMutationRateExpr:
    """Test the regional_mutation_rate_expr function."""

    def test_matches_the_fitted_model(self) -> None:
        """
        Test that the Hail rate equals the numpy prediction divided by ave.

        The model is read back from its Table, as a pipeline would. The two features
        have different means and sds, so the test covers the mean subtraction, each
        feature's sd and coefficient, and the feature order in both predictors.
        """
        models = fit_regional_mutation_rate_models(
            _training_df().assign(context="AAA"), {"AAA": ["a", "b"]}, 0.3
        )
        model = regional_mutation_rate_models_to_ht(models).collect()[0]
        values = [(-1.0, 9.0), (0.5, 5.0), (2.0, 1.0), (None, 1.0)]
        ht = _features_ht(
            [(1 + 1000 * i, 1001 + 1000 * i, a, b) for i, (a, b) in enumerate(values)]
        )
        ht = ht.annotate(rate=regional_mutation_rate_expr(ht.row, model))
        rates = [r.rate for r in ht.collect()]
        windows = pd.DataFrame(values[:3], columns=["a", "b"])
        expected = (
            predict_regional_mutation_rate_probability(models[0], windows)
            / models[0]["ave"]
        )
        assert rates[:3] == pytest.approx(list(expected))
        assert not np.allclose(expected, 1.0)
        assert rates[3] is None


class TestComputeRegionalMutationRateHt:
    """Test the compute_regional_mutation_rate_ht function."""

    def test_rates(self) -> None:
        """
        Test the rate per window, the 1.0 fill, and which windows are kept.

        A window missing a context's feature gets 1.0 for that context, contexts
        without a model get 1.0, and a window with no rate for any modeled context
        is removed.
        """
        ht = _features_ht(
            [(1, 1001, 0.0, 2.0), (1001, 2001, 1.0, None), (2001, 3001, None, None)]
        )
        ht = compute_regional_mutation_rate_ht(ht, MODELS, ["AAA", "ACG", "CCC"])
        assert list(ht.row) == ["interval", "adj_r"]
        rows = {r.interval.start.position: r.adj_r for r in ht.collect()}
        assert set(rows) == {1, 1001}
        assert set(rows[1]) == {"AAA", "ACG", "CCC"}
        assert rows[1]["AAA"] == pytest.approx(1.0)
        assert rows[1]["ACG"] == pytest.approx(_expit(-1.0) / 0.5)
        assert rows[1]["CCC"] == 1.0
        assert rows[1001]["AAA"] == pytest.approx(_expit(1.0) / 0.5)
        assert rows[1001]["ACG"] == 1.0

    def test_no_models(self) -> None:
        """Test that no models give an empty Table with the same fields."""
        ht = compute_regional_mutation_rate_ht(
            _features_ht([(1, 1001, 0.0, 2.0)]), [], ["AAA"]
        )
        assert ht.row.dtype == hl.tstruct(
            interval=INTERVAL_TYPE, adj_r=hl.tdict(hl.tstr, hl.tfloat64)
        )
        assert ht.count() == 0


def _rate_ht() -> hl.Table:
    """Return two windows with rates for AAA, ACG, and CCC."""
    return hl.Table.parallelize(
        [
            {
                "interval": _window(1, 100),
                "adj_r": {"AAA": 2.0, "ACG": 0.5, "CCC": 3.0},
            },
            {
                "interval": _window(100, 200),
                "adj_r": {"AAA": 1.0, "ACG": 1.0, "CCC": 1.0},
            },
        ],
        hl.tstruct(interval=INTERVAL_TYPE, adj_r=hl.tdict(hl.tstr, hl.tfloat64)),
        key="interval",
    )


def _sites_ht() -> hl.Table:
    """
    Return weighted sites, with one outside the windows and some without a weight.

    The site with a missing weight and the one with weight 0 are in the first
    window, which has rate 2.0 for AAA, so counting either changes the factor.
    """
    ht = hl.Table.parallelize(
        [
            {"pos": 10, "context": "AAA", "weight": 1.0},
            {"pos": 150, "context": "AAA", "weight": 3.0},
            {"pos": 20, "context": "ACG", "weight": 2.0},
            {"pos": 500, "context": "AAA", "weight": 10.0},
            {"pos": 30, "context": "AAA", "weight": None},
            {"pos": 40, "context": "AAA", "weight": 0.0},
        ],
        hl.tstruct(pos=hl.tint32, context=hl.tstr, weight=hl.tfloat64),
    )
    return ht.key_by(locus=hl.locus("chr1", ht.pos, "GRCh38"))


def _normalization_factors(sites_ht: hl.Table, rate_ht: hl.Table) -> dict:
    """Return the normalization factors of the sites' rates, weighted by ``weight``."""
    return sites_ht.aggregate(
        regional_mutation_rate_normalization_factors_agg_expr(
            sites_ht.context,
            rate_ht[sites_ht.locus].adj_r.get(sites_ht.context),
            sites_ht.weight,
        )
    )


class TestRegionalMutationRateNormalizationFactorsAggExpr:
    """Test the regional_mutation_rate_normalization_factors_agg_expr function."""

    def test_factors(self) -> None:
        """
        Test that a factor is the weighted mean rate of the sites in a window.

        Sites outside every window, with a missing weight, or with weight 0 are
        left out, and contexts without sites get no factor.
        """
        factors = _normalization_factors(_sites_ht(), _rate_ht())
        assert factors == pytest.approx({"AAA": (1 * 2 + 3 * 1) / 4, "ACG": 0.5})


class TestNormalizeRegionalMutationRateExpr:
    """Test the normalize_regional_mutation_rate_expr function."""

    def test_normalized_mean_is_one(self) -> None:
        """Test that normalized rates average 1, and contexts without a factor stay."""
        sites_ht = _sites_ht()
        rate_ht = _rate_ht()
        factors = _normalization_factors(sites_ht, rate_ht)
        rate_ht = rate_ht.annotate(
            adj_r=normalize_regional_mutation_rate_expr(rate_ht.adj_r, factors)
        )
        rates = {r.interval.start.position: r.adj_r for r in rate_ht.collect()}
        assert (1 * rates[1]["AAA"] + 3 * rates[100]["AAA"]) / 4 == pytest.approx(1.0)
        assert rates[1]["ACG"] == pytest.approx(1.0)
        assert rates[1]["CCC"] == 3.0
