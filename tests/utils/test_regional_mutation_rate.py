"""Tests for the regional mutation rate utility module."""

import logging
from typing import Dict, List, Optional

import hail as hl
import numpy as np
import pytest

from gnomad.utils.constraint import (
    choose_regional_mutation_rate_penalty,
    compute_regional_mutation_rate_ht,
    fit_regional_mutation_rate_models,
    normalize_regional_mutation_rate_expr,
    regional_mutation_rate_normalization_factors_agg_expr,
    tune_regional_mutation_rate_penalty,
)
from gnomad.utils.constraint.regional_mutation_rate import MODELS_HT_SCHEMA

INTERVAL_TYPE = hl.tinterval(hl.tlocus("GRCh38"))
"""Type of the window intervals in the test feature and rate Tables."""

ROW_TYPE = hl.tstruct(
    y=hl.tint32,
    a=hl.tfloat64,
    b=hl.tfloat64,
    c=hl.tfloat64,
    context=hl.tstr,
    fold=hl.tstr,
)
"""Row type of the test training Tables."""


def _expit(x: np.ndarray) -> np.ndarray:
    """Return the logistic function of ``x``."""
    return 1.0 / (1.0 + np.exp(-x))


def _window(start: int, end: int) -> hl.Interval:
    """Return the chr1 interval ``[start, end)``."""
    return hl.Interval(
        hl.Locus("chr1", start, "GRCh38"), hl.Locus("chr1", end, "GRCh38")
    )


def _training_rows(n: int = 2000, seed: int = 0, context: str = "AAA") -> List[Dict]:
    """
    Return training rows whose labels follow a logistic model of ``a`` and ``b``.

    ``c`` is 0, and the rows alternate between the folds ``even`` and ``odd``.
    """
    rng = np.random.default_rng(seed)
    a = rng.normal(size=n)
    b = 3 * rng.normal(size=n) + 5
    y = rng.random(n) < _expit(-2 + 0.8 * a - 0.5 * (b - 5) / 3)
    return [
        {
            "y": int(y[i]),
            "a": float(a[i]),
            "b": float(b[i]),
            "c": 0.0,
            "context": context,
            "fold": "odd" if i % 2 else "even",
        }
        for i in range(n)
    ]


def _training_ht(rows: List[Dict]) -> hl.Table:
    """Return a training Table of the rows."""
    return hl.Table.parallelize(rows, ROW_TYPE)


def _arrays(rows: List[Dict], features: List[str]) -> tuple:
    """Return the labels and the feature matrix of the rows, for numpy checks."""
    return (
        np.array([r["y"] for r in rows], dtype=float),
        np.array([[r[f] for f in features] for r in rows], dtype=float),
    )


def _models(ht: hl.Table) -> Dict[str, hl.Struct]:
    """Return the models of a models Table by context."""
    return {m.context: m for m in ht.collect()}


def _predict(model: hl.Struct, x: np.ndarray) -> np.ndarray:
    """Return a model's predicted probabilities for a feature matrix, with numpy."""
    z = (x - np.array(model.mean)) / np.array(model.sd)
    return _expit(model.intercept + z @ np.array(model.coefs))


def _log_loss(y: np.ndarray, p: np.ndarray) -> float:
    """Return the summed log loss of the labels under the predicted probabilities."""
    return float(-np.sum(y * np.log(p) + (1 - y) * np.log(1 - p)))


class TestFitRegionalMutationRateModels:
    """Test the fit_regional_mutation_rate_models function."""

    @staticmethod
    def _gradient(rows: List[Dict], model: hl.Struct, penalty: float) -> np.ndarray:
        """Return the gradient of the penalized mean log loss at a model, with numpy."""
        y, x = _arrays(rows, model.features)
        z = (x - np.array(model.mean)) / np.array(model.sd)
        residual = _predict(model, x) - y
        coefs = np.array(model.coefs)
        return np.r_[residual.mean(), z.T @ residual / len(y) + penalty * coefs]

    @pytest.mark.parametrize("penalty", [0.3, 0.0])
    def test_minimizes_the_penalized_log_loss(self, penalty: float) -> None:
        """
        Test that the fit minimizes mean log loss + penalty / 2 * sum(coef ** 2).

        The gradient of that objective, with the intercept unpenalized and the
        features standardized with the stored ``mean`` and ``sd`` (ddof=0), is 0 at
        the fit. A penalty of 0 gives the unpenalized fit.
        """
        rows = _training_rows()
        model = _models(
            fit_regional_mutation_rate_models(
                _training_ht(rows), {"AAA": ["a", "b"]}, penalty
            )
        )["AAA"]
        assert self._gradient(rows, model, penalty) == pytest.approx(
            [0, 0, 0], abs=1e-8
        )
        _, x = _arrays(rows, ["a", "b"])
        assert model.mean == pytest.approx(list(x.mean(axis=0)))
        assert model.sd == pytest.approx(list(x.std(axis=0)))
        assert model.converged
        assert model.ave == pytest.approx(_expit(model.intercept))

    def test_counts_and_expressions(self) -> None:
        """Test the label counts, and labels and contexts from other expressions."""
        rows = _training_rows()
        ht = _training_ht(rows).rename({"y": "label", "context": "ctx"})
        model = _models(
            fit_regional_mutation_rate_models(
                ht,
                {"AAA": ["a", "b"]},
                0.3,
                label_expr=ht.label,
                context_expr=ht.ctx,
            )
        )["AAA"]
        n_mutated = sum(r["y"] for r in rows)
        assert (model.n_mutated, model.n_background) == (n_mutated, len(rows) - n_mutated)

    @pytest.mark.parametrize("value", [0.0, 0.4123])
    def test_constant_feature_is_left_out(self, value: float) -> None:
        """
        Test that a feature with a single value is left out of the model.

        Z-scoring a constant feature gives NaN, and a constant like 0.4123 can have a
        computed sd of about 1e-17 rather than 0. The kept features' mean, sd, and
        coefficients stay aligned with them.
        """
        rows = [{**r, "c": value} for r in _training_rows()]
        ht = _training_ht(rows)
        model = _models(fit_regional_mutation_rate_models(ht, {"AAA": ["a", "c", "b"]}, 0.3))[
            "AAA"
        ]
        expected = _models(fit_regional_mutation_rate_models(ht, {"AAA": ["a", "b"]}, 0.3))[
            "AAA"
        ]
        assert model.features == ["a", "b"]
        for field in ["mean", "sd", "coefs"]:
            assert model[field] == pytest.approx(expected[field])

    def test_rows_with_missing_features_are_dropped(self) -> None:
        """
        Test that rows missing a model feature aren't counted or used.

        A missing value in a field that isn't a feature doesn't drop the row.
        """
        rows = _training_rows()
        rows = [
            {**r, "a": None} if i < 10 else {**r, "c": None}
            for i, r in enumerate(rows)
        ]
        model = _models(
            fit_regional_mutation_rate_models(_training_ht(rows), {"AAA": ["a", "b"]}, 0.3)
        )["AAA"]
        assert model.n_mutated + model.n_background == len(rows) - 10
        assert model.n_mutated == sum(r["y"] for r in rows[10:])

    @pytest.mark.parametrize(
        "case", ["no_mutated", "no_background", "no_varying_feature"]
    )
    def test_no_model(self, case: str) -> None:
        """Test that there's no model without both labels or a varying feature."""
        rows = _training_rows()
        if case == "no_mutated":
            rows = [r for r in rows if r["y"] == 0]
        elif case == "no_background":
            rows = [r for r in rows if r["y"] == 1]
        else:
            rows = [{**r, "a": 1.0, "b": 2.0} for r in rows]
        ht = fit_regional_mutation_rate_models(_training_ht(rows), {"AAA": ["a", "b"]}, 0.3)
        assert ht.count() == 0

    def test_not_converged(self) -> None:
        """Test that a fit stopped at ``max_iter`` is flagged as not converged."""
        model = _models(
            fit_regional_mutation_rate_models(
                _training_ht(_training_rows()), {"AAA": ["a", "b"]}, 0.3, max_iter=1
            )
        )["AAA"]
        assert not model.converged

    def test_models_per_context(self) -> None:
        """
        Test that each context gets a model of its own rows and features.

        A context with no rows gets no model, and a feature constant in one context
        is left out of that context's model only.
        """
        aaa = _training_rows(seed=1)
        acg = [{**r, "b": 0.0} for r in _training_rows(seed=2, context="ACG")]
        models = _models(
            fit_regional_mutation_rate_models(
                _training_ht(aaa + acg),
                {"AAA": ["a", "b"], "ACG": ["a", "b"], "CCC": ["a"]},
                0.3,
            )
        )
        assert sorted(models) == ["AAA", "ACG"]
        assert models["AAA"].features == ["a", "b"]
        assert models["ACG"].features == ["a"]
        expected = _models(
            fit_regional_mutation_rate_models(_training_ht(aaa), {"AAA": ["a", "b"]}, 0.3)
        )["AAA"]
        assert models["AAA"].coefs == pytest.approx(expected.coefs)


class TestTuneRegionalMutationRatePenalty:
    """Test the tune_regional_mutation_rate_penalty function."""

    @staticmethod
    def _expected_loss(
        rows: List[Dict], fold: str, penalty: float, features: List[str]
    ) -> float:
        """Return a fold's summed log loss under a model fit on the other rows."""
        fit_rows = [r for r in rows if r["fold"] != fold]
        model = _models(
            fit_regional_mutation_rate_models(
                _training_ht(fit_rows), {rows[0]["context"]: features}, penalty
            )
        )[rows[0]["context"]]
        y, x = _arrays([r for r in rows if r["fold"] == fold], features)
        return _log_loss(y, _predict(model, x))

    def test_held_out_log_loss(self) -> None:
        """
        Test that each fold is scored by a model fit on the other folds.

        Rows with a missing fold are used for fitting but not scored.
        """
        rows = _training_rows(seed=2)
        rows = [{**r, "fold": None} if i < 200 else r for i, r in enumerate(rows)]
        tuning = {
            (r.fold, r.context): r
            for r in tune_regional_mutation_rate_penalty(
                _training_ht(rows), {"AAA": ["a", "b"]}, [0.3]
            ).collect()
        }
        assert sorted(tuning) == [("even", "AAA"), ("odd", "AAA")]
        assert sum(r.n_test for r in tuning.values()) == len(rows) - 200
        assert tuning[("odd", "AAA")].log_loss == pytest.approx(
            self._expected_loss(rows, "odd", 0.3, ["a", "b"]), rel=1e-6
        )

    def test_contexts_and_expressions(self) -> None:
        """
        Test that each context is tuned on its own rows, with other expressions.

        A held-out fold with only one label is still scored.
        """
        aaa = _training_rows(seed=2)
        acg = _training_rows(seed=3, context="ACG")
        only_background = [
            {**r, "fold": "only_background"} for r in acg if r["y"] == 0
        ][:100]
        ht = _training_ht(aaa + acg + only_background)
        ht = ht.rename({"y": "label", "context": "ctx", "fold": "split"})
        tuning = {
            (r.context, r.fold): r
            for r in tune_regional_mutation_rate_penalty(
                ht,
                {"AAA": ["a", "b"], "ACG": ["a", "b"]},
                [0.3],
                label_expr=ht.label,
                context_expr=ht.ctx,
                fold_expr=ht.split,
            ).collect()
        }
        assert sorted(tuning) == [
            ("AAA", "even"),
            ("AAA", "odd"),
            ("ACG", "even"),
            ("ACG", "odd"),
            ("ACG", "only_background"),
        ]
        assert tuning[("AAA", "odd")].log_loss == pytest.approx(
            self._expected_loss(aaa, "odd", 0.3, ["a", "b"]), rel=1e-6
        )
        acg_rows = acg + only_background
        row = tuning[("ACG", "only_background")]
        assert row.n_test == 100
        assert row.log_loss == pytest.approx(
            self._expected_loss(acg_rows, "only_background", 0.3, ["a", "b"]),
            rel=1e-6,
        )

    def test_one_fold(self) -> None:
        """Test that nothing is scored when every row is in one fold."""
        rows = [{**r, "fold": "even"} for r in _training_rows()]
        ht = tune_regional_mutation_rate_penalty(
            _training_ht(rows), {"AAA": ["a", "b"]}, [0.01, 1000.0]
        )
        assert ht.count() == 0
        assert list(ht.row) == ["penalty", "context", "fold", "n_test", "log_loss"]

    def test_too_strong_penalty_loses(self) -> None:
        """Test that a penalty that flattens a real effect has a higher held-out loss."""
        ht = tune_regional_mutation_rate_penalty(
            _training_ht(_training_rows(seed=2)), {"AAA": ["a", "b"]}, [0.01, 1000.0]
        )
        loss = ht.aggregate(hl.agg.group_by(ht.penalty, hl.agg.sum(ht.log_loss)))
        assert loss[0.01] < loss[1000.0]


class TestChooseRegionalMutationRatePenalty:
    """Test the choose_regional_mutation_rate_penalty function."""

    ROWS = [
        (0.1, "AAA", 5.0),
        (0.1, "ACG", 3.0),
        (0.3, "AAA", 4.0),
        (0.3, "ACG", 3.5),
        (1.0, "AAA", 4.0),
        (1.0, "ACG", 4.0),
    ]
    """Summed losses of 8.0, 7.5, and 8.0 for penalties 0.1, 0.3, and 1.0."""

    @classmethod
    def _tuning_ht(cls, penalties: Optional[List[float]] = None) -> hl.Table:
        """Return a tuning Table with the rows of the given penalties."""
        return hl.Table.parallelize(
            [
                {"penalty": p, "context": c, "fold": "odd", "n_test": 10, "log_loss": l}
                for p, c, l in cls.ROWS
                if penalties is None or p in penalties
            ],
            hl.tstruct(
                penalty=hl.tfloat64,
                context=hl.tstr,
                fold=hl.tstr,
                n_test=hl.tint64,
                log_loss=hl.tfloat64,
            ),
            key=["penalty", "context", "fold"],
        )

    @staticmethod
    def _warnings(caplog: pytest.LogCaptureFixture) -> list:
        """Return the warning messages this module logged."""
        return [
            r.getMessage() for r in caplog.records if r.name == "regional_mutation_rate"
        ]

    def test_lowest_summed_loss(self, caplog: pytest.LogCaptureFixture) -> None:
        """Test that the penalty with the lowest loss summed over contexts is chosen."""
        with caplog.at_level(logging.WARNING):
            assert choose_regional_mutation_rate_penalty(self._tuning_ht(), 1.0) == 0.3
        assert not self._warnings(caplog)

    @pytest.mark.parametrize("penalties", [[0.1, 0.3], [0.3, 1.0]])
    def test_end_of_grid_warns(
        self, caplog: pytest.LogCaptureFixture, penalties: List[float]
    ) -> None:
        """Test that a chosen penalty at either end of those tried logs a warning."""
        with caplog.at_level(logging.WARNING):
            assert (
                choose_regional_mutation_rate_penalty(self._tuning_ht(penalties), 1.0)
                == 0.3
            )
        assert any("at the end of those tried" in m for m in self._warnings(caplog))

    def test_one_penalty(self, caplog: pytest.LogCaptureFixture) -> None:
        """Test that a single penalty tried is returned without a warning."""
        with caplog.at_level(logging.WARNING):
            assert (
                choose_regional_mutation_rate_penalty(self._tuning_ht([0.3]), 1.0) == 0.3
            )
        assert not self._warnings(caplog)

    def test_empty(self, caplog: pytest.LogCaptureFixture) -> None:
        """Test that the default penalty is returned, with a warning, when empty."""
        with caplog.at_level(logging.WARNING):
            assert choose_regional_mutation_rate_penalty(self._tuning_ht([]), 0.3) == 0.3
        assert any("No held-out rows" in m for m in self._warnings(caplog))


def _features_ht(rows: list) -> hl.Table:
    """Return a feature Table keyed by window with float features ``a`` and ``b``."""
    return hl.Table.parallelize(
        [{"interval": _window(s, e), "a": a, "b": b} for s, e, a, b in rows],
        hl.tstruct(interval=INTERVAL_TYPE, a=hl.tfloat64, b=hl.tfloat64),
        key="interval",
    )


def _models_ht() -> hl.Table:
    """Return two one-feature models with hand-computable rates."""
    common = {"intercept": 0.0, "ave": 0.5, "n_mutated": 10, "n_background": 10}
    return hl.Table.parallelize(
        [
            {
                "context": "AAA",
                "features": ["a"],
                "mean": [0.0],
                "sd": [1.0],
                "coefs": [1.0],
                "converged": True,
                **common,
            },
            {
                "context": "ACG",
                "features": ["b"],
                "mean": [0.0],
                "sd": [2.0],
                "coefs": [-1.0],
                "converged": True,
                **common,
            },
        ],
        MODELS_HT_SCHEMA,
        key="context",
    )


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
            [(1, 1001, 0.5, 2.0), (1001, 2001, 1.0, None), (2001, 3001, None, None)]
        )
        ht = compute_regional_mutation_rate_ht(ht, _models_ht(), ["AAA", "ACG", "CCC"])
        assert list(ht.row) == ["interval", "adj_r"]
        assert list(ht.key) == ["interval"]
        rows = {r.interval.start.position: r.adj_r for r in ht.collect()}
        assert set(rows) == {1, 1001}
        assert set(rows[1]) == {"AAA", "ACG", "CCC"}
        assert rows[1]["AAA"] == pytest.approx(_expit(0.5) / 0.5)
        assert rows[1]["ACG"] == pytest.approx(_expit(-1.0) / 0.5)
        assert rows[1]["CCC"] == 1.0
        assert rows[1001]["AAA"] == pytest.approx(_expit(1.0) / 0.5)
        assert rows[1001]["ACG"] == 1.0

    def test_only_requested_contexts(self) -> None:
        """
        Test that models of contexts that weren't asked for are ignored.

        The second window has an AAA rate but no ACG rate, so with only ACG asked
        for it is removed rather than kept with ACG filled with 1.0.
        """
        ht = _features_ht([(1, 1001, 0.5, 2.0), (1001, 2001, 1.0, None)])
        ht = compute_regional_mutation_rate_ht(ht, _models_ht(), ["ACG"])
        rows = {r.interval.start.position: r.adj_r for r in ht.collect()}
        assert rows == {1: pytest.approx({"ACG": _expit(-1.0) / 0.5})}

    def test_no_models(self) -> None:
        """Test that no models give an empty Table with the same fields and key."""
        ht = compute_regional_mutation_rate_ht(
            _features_ht([(1, 1001, 0.0, 2.0)]), _models_ht().filter(False), ["AAA"]
        )
        assert ht.row.dtype == hl.tstruct(
            interval=INTERVAL_TYPE, adj_r=hl.tdict(hl.tstr, hl.tfloat64)
        )
        assert list(ht.key) == ["interval"]
        assert ht.count() == 0

    def test_matches_the_fitted_model(self) -> None:
        """
        Test that the rate is the fitted model's prediction divided by ``ave``.

        The two features have different means and sds, so the test covers the mean
        subtraction, each feature's sd and coefficient, and the feature order. A
        window at the training means has rate 1, and a window missing a feature has
        no rate, so it is removed.
        """
        models_ht = fit_regional_mutation_rate_models(
            _training_ht(_training_rows()), {"AAA": ["a", "b"]}, 0.3
        )
        model = _models(models_ht)["AAA"]
        values = [(-1.0, 9.0), (0.5, 5.0), (2.0, 1.0), tuple(model.mean), (None, 1.0)]
        ht = _features_ht(
            [(1 + 1000 * i, 1001 + 1000 * i, a, b) for i, (a, b) in enumerate(values)]
        )
        rates = [
            r.adj_r["AAA"]
            for r in compute_regional_mutation_rate_ht(ht, models_ht, ["AAA"]).collect()
        ]
        expected = _predict(model, np.array(values[:4])) / model.ave
        assert rates == pytest.approx(list(expected))
        assert rates[3] == pytest.approx(1.0)
        assert not np.allclose(expected[:3], 1.0)


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

    The CCC sites are in the first window, but one has weight 0 and the other no
    weight, so CCC has no positively weighted site to normalize on.
    """
    ht = hl.Table.parallelize(
        [
            {"pos": 10, "context": "AAA", "weight": 1.0},
            {"pos": 150, "context": "AAA", "weight": 3.0},
            {"pos": 20, "context": "ACG", "weight": 2.0},
            {"pos": 500, "context": "AAA", "weight": 10.0},
            {"pos": 30, "context": "AAA", "weight": None},
            {"pos": 40, "context": "AAA", "weight": 0.0},
            {"pos": 50, "context": "CCC", "weight": 0.0},
            {"pos": 60, "context": "CCC", "weight": None},
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
        left out. A context with no positively weighted site in a window (CCC) gets
        no factor rather than 0 / 0.
        """
        factors = _normalization_factors(_sites_ht(), _rate_ht())
        assert factors == pytest.approx({"AAA": (1 * 2 + 3 * 1) / 4, "ACG": 0.5})


class TestNormalizeRegionalMutationRateExpr:
    """Test the normalize_regional_mutation_rate_expr function."""

    @pytest.mark.parametrize("as_expr", [False, True])
    def test_normalized_mean_is_one(self, as_expr: bool) -> None:
        """
        Test that normalized rates average 1, and contexts without a factor stay.

        The factors can be a dict or a dict expression.
        """
        sites_ht = _sites_ht()
        rate_ht = _rate_ht()
        factors = _normalization_factors(sites_ht, rate_ht)
        if as_expr:
            factors = hl.literal(factors, hl.tdict(hl.tstr, hl.tfloat64))
        rate_ht = rate_ht.annotate(
            adj_r=normalize_regional_mutation_rate_expr(rate_ht.adj_r, factors)
        )
        rates = {r.interval.start.position: r.adj_r for r in rate_ht.collect()}
        assert (1 * rates[1]["AAA"] + 3 * rates[100]["AAA"]) / 4 == pytest.approx(1.0)
        assert rates[1]["ACG"] == pytest.approx(1.0)
        assert rates[1]["CCC"] == 3.0
