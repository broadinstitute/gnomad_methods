"""Tests for the regional mutation rate utility module."""

import logging
from typing import Dict, List, Optional, Union

import hail as hl
import numpy as np
import pytest

from gnomad.utils.constraint import (
    MODELS_HT_SCHEMA,
    choose_regional_mutation_rate_penalty,
    fit_regional_mutation_rate_models,
    normalize_regional_mutation_rate_expr,
    prepare_regional_mutation_rate_training_ht,
    regional_mutation_rate_expr,
    regional_mutation_rate_normalization_factors_agg_expr,
    tune_regional_mutation_rate_penalty,
)
from gnomad.utils.constraint.regional_mutation_rate import _newton_update

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


def _heavy_tailed_rows() -> List[Dict]:
    """
    Return 50 rows of three lognormal features, with 2 of them mutated.

    The features have an sd of 2.5 on the log scale, so a few standardized values are
    far out. At a penalty of 0.01, full Newton steps from the start overshoot on
    these rows and never converge.
    """
    rng = np.random.default_rng(23)
    x = np.exp(2.5 * rng.normal(size=(50, 3)))
    y = rng.random(50) < 0.05
    return [
        {
            "y": int(y[i]),
            "a": float(x[i, 0]),
            "b": float(x[i, 1]),
            "c": float(x[i, 2]),
            "context": "AAA",
            "fold": None,
        }
        for i in range(50)
    ]


def _training_ht(rows: List[Dict]) -> hl.Table:
    """Return a training Table of the rows."""
    return hl.Table.parallelize(rows, ROW_TYPE)


def _prepared(
    rows: List[Dict],
    features: Union[List[str], Dict[str, List[str]]],
    fold: bool = False,
) -> hl.Table:
    """Return the prepared training Table of the rows, with their fold if asked."""
    ht = _training_ht(rows)
    return prepare_regional_mutation_rate_training_ht(
        ht, features, fold_expr=ht.fold if fold else None
    )


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


def _warnings(caplog: pytest.LogCaptureFixture) -> list:
    """Return the warning messages the constraint utilities logged."""
    return [
        r.getMessage()
        for r in caplog.records
        if r.name == "constraint_utils" and r.levelno == logging.WARNING
    ]


class TestPrepareRegionalMutationRateTrainingHt:
    """Test the prepare_regional_mutation_rate_training_ht function."""

    ROWS = [
        {"y": 1, "a": 1.0, "b": 2.0, "c": 0.0, "context": "AAA", "fold": "odd"},
        {"y": 0, "a": None, "b": 3.0, "c": 0.0, "context": "AAA", "fold": None},
        {"y": 0, "a": 1.0, "b": 4.0, "c": 0.0, "context": "ACG", "fold": "even"},
        {"y": 1, "a": 1.0, "b": 5.0, "c": 0.0, "context": "CCC", "fold": "odd"},
    ]
    """Rows of three contexts, one with a missing feature and one with no fold."""

    def test_features_per_context(self) -> None:
        """
        Test each row's label, context, fold, and context's features, and the globals.

        Rows with a missing value are kept, and a context without features has no
        ``x``.
        """
        features = {"AAA": ["b", "a"], "ACG": ["b"]}
        ht = _prepared(self.ROWS, features, fold=True)
        assert list(ht.row) == ["label", "context", "fold", "x"]
        assert hl.eval(ht.features) == features
        assert [(r.label, r.context, r.fold, r.x) for r in ht.collect()] == [
            (1.0, "AAA", "odd", [2.0, 1.0]),
            (0.0, "AAA", None, [3.0, None]),
            (0.0, "ACG", "even", [4.0]),
            (1.0, "CCC", "odd", None),
        ]

    def test_one_feature_list(self) -> None:
        """Test that one list of features gives every context the same features."""
        ht = _prepared(self.ROWS, ["b", "a"])
        assert list(ht.row) == ["label", "context", "x"]
        assert hl.eval(ht.features) == ["b", "a"]
        assert [r.x for r in ht.collect()] == [
            [2.0, 1.0],
            [3.0, None],
            [4.0, 1.0],
            [5.0, 1.0],
        ]

    def test_feature_types(self) -> None:
        """
        Test that int and boolean features become float64 and string features raise.

        ``hl.float64`` would parse a string field, so it is rejected up front.
        """
        ht = hl.Table.parallelize(
            [{"y": True, "n": 2, "flag": True, "s": "1.5", "context": "AAA"}],
            hl.tstruct(
                y=hl.tbool, n=hl.tint32, flag=hl.tbool, s=hl.tstr, context=hl.tstr
            ),
        )
        prepared = prepare_regional_mutation_rate_training_ht(ht, ["n", "flag"])
        assert prepared.x.dtype == hl.tarray(hl.tfloat64)
        assert [(r.label, r.x) for r in prepared.collect()] == [(1.0, [2.0, 1.0])]
        with pytest.raises(ValueError, match="these are not: s"):
            prepare_regional_mutation_rate_training_ht(ht, {"AAA": ["n", "s"]})

    def test_fit_and_tune_need_a_prepared_table(self) -> None:
        """Test that the fit and the tuning say how to prepare a Table they can't use."""
        with pytest.raises(ValueError, match="prepare_regional_mutation_rate"):
            fit_regional_mutation_rate_models(_training_ht(_training_rows()), 0.3)
        with pytest.raises(ValueError, match="no fold.*with a fold_expr"):
            tune_regional_mutation_rate_penalty(
                _prepared(_training_rows(), {"AAA": ["a", "b"]}), [0.3]
            )


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

    @pytest.mark.parametrize("penalty", [0.3, 0.01])
    def test_minimizes_the_penalized_log_loss(self, penalty: float) -> None:
        """
        Test that the fit minimizes mean log loss + penalty / 2 * sum(coef ** 2).

        The gradient of that objective, with the intercept unpenalized and the
        features standardized with the stored ``mean`` and ``sd`` (ddof=0), is 0 at
        the fit.
        """
        rows = _training_rows()
        model = _models(
            fit_regional_mutation_rate_models(
                _prepared(rows, {"AAA": ["a", "b"]}), penalty
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

    def test_step_halving(self) -> None:
        """
        Test that the fit reaches the minimum where full Newton steps never converge.

        On these heavy-tailed rows, full steps overshoot and cycle, so a step that
        raises the objective has to be halved.
        """
        rows = _heavy_tailed_rows()
        model = _models(
            fit_regional_mutation_rate_models(
                _prepared(rows, {"AAA": ["a", "b", "c"]}), 0.01
            )
        )["AAA"]
        assert model.converged
        assert self._gradient(rows, model, 0.01) == pytest.approx(
            [0, 0, 0, 0], abs=1e-8
        )

    def test_penalty_must_be_positive(self) -> None:
        """
        Test that the fit and the tuning reject a penalty of 0.

        With no penalty, full Newton steps on a heavy-tailed feature ran off until
        every prediction was exactly 0 or 1, where the fit looked converged.
        """
        ht = _prepared(_training_rows(), {"AAA": ["a", "b"]}, fold=True)
        with pytest.raises(ValueError, match="greater than 0"):
            fit_regional_mutation_rate_models(ht, 0.0)
        with pytest.raises(ValueError, match="greater than 0"):
            tune_regional_mutation_rate_penalty(ht, [0.3, 0.0])

    def test_counts_and_expressions(self) -> None:
        """Test the label counts, and labels and contexts from other expressions."""
        rows = _training_rows()
        ht = _training_ht(rows).rename({"y": "label", "context": "ctx"})
        model = _models(
            fit_regional_mutation_rate_models(
                prepare_regional_mutation_rate_training_ht(
                    ht, {"AAA": ["a", "b"]}, label_expr=ht.label, context_expr=ht.ctx
                ),
                0.3,
            )
        )["AAA"]
        n_mutated = sum(r["y"] for r in rows)
        assert (model.n_mutated, model.n_background) == (
            n_mutated,
            len(rows) - n_mutated,
        )

    @pytest.mark.parametrize("value", [0.0, 0.4123])
    def test_constant_feature_is_left_out(self, value: float) -> None:
        """
        Test that a feature with a single value is left out of the model.

        Standardizing a constant feature divides by an sd of 0, or, for a constant
        like 0.4123 whose computed mean differs from it by rounding, by a tiny sd.
        The kept features' mean, sd, and coefficients stay aligned with them.
        """
        rows = [{**r, "c": value} for r in _training_rows()]
        model = _models(
            fit_regional_mutation_rate_models(
                _prepared(rows, {"AAA": ["a", "c", "b"]}), 0.3
            )
        )["AAA"]
        expected = _models(
            fit_regional_mutation_rate_models(_prepared(rows, {"AAA": ["a", "b"]}), 0.3)
        )["AAA"]
        assert model.features == ["a", "b"]
        for field in ["mean", "sd", "coefs"]:
            assert model[field] == pytest.approx(expected[field])

    def test_nearly_constant_feature(self) -> None:
        """
        Test that a feature that varies little next to its mean gets its exact sd.

        An sd from the sum of squares, in the same pass as the mean, is NaN for these
        values, which made the fit fail for every context.
        """
        rows = [
            {**r, "c": 0.4123 + (1e-9 if i == 0 else 0.0)}
            for i, r in enumerate(_training_rows())
        ]
        model = _models(
            fit_regional_mutation_rate_models(
                _prepared(rows, {"AAA": ["a", "b", "c"]}), 0.3
            )
        )["AAA"]
        _, x = _arrays(rows, ["a", "b", "c"])
        assert model.features == ["a", "b", "c"]
        assert model.sd == pytest.approx(list(x.std(axis=0)), rel=1e-6)
        assert model.converged

    def test_rows_with_missing_values_are_dropped(self) -> None:
        """
        Test that rows missing the label or a model feature aren't counted or used.

        A missing value in a field that isn't a feature doesn't drop the row.
        """
        rows = _training_rows()
        for r in rows[:10]:
            r["a"] = None
        for r in rows[10:15]:
            r["y"] = None
        for r in rows[15:]:
            r["c"] = None
        model = _models(
            fit_regional_mutation_rate_models(_prepared(rows, {"AAA": ["a", "b"]}), 0.3)
        )["AAA"]
        assert model.n_mutated + model.n_background == len(rows) - 15
        assert model.n_mutated == sum(r["y"] for r in rows[15:])

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
        ht = fit_regional_mutation_rate_models(
            _prepared(rows, {"AAA": ["a", "b"]}), 0.3
        )
        assert ht.count() == 0

    def test_not_converged(self, caplog: pytest.LogCaptureFixture) -> None:
        """Test that a fit stopped at ``max_iter`` is flagged and logged."""
        with caplog.at_level(logging.WARNING):
            model = _models(
                fit_regional_mutation_rate_models(
                    _prepared(_training_rows(), {"AAA": ["a", "b"]}), 0.3, max_iter=1
                )
            )["AAA"]
        assert not model.converged
        assert any("AAA did not converge" in m for m in _warnings(caplog))

    def test_one_feature_list(self) -> None:
        """
        Test that one list of features fits every context with rows on those features.

        The models equal those from a dict giving each context the same list.
        """
        rows = _training_rows(seed=1) + _training_rows(seed=2, context="ACG")
        shared = _models(
            fit_regional_mutation_rate_models(_prepared(rows, ["a", "b"]), 0.3)
        )
        per_context = _models(
            fit_regional_mutation_rate_models(
                _prepared(rows, {"AAA": ["a", "b"], "ACG": ["a", "b"]}), 0.3
            )
        )
        assert sorted(shared) == sorted(per_context) == ["AAA", "ACG"]
        for context, model in shared.items():
            assert model.features == per_context[context].features
            for field in ["mean", "sd", "coefs", "intercept"]:
                assert model[field] == pytest.approx(per_context[context][field])

    def test_models_per_context(self) -> None:
        """
        Test that contexts fit together each get the model of their own rows.

        The contexts have different numbers and orders of features. A feature
        constant in one context is left out of that context's model only, a context
        with no rows gets no model, and rows of a context without features are
        ignored.
        """
        aaa = _training_rows(seed=1)
        acg = [
            {**r, "c": r["a"] * r["b"]} for r in _training_rows(seed=2, context="ACG")
        ]
        ggg = _training_rows(seed=3, context="GGG")
        features = {"AAA": ["a", "b", "c"], "ACG": ["c", "b", "a"], "CCC": ["a"]}
        models = _models(
            fit_regional_mutation_rate_models(_prepared(aaa + acg + ggg, features), 0.3)
        )
        assert sorted(models) == ["AAA", "ACG"]
        assert models["AAA"].features == ["a", "b"]
        assert models["ACG"].features == ["c", "b", "a"]
        for context, rows in [("AAA", aaa), ("ACG", acg)]:
            alone = _models(
                fit_regional_mutation_rate_models(
                    _prepared(rows, {context: features[context]}), 0.3
                )
            )[context]
            for field in ["mean", "sd", "coefs", "intercept"]:
                assert models[context][field] == pytest.approx(alone[field])


class TestNewtonUpdate:
    """Test the step halving in _newton_update."""

    @staticmethod
    def _model() -> dict:
        """Return a model one step of [0, 0.2] past its last accepted point."""
        return {
            "intercept": 0.0,
            "coefs": [1.0],
            "penalty": 0.1,
            "n": 1000,
            "objective": 0.5,
            "start": np.array([0.0, 1.2]),
            "step": np.array([0.0, 0.2]),
            "converged": False,
        }

    @staticmethod
    def _sums(objective: float) -> hl.Struct:
        """Return sums whose loss gives the model this objective at its coefficients."""
        # The penalty term at coefficient 1.0 is 0.1 / 2.
        return hl.Struct(
            g=[1.0, -2.0], h=[300.0, 0.0, 0.0, 400.0], loss=1000 * (objective - 0.05)
        )

    def test_rise_halves_the_step(self) -> None:
        """Test that a step that raised the objective is halved from its start."""
        model = self._model()
        _newton_update(model, self._sums(0.6), 1e-8)
        assert list(model["step"]) == [0.0, 0.1]
        assert [model["intercept"], *model["coefs"]] == pytest.approx([0.0, 1.1])
        assert not model["converged"]

    def test_rise_within_rounding_is_kept(self) -> None:
        """
        Test that a rise within the rounding of the summed loss isn't an overshoot.

        Near the minimum the rounding can exceed a step's real change. Treating it as
        an overshoot halved the step over and over on real data, so fits that full
        Newton steps converged on never converged.
        """
        model = self._model()
        _newton_update(model, self._sums(0.5 + 1e-15), 1e-8)
        assert list(model["start"]) == [0.0, 1.0]
        assert model["objective"] == pytest.approx(0.5 + 1e-15)

    def test_step_halved_below_tol_has_converged(self) -> None:
        """Test that a model whose halved step is below ``tol`` has converged."""
        model = {**self._model(), "step": np.array([0.0, 1.5e-8])}
        _newton_update(model, self._sums(0.6), 1e-8)
        assert model["converged"]


class TestTuneRegionalMutationRatePenalty:
    """Test the tune_regional_mutation_rate_penalty function."""

    @staticmethod
    def _expected_loss(
        rows: List[Dict], fold: str, penalty: float, features: List[str]
    ) -> float:
        """Return a fold's summed log loss under a model fit on the other rows."""
        rows = [r for r in rows if all(r[f] is not None for f in ["y", *features])]
        fit_rows = [r for r in rows if r["fold"] != fold]
        model = _models(
            fit_regional_mutation_rate_models(
                _prepared(fit_rows, {rows[0]["context"]: features}), penalty
            )
        )[rows[0]["context"]]
        y, x = _arrays([r for r in rows if r["fold"] == fold], features)
        return _log_loss(y, _predict(model, x))

    @pytest.mark.parametrize("features", [["a", "b"], {"AAA": ["a", "b"]}])
    def test_held_out_log_loss(
        self, features: Union[List[str], Dict[str, List[str]]]
    ) -> None:
        """
        Test that each fold is scored, at each penalty, by a model of the other folds.

        Rows with a missing fold are used for fitting but not scored, and rows
        missing the label or a feature are neither. The features can be one list or
        a dict per context.
        """
        rows = _training_rows(seed=2)
        for r in rows[:200]:
            r["fold"] = None
        for r in rows[200:220]:
            r["a"] = None
        for r in rows[220:230]:
            r["y"] = None
        penalties = [0.3, 30.0]
        tuning = {
            (r.penalty, r.fold): r
            for r in tune_regional_mutation_rate_penalty(
                _prepared(rows, features, fold=True), penalties
            ).collect()
        }
        assert sorted(tuning) == [(p, f) for p in penalties for f in ["even", "odd"]]
        assert sum(r.n_test for r in tuning.values()) == 2 * (len(rows) - 230)
        for (penalty, fold), r in tuning.items():
            assert r.log_loss == pytest.approx(
                self._expected_loss(rows, fold, penalty, ["a", "b"]), rel=1e-6
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
                prepare_regional_mutation_rate_training_ht(
                    ht,
                    {"AAA": ["a", "b"], "ACG": ["a", "b"]},
                    label_expr=ht.label,
                    context_expr=ht.ctx,
                    fold_expr=ht.split,
                ),
                [0.3],
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
            _prepared(rows, {"AAA": ["a", "b"]}, fold=True), [0.01, 1000.0]
        )
        assert ht.count() == 0
        assert list(ht.row) == ["penalty", "context", "fold", "n_test", "log_loss"]

    def test_repeated_penalty_is_tried_once(self) -> None:
        """
        Test that a penalty listed twice gets one row per context and fold.

        Two rows would double its summed loss when the penalty is chosen.
        """
        ht = tune_regional_mutation_rate_penalty(
            _prepared(_training_rows(seed=2), {"AAA": ["a", "b"]}, fold=True),
            [0.3, 0.3],
        )
        assert ht.aggregate(hl.agg.counter(ht.penalty)) == {0.3: 2}

    def test_too_strong_penalty_loses(self) -> None:
        """Test that a penalty that flattens a real effect has a higher held-out loss."""
        ht = tune_regional_mutation_rate_penalty(
            _prepared(_training_rows(seed=2), {"AAA": ["a", "b"]}, fold=True),
            [0.01, 1000.0],
        )
        loss = {(r.penalty, r.fold): r.log_loss for r in ht.collect()}
        for fold in ["even", "odd"]:
            assert loss[(0.01, fold)] < loss[(1000.0, fold)]


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

    def test_lowest_summed_loss(self, caplog: pytest.LogCaptureFixture) -> None:
        """Test that the penalty with the lowest loss summed over contexts is chosen."""
        with caplog.at_level(logging.WARNING):
            assert choose_regional_mutation_rate_penalty(self._tuning_ht(), 1.0) == 0.3
        assert not _warnings(caplog)

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
        assert any("at the end of those tried" in m for m in _warnings(caplog))

    def test_tie_goes_to_the_smaller_penalty(self) -> None:
        """Test that of penalties with the same summed loss, the smaller is chosen."""
        assert (
            choose_regional_mutation_rate_penalty(self._tuning_ht([0.1, 1.0]), 0.3)
            == 0.1
        )

    def test_one_penalty(self, caplog: pytest.LogCaptureFixture) -> None:
        """Test that a single penalty tried is returned without a warning."""
        with caplog.at_level(logging.WARNING):
            assert (
                choose_regional_mutation_rate_penalty(self._tuning_ht([0.3]), 1.0)
                == 0.3
            )
        assert not _warnings(caplog)

    def test_empty(self, caplog: pytest.LogCaptureFixture) -> None:
        """Test that the default penalty is returned, with a warning, when empty."""
        with caplog.at_level(logging.WARNING):
            assert (
                choose_regional_mutation_rate_penalty(self._tuning_ht([]), 0.3) == 0.3
            )
        assert any("No held-out rows" in m for m in _warnings(caplog))


def _features_ht(rows: list) -> hl.Table:
    """Return a feature Table keyed by window with float features ``a`` and ``b``."""
    return hl.Table.parallelize(
        [{"interval": _window(s, e), "a": a, "b": b} for s, e, a, b in rows],
        hl.tstruct(interval=INTERVAL_TYPE, a=hl.tfloat64, b=hl.tfloat64),
        key="interval",
    )


def _models_ht() -> hl.Table:
    """
    Return two one-feature models with hand-computable rates.

    The AAA model's rate is ``expit(0.5 + (a - 1)) / 0.6``, and the ACG model's is
    ``expit(-b / 2) / 0.5``.
    """
    common = {"n_mutated": 10, "n_background": 10, "converged": True}
    return hl.Table.parallelize(
        [
            {
                "context": "AAA",
                "features": ["a"],
                "mean": [1.0],
                "sd": [1.0],
                "coefs": [1.0],
                "intercept": 0.5,
                "ave": 0.6,
                **common,
            },
            {
                "context": "ACG",
                "features": ["b"],
                "mean": [0.0],
                "sd": [2.0],
                "coefs": [-1.0],
                "intercept": 0.0,
                "ave": 0.5,
                **common,
            },
        ],
        MODELS_HT_SCHEMA,
        key="context",
    )


def _rates(
    features_ht: hl.Table,
    models_ht: hl.Table,
    contexts: List[str],
    missing_feature_rate: Optional[float] = 1.0,
) -> Dict[int, Optional[Dict[str, float]]]:
    """Return each window's rates by its start, from regional_mutation_rate_expr."""
    ht = features_ht.annotate(
        adj_r=regional_mutation_rate_expr(
            features_ht.row, models_ht, contexts, missing_feature_rate
        )
    )
    assert ht.adj_r.dtype == hl.tdict(hl.tstr, hl.tfloat64)
    return {r.interval.start.position: r.adj_r for r in ht.collect()}


class TestRegionalMutationRateExpr:
    """Test the regional_mutation_rate_expr function."""

    def test_rates(self) -> None:
        """
        Test the rate per window and the 1.0 fill.

        A window missing a context's feature gets 1.0 for that context, contexts
        without a model get 1.0, and a window with no rate for any modeled context
        has no rates.
        """
        rates = _rates(
            _features_ht(
                [(1, 1001, 0.5, 2.0), (1001, 2001, 1.0, None), (2001, 3001, None, None)]
            ),
            _models_ht(),
            ["AAA", "ACG", "CCC"],
        )
        assert rates[1] == pytest.approx(
            {"AAA": _expit(0.0) / 0.6, "ACG": _expit(-1.0) / 0.5, "CCC": 1.0}
        )
        assert rates[1001] == pytest.approx(
            {"AAA": _expit(0.5) / 0.6, "ACG": 1.0, "CCC": 1.0}
        )
        assert rates[2001] is None

    def test_only_requested_contexts(self) -> None:
        """
        Test that models of contexts that weren't asked for are ignored.

        The second window has an AAA rate but no ACG rate, so with only ACG asked
        for it has no rates rather than ACG filled with 1.0.
        """
        rates = _rates(
            _features_ht([(1, 1001, 0.5, 2.0), (1001, 2001, 1.0, None)]),
            _models_ht(),
            ["ACG"],
        )
        assert rates[1] == pytest.approx({"ACG": _expit(-1.0) / 0.5})
        assert rates[1001] is None

    def test_missing_feature_rate(self) -> None:
        """
        Test that a context missing a model feature gets ``missing_feature_rate``.

        With None, its rate is left missing, while a context without a model still
        gets 1.0.
        """
        features_ht = _features_ht([(1, 1001, 1.0, None)])
        assert _rates(features_ht, _models_ht(), ["ACG", "AAA", "CCC"])[1] == (
            pytest.approx({"AAA": _expit(0.5) / 0.6, "ACG": 1.0, "CCC": 1.0})
        )
        rates = _rates(features_ht, _models_ht(), ["ACG", "AAA", "CCC"], None)[1]
        assert rates["ACG"] is None
        assert rates["AAA"] == pytest.approx(_expit(0.5) / 0.6)
        assert rates["CCC"] == 1.0

    def test_no_models(self) -> None:
        """Test that with no models every window's rates are missing."""
        rates = _rates(
            _features_ht([(1, 1001, 0.0, 2.0)]), _models_ht().filter(False), ["AAA"]
        )
        assert rates == {1: None}

    def test_matches_the_fitted_model(self) -> None:
        """
        Test that the rate is the fitted model's prediction divided by ``ave``.

        The two features have different means and sds, so the test covers the mean
        subtraction, each feature's sd and coefficient, and the feature order. A
        window at the training means has rate 1, and a window missing a feature has
        no rate.
        """
        models_ht = fit_regional_mutation_rate_models(
            _prepared(_training_rows(), {"AAA": ["a", "b"]}), 0.3
        )
        model = _models(models_ht)["AAA"]
        values = [(-1.0, 9.0), (0.5, 5.0), (2.0, 1.0), tuple(model.mean), (None, 1.0)]
        rates = _rates(
            _features_ht(
                [
                    (1 + 1000 * i, 1001 + 1000 * i, a, b)
                    for i, (a, b) in enumerate(values)
                ]
            ),
            models_ht,
            ["AAA"],
        )
        expected = _predict(model, np.array(values[:4])) / model.ave
        assert [rates[1 + 1000 * i]["AAA"] for i in range(4)] == pytest.approx(
            list(expected)
        )
        assert rates[3001]["AAA"] == pytest.approx(1.0)
        assert rates[4001] is None
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
        assert rates[100]["ACG"] == pytest.approx(2.0)
        assert rates[1]["CCC"] == 3.0
