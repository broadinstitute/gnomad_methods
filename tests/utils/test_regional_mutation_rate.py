"""Tests for the regional mutation rate utility module."""

import logging
from typing import Dict, List, Optional, Union

import hail as hl
import numpy as np
import pytest
import statsmodels.api as sm

from gnomad.utils.constraint import (
    MODELS_HT_SCHEMA,
    choose_regional_mutation_rate_penalty,
    fit_regional_mutation_rate_models,
    normalize_regional_mutation_rate_expr,
    prepare_regional_mutation_rate_training_ht,
    regional_mutation_rate_expr,
    regional_mutation_rate_normalization_factors_agg_expr,
    screen_regional_mutation_rate_features,
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


def _lognormal_rows(x: np.ndarray, y: np.ndarray) -> List[Dict]:
    """Return AAA training rows with labels ``y`` and features ``a``, ``b``, ``c``."""
    return [
        {
            "y": int(y[i]),
            "a": float(x[i, 0]),
            "b": float(x[i, 1]),
            "c": float(x[i, 2]),
            "context": "AAA",
            "fold": None,
        }
        for i in range(len(y))
    ]


def _heavy_tailed_rows() -> List[Dict]:
    """
    Return 200 rows of three lognormal features, labeled by a logistic model of ``a``.

    The features have an sd of 2 on the log scale, so a few standardized values are
    far out. The second full Newton step from the start overshoots on these rows.
    """
    rng = np.random.default_rng(52)
    log_x = rng.normal(size=(200, 3))
    y = rng.random(200) < _expit(np.log(0.1 / 0.9) + 0.5 * log_x[:, 0])
    return _lognormal_rows(np.exp(2.0 * log_x), y)


def _rounding_rows() -> List[Dict]:
    """
    Return 50 rows of three lognormal features, with 2 of them mutated.

    The features have an sd of 2.5 on the log scale, so a few standardized values are
    far out. The log loss of these rows has a minimum, but at it 4 rows' fitted
    probabilities round to 0. At a penalty of 0.01, full Newton steps from the start
    overshoot on these rows and don't converge.
    """
    rng = np.random.default_rng(23)
    x = np.exp(2.5 * rng.normal(size=(50, 3)))
    y = rng.random(50) < 0.05
    return _lognormal_rows(x, y)


def _separated_rows() -> List[Dict]:
    """
    Return training rows where ``c`` is 1 on 5% of the non-mutated rows and 0 elsewhere.

    ``c`` separates the labels in part: the log loss falls without bound as its
    coefficient goes to minus infinity.
    """
    rng = np.random.default_rng(3)
    return [
        {**r, "c": float(r["y"] == 0 and rng.random() < 0.05)} for r in _training_rows()
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

    def test_fold_must_be_a_string(self) -> None:
        """
        Test that a fold that isn't a string raises.

        A boolean or integer fold used to be kept, and the tuning then failed with an
        unrelated type error after two passes over the Table.
        """
        ht = _training_ht(_training_rows())
        with pytest.raises(ValueError, match="fold_expr must be a string"):
            prepare_regional_mutation_rate_training_ht(
                ht, ["a", "b"], fold_expr=ht.fold == "odd"
            )

    def test_fit_and_tune_need_a_prepared_table(self) -> None:
        """Test that the fit and the tuning say how to prepare a Table they can't use."""
        with pytest.raises(ValueError, match="prepare_regional_mutation_rate"):
            fit_regional_mutation_rate_models(_training_ht(_training_rows()))
        with pytest.raises(ValueError, match="no fold.*with a fold_expr"):
            tune_regional_mutation_rate_penalty(
                _prepared(_training_rows(), {"AAA": ["a", "b"]}), [0.3]
            )


class TestFitRegionalMutationRateModels:
    """Test the fit_regional_mutation_rate_models function."""

    @staticmethod
    def _gradient(
        rows: List[Dict], model: hl.Struct, penalty: float = 0.0
    ) -> np.ndarray:
        """Return the gradient of the penalized mean log loss at a model, with numpy."""
        y, x = _arrays(rows, model.features)
        z = (x - np.array(model.mean)) / np.array(model.sd)
        residual = _predict(model, x) - y
        coefs = np.array(model.coefs)
        return np.r_[residual.mean(), z.T @ residual / len(y) + penalty * coefs]

    @pytest.mark.parametrize("penalty", [None, 0.01, 0.3])
    def test_minimizes_the_penalized_log_loss(self, penalty: Optional[float]) -> None:
        """
        Test that the fit minimizes mean log loss + penalty / 2 * sum(coef ** 2).

        The gradient of that objective, with the intercept unpenalized and the
        features standardized with the stored ``mean`` and ``sd`` (ddof=0), is 0 at
        the fit. With no penalty given, there is none. The models Table records the
        penalty.
        """
        rows = _training_rows()
        ht = _prepared(rows, {"AAA": ["a", "b"]})
        models_ht = (
            fit_regional_mutation_rate_models(ht)
            if penalty is None
            else fit_regional_mutation_rate_models(ht, penalty)
        )
        assert hl.eval(models_ht.ridge_penalty) == (penalty or 0.0)
        model = _models(models_ht)["AAA"]
        assert self._gradient(rows, model, penalty or 0.0) == pytest.approx(
            [0, 0, 0], abs=1e-8
        )
        _, x = _arrays(rows, ["a", "b"])
        assert model.mean == pytest.approx(list(x.mean(axis=0)))
        assert model.sd == pytest.approx(list(x.std(axis=0)))
        assert model.converged
        assert model.ave == pytest.approx(_expit(model.intercept))

    def test_step_halving(self) -> None:
        """
        Test that the fit reaches the minimum where full Newton steps overshoot.

        On these heavy-tailed rows, a full step raises the loss, so it has to be
        halved.
        """
        rows = _heavy_tailed_rows()
        model = _models(
            fit_regional_mutation_rate_models(_prepared(rows, {"AAA": ["a", "b", "c"]}))
        )["AAA"]
        assert model.converged
        assert self._gradient(rows, model) == pytest.approx([0, 0, 0, 0], abs=1e-8)

    def test_step_halving_with_a_penalty(self) -> None:
        """
        Test that a penalized fit reaches the minimum where full Newton steps don't.

        Full steps overshoot on these rows at a penalty of 0.01 and don't converge,
        so a step that raises the objective has to be halved.
        """
        rows = _rounding_rows()
        model = _models(
            fit_regional_mutation_rate_models(
                _prepared(rows, {"AAA": ["a", "b", "c"]}), 0.01
            )
        )["AAA"]
        assert model.converged
        assert self._gradient(rows, model, 0.01) == pytest.approx(
            [0, 0, 0, 0], abs=1e-8
        )

    def test_penalized_fit_can_round_probabilities(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        """
        Test that a penalized fit isn't logged for fitted probabilities of 0 or 1.

        At a penalty of 1e-5 a fitted probability of these rows still rounds to 0 at
        the minimum. The penalty gives the objective a minimum, so unlike an
        unpenalized fit, the rounding can't mean the features separate the labels.
        """
        rows = _rounding_rows()
        with caplog.at_level(logging.WARNING):
            model = _models(
                fit_regional_mutation_rate_models(
                    _prepared(rows, {"AAA": ["a", "b", "c"]}), 1e-5
                )
            )["AAA"]
        assert model.converged
        assert self._gradient(rows, model, 1e-5) == pytest.approx(
            [0, 0, 0, 0], abs=1e-8
        )
        assert not any("fitted probability" in m for m in _warnings(caplog))

    @pytest.mark.parametrize("penalty", [-1.0, float("nan"), float("inf")])
    def test_fit_and_tune_reject_bad_penalties(self, penalty: float) -> None:
        """
        Test that the fit and the tuning reject a negative, NaN, or infinite penalty.

        A negative penalty rewards large coefficients, and a NaN penalty makes the
        models NaN. An infinite one passed the check and then failed in the
        least-squares solve after the training passes. Tuning with no penalties
        raises too.
        """
        ht = _prepared(_training_rows(), {"AAA": ["a", "b"]}, fold=True)
        with pytest.raises(ValueError, match="finite and at least 0"):
            fit_regional_mutation_rate_models(ht, penalty)
        with pytest.raises(ValueError, match="finite and at least 0"):
            tune_regional_mutation_rate_penalty(ht, [0.3, penalty])
        with pytest.raises(ValueError, match="No ridge penalties"):
            tune_regional_mutation_rate_penalty(ht, [])

    @pytest.mark.parametrize("value", [float("nan"), float("inf")])
    def test_non_finite_features_raise(self, value: float) -> None:
        """
        Test that a NaN or infinite feature value makes the fit and the tuning raise.

        NaN counts as defined in Hail, so such a row survived the missing-value
        filter and made its context's model NaN.
        """
        rows = _training_rows()
        rows[3]["a"] = value
        ht = _prepared(rows, {"AAA": ["a", "b"]}, fold=True)
        with pytest.raises(ValueError, match="1 training rows have a NaN or infinite"):
            fit_regional_mutation_rate_models(ht)
        with pytest.raises(ValueError, match="NaN or infinite"):
            tune_regional_mutation_rate_penalty(ht, [0.3])

    @pytest.mark.parametrize("label", [2.0, float("nan")])
    def test_labels_must_be_0_or_1(self, label: float) -> None:
        """
        Test that a label other than 0 or 1, including NaN, makes the fit raise.

        A NaN label counts as defined in Hail, and crashed the fit with an unrelated
        error; a label of 2 was counted as mutated.
        """
        ht = _training_ht(_training_rows())
        ht = prepare_regional_mutation_rate_training_ht(
            ht,
            {"AAA": ["a", "b"]},
            label_expr=hl.if_else(ht.a > 2, label, hl.float64(ht.y)),
        )
        with pytest.raises(ValueError, match="have a label other than 0 or 1"):
            fit_regional_mutation_rate_models(ht)

    @pytest.mark.parametrize("max_iter", [2, 3])
    def test_stopped_early_keeps_the_last_accepted_point(self, max_iter: int) -> None:
        """
        Test that a fit that runs out of steps keeps its last accepted point.

        On the heavy-tailed rows the first Newton step is accepted and the second
        overshoots. Stopped after the second step is taken (2 passes) or while it is
        halved (3 passes), the fit returns the coefficients it last accepted, one
        Newton step from the start, rather than a point it never evaluated.
        """
        rows = _heavy_tailed_rows()
        model = _models(
            fit_regional_mutation_rate_models(
                _prepared(rows, {"AAA": ["a", "b", "c"]}), max_iter=max_iter
            )
        )["AAA"]
        assert not model.converged
        # One Newton step from the start: the label rate's intercept, coefficients 0.
        # solve(H, g) on the summed terms is the module's lstsq step on the means.
        y, x = _arrays(rows, ["a", "b", "c"])
        v = np.c_[np.ones(len(y)), (x - x.mean(axis=0)) / x.std(axis=0)]
        start = np.r_[np.log(y.sum() / (len(y) - y.sum())), np.zeros(3)]
        p = _expit(v @ start)
        step = np.linalg.solve((v.T * p * (1 - p)) @ v, v.T @ (y - p))
        assert [model.intercept, *model.coefs] == pytest.approx(list(start + step))

    def test_collinear_features(self) -> None:
        """
        Test that exactly collinear features still fit, with the same predictions.

        ``c`` is a linear function of ``a``, so the Hessian is singular. The
        least-squares steps split ``a``'s effect between ``a`` and ``c``, and predict
        what the model without ``c`` does.
        """
        rows = [{**r, "c": 2 * r["a"] + 1} for r in _training_rows()]
        model = _models(
            fit_regional_mutation_rate_models(_prepared(rows, {"AAA": ["a", "b", "c"]}))
        )["AAA"]
        without_c = _models(
            fit_regional_mutation_rate_models(_prepared(rows, {"AAA": ["a", "b"]}))
        )["AAA"]
        _, x = _arrays(rows, ["a", "b", "c"])
        assert model.converged
        assert list(_predict(model, x)) == pytest.approx(
            list(_predict(without_c, x[:, :2]))
        )

    def test_separated_labels(self, caplog: pytest.LogCaptureFixture) -> None:
        """
        Test that a fit whose features separate the labels is flagged and logged.

        The log loss of these rows has no minimum. As ``c``'s coefficient grows, the
        curvature along it vanishes, the least-squares steps stop moving it, and the
        fit looked converged at an arbitrary coefficient.
        """
        with caplog.at_level(logging.WARNING):
            model = _models(
                fit_regional_mutation_rate_models(
                    _prepared(_separated_rows(), {"AAA": ["a", "b", "c"]})
                )
            )["AAA"]
        assert not model.converged
        assert any(
            "AAA did not converge" in m and "separate the labels" in m
            for m in _warnings(caplog)
        )

    @pytest.mark.parametrize("case", ["separated", "collinear"])
    def test_penalized_fit_has_a_minimum(self, case: str) -> None:
        """
        Test that a penalty fits rows on which an unpenalized fit has no unique minimum.

        The penalty adds curvature along every coefficient, so the Hessian keeps its
        rank: the fit of the separated rows converges instead of being marked
        separated, and collinear features get unique coefficients.
        """
        rows = (
            _separated_rows()
            if case == "separated"
            else [{**r, "c": 2 * r["a"] + 1} for r in _training_rows()]
        )
        model = _models(
            fit_regional_mutation_rate_models(
                _prepared(rows, {"AAA": ["a", "b", "c"]}), 0.3
            )
        )["AAA"]
        assert model.converged
        assert self._gradient(rows, model, 0.3) == pytest.approx([0, 0, 0, 0], abs=1e-8)

    def test_rounded_probabilities_are_logged(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        """
        Test that fitted probabilities of 0 or 1 at a real minimum are only logged.

        At these rows' minimum (the gradient is 0), 4 fitted probabilities round to 0.
        That marked the fit not converged, as if its features separated the labels;
        R's glm.fit only warns.
        """
        rows = _rounding_rows()
        with caplog.at_level(logging.WARNING):
            model = _models(
                fit_regional_mutation_rate_models(
                    _prepared(rows, {"AAA": ["a", "b", "c"]})
                )
            )["AAA"]
        assert model.converged
        assert self._gradient(rows, model) == pytest.approx([0, 0, 0, 0], abs=1e-8)
        assert any(
            "AAA: 4 training rows have a fitted probability numerically 0 or 1" in m
            for m in _warnings(caplog)
        )

    def test_counts_and_expressions(self) -> None:
        """Test the label counts, and labels and contexts from other expressions."""
        rows = _training_rows()
        ht = _training_ht(rows).rename({"y": "label", "context": "ctx"})
        model = _models(
            fit_regional_mutation_rate_models(
                prepare_regional_mutation_rate_training_ht(
                    ht, {"AAA": ["a", "b"]}, label_expr=ht.label, context_expr=ht.ctx
                )
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
            fit_regional_mutation_rate_models(_prepared(rows, {"AAA": ["a", "c", "b"]}))
        )["AAA"]
        expected = _models(
            fit_regional_mutation_rate_models(_prepared(rows, {"AAA": ["a", "b"]}))
        )["AAA"]
        assert model.features == ["a", "b"]
        for field in ["mean", "sd", "coefs"]:
            assert model[field] == pytest.approx(expected[field])

    def test_nearly_constant_feature(self) -> None:
        """
        Test that a feature that varies little next to its mean gets its exact sd.

        An sd from the sum of squares, in the same pass as the mean, is NaN for these
        values, which made the fit fail for every context. Half the rows have the
        larger value, since without a penalty a feature that sets one row apart
        separates its label.
        """
        rows = [
            {**r, "c": 0.4123 + 1e-9 * (i % 2)} for i, r in enumerate(_training_rows())
        ]
        model = _models(
            fit_regional_mutation_rate_models(_prepared(rows, {"AAA": ["a", "b", "c"]}))
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
            fit_regional_mutation_rate_models(_prepared(rows, {"AAA": ["a", "b"]}))
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
        ht = fit_regional_mutation_rate_models(_prepared(rows, {"AAA": ["a", "b"]}))
        assert ht.count() == 0

    def test_not_converged(self, caplog: pytest.LogCaptureFixture) -> None:
        """Test that a fit stopped at ``max_iter`` is flagged and logged."""
        with caplog.at_level(logging.WARNING):
            model = _models(
                fit_regional_mutation_rate_models(
                    _prepared(_training_rows(), {"AAA": ["a", "b"]}), max_iter=1
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
        shared = _models(fit_regional_mutation_rate_models(_prepared(rows, ["a", "b"])))
        per_context = _models(
            fit_regional_mutation_rate_models(
                _prepared(rows, {"AAA": ["a", "b"], "ACG": ["a", "b"]})
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
            fit_regional_mutation_rate_models(_prepared(aaa + acg + ggg, features))
        )
        assert sorted(models) == ["AAA", "ACG"]
        assert models["AAA"].features == ["a", "b"]
        assert models["ACG"].features == ["c", "b", "a"]
        for context, rows in [("AAA", aaa), ("ACG", acg)]:
            alone = _models(
                fit_regional_mutation_rate_models(
                    _prepared(rows, {context: features[context]})
                )
            )[context]
            for field in ["mean", "sd", "coefs", "intercept"]:
                assert models[context][field] == pytest.approx(alone[field])


class TestScreenRegionalMutationRateFeatures:
    """Test the screen_regional_mutation_rate_features function."""

    @staticmethod
    def _screen(ht: hl.Table) -> Dict[tuple, hl.Struct]:
        """Return the screen of a prepared Table by context and feature."""
        return {
            (r.context, r.feature): r
            for r in screen_regional_mutation_rate_features(ht).collect()
        }

    @staticmethod
    def _check_statsmodels(r: hl.Struct, rows: List[Dict], feature: str) -> None:
        """
        Check a screen row against statsmodels' fit of the label on the feature.

        The fit is the univariate unpenalized logistic fit of the label on the
        feature standardized with ddof=0, the screen Gnocchi ran with statsmodels.
        Gnocchi used ``fit_regularized`` with no penalty, which gives the same fit
        up to its tolerance. The p-value is compared on the log scale, since the
        p-values of the real features are around 1e-40.
        """
        y, x = _arrays(rows, [feature])
        z = (x - x.mean()) / x.std()
        fit = sm.Logit(y, sm.add_constant(z)).fit(disp=False)
        assert r.converged
        assert r.coef == pytest.approx(fit.params[1], rel=1e-6)
        assert r.se == pytest.approx(fit.bse[1], rel=1e-6)
        assert r.coef / r.se == pytest.approx(fit.tvalues[1], rel=1e-6)
        assert np.log(r.p_value) == pytest.approx(np.log(fit.pvalues[1]), rel=1e-6)
        assert (r.n, r.n_mutated) == (len(y), int(y.sum()))

    def test_matches_statsmodels(self) -> None:
        """
        Test each feature's coefficient, standard error, and p-value against statsmodels.

        ``c`` is constant, so it has no row.
        """
        rows = _training_rows()
        screen = self._screen(_prepared(rows, {"AAA": ["a", "b", "c"]}))
        assert sorted(screen) == [("AAA", "a"), ("AAA", "b")]
        for feature in ["a", "b"]:
            self._check_statsmodels(screen[("AAA", feature)], rows, feature)

    def test_noise_feature(self) -> None:
        """
        Test a feature with no association with the label against statsmodels.

        Its p-value is far above any selection cutoff.
        """
        rng = np.random.default_rng(5)
        rows = [{**r, "a": float(rng.normal())} for r in _training_rows()]
        r = self._screen(_prepared(rows, {"AAA": ["a"]}))[("AAA", "a")]
        self._check_statsmodels(r, rows, "a")
        assert r.p_value > 0.05

    def test_rows_missing_another_feature_count(self) -> None:
        """
        Test that each feature is fit on the rows where it is defined.

        A row missing ``a`` still counts for ``b``, as each feature was screened on
        its own rows in Gnocchi, so ``b``'s screen matches the fit on all the rows.
        """
        rows = _training_rows()
        for r in rows[:100]:
            r["a"] = None
        screen = self._screen(_prepared(rows, {"AAA": ["a", "b"]}))
        assert screen[("AAA", "a")].n == len(rows) - 100
        assert screen[("AAA", "a")].n_mutated == sum(r["y"] for r in rows[100:])
        assert screen[("AAA", "b")].n == len(rows)
        all_rows = self._screen(_prepared(_training_rows(), {"AAA": ["b"]}))
        assert screen[("AAA", "b")].p_value == pytest.approx(
            all_rows[("AAA", "b")].p_value
        )

    def test_separating_feature_has_no_p_value(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        """
        Test that a fit that doesn't converge has no se or p-value, and is logged.

        ``c`` separates the labels in part, so its log loss has no minimum. The
        other features' fits are unaffected.
        """
        with caplog.at_level(logging.WARNING):
            screen = self._screen(_prepared(_separated_rows(), {"AAA": ["a", "c"]}))
        c = screen[("AAA", "c")]
        assert not c.converged
        assert (c.se, c.p_value) == (None, None)
        assert screen[("AAA", "a")].converged
        assert screen[("AAA", "a")].p_value < 1e-10
        assert any(
            "1 of 2 screen fits" in m and "AAA|c" in m for m in _warnings(caplog)
        )

    def test_max_iter_and_tol(self, caplog: pytest.LogCaptureFixture) -> None:
        """
        Test that ``max_iter`` and ``tol`` reach the fits.

        One Newton step from 0 doesn't converge at the default tolerance, so the fit
        has no p-value, and is logged. With a tolerance larger than that step, the
        same step converges, short of the full fit.
        """
        ht = _prepared(_training_rows(), {"AAA": ["a"]})
        with caplog.at_level(logging.WARNING):
            one_step = screen_regional_mutation_rate_features(ht, max_iter=1).collect()
        assert not one_step[0].converged
        assert (one_step[0].se, one_step[0].p_value) == (None, None)
        assert any("1 of 1 screen fits" in m for m in _warnings(caplog))
        loose = screen_regional_mutation_rate_features(ht, max_iter=1, tol=10.0)
        full = screen_regional_mutation_rate_features(ht)
        loose, full = loose.collect()[0], full.collect()[0]
        assert loose.converged
        assert loose.coef != pytest.approx(full.coef, rel=1e-3)

    def test_pairs_with_no_fit(self, caplog: pytest.LogCaptureFixture) -> None:
        """
        Test that (context, feature) pairs with no fit have no row, and are logged.

        ``a`` is missing in every row, ``c`` is constant, and GGG has no rows. With no
        fit at all, the screen is an empty Table that still records its features.
        """
        rows = [{**r, "a": None} for r in _training_rows()]
        with caplog.at_level(logging.WARNING):
            screen = self._screen(
                _prepared(rows, {"AAA": ["a", "b", "c"], "GGG": ["b"]})
            )
        assert sorted(screen) == [("AAA", "b")]
        assert any(
            "3 of 4 (context, feature) pairs" in m and "AAA|a, AAA|c, GGG|b" in m
            for m in _warnings(caplog)
        )
        empty = screen_regional_mutation_rate_features(
            _prepared(rows, {"AAA": ["a", "c"]})
        )
        assert empty.collect() == []
        assert hl.eval(empty.features) == {"AAA": ["a", "c"]}

    def test_contexts(self) -> None:
        """
        Test that each context is screened on its own features and rows.

        CCC's rows all have label 0, so it has no rows. With one feature list,
        every context the rows have is screened on it. The ``features`` global
        records what was screened, so a selection from a screen of other features
        can be caught.
        """
        rows = (
            _training_rows(context="AAA")
            + _training_rows(seed=1, context="ACG")
            + [{**r, "y": 0, "context": "CCC"} for r in _training_rows(n=50, seed=2)]
        )
        features = {"AAA": ["a"], "ACG": ["b"], "CCC": ["a"]}
        per_context_ht = screen_regional_mutation_rate_features(
            _prepared(rows, features)
        )
        assert hl.eval(per_context_ht.features) == features
        per_context = {(r.context, r.feature): r for r in per_context_ht.collect()}
        assert sorted(per_context) == [("AAA", "a"), ("ACG", "b")]
        aaa_only = self._screen(_prepared(rows[:2000], {"AAA": ["a"]}))
        assert per_context[("AAA", "a")].p_value == pytest.approx(
            aaa_only[("AAA", "a")].p_value
        )
        one_list_ht = screen_regional_mutation_rate_features(_prepared(rows, ["a"]))
        assert hl.eval(one_list_ht.features) == {
            "AAA": ["a"],
            "ACG": ["a"],
            "CCC": ["a"],
        }
        one_list = {(r.context, r.feature): r for r in one_list_ht.collect()}
        assert sorted(one_list) == [("AAA", "a"), ("ACG", "a")]

    def test_bad_inputs_raise(self) -> None:
        """
        Test that an unprepared Table, repeated features, NaN, and bad labels raise.

        A feature listed twice for a context would have each of its rows fit twice,
        shrinking its standard error by a factor of sqrt(2). A bad row is counted
        once, not once per feature.
        """
        with pytest.raises(ValueError, match="prepare_regional_mutation_rate"):
            screen_regional_mutation_rate_features(_training_ht(_training_rows()))
        with pytest.raises(ValueError, match=r"more than once.*'AAA': \['a'\]"):
            screen_regional_mutation_rate_features(
                _prepared(_training_rows(), {"AAA": ["a", "b", "a"]})
            )
        rows = _training_rows()
        rows[3]["a"] = float("nan")
        rows[3]["b"] = float("inf")
        with pytest.raises(ValueError, match="^1 training rows have a NaN or infinite"):
            screen_regional_mutation_rate_features(_prepared(rows, {"AAA": ["a", "b"]}))
        ht = _training_ht(_training_rows())
        ht = prepare_regional_mutation_rate_training_ht(
            ht,
            {"AAA": ["a", "b"]},
            label_expr=hl.if_else(ht.a > 2, 2.0, hl.float64(ht.y)),
        )
        n_bad = sum(r["a"] > 2 for r in _training_rows())
        with pytest.raises(
            ValueError, match=f"^{n_bad} training rows have a label other than 0 or 1"
        ):
            screen_regional_mutation_rate_features(ht)


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
    def _sums(
        model: dict, objective: float, h: Optional[List[float]] = None
    ) -> hl.Struct:
        """Return sums whose loss gives the model this objective at its coefficients."""
        penalty_term = model["penalty"] / 2 * sum(c**2 for c in model["coefs"])
        return hl.Struct(
            g=[1.0, -2.0],
            h=h or [300.0, 0.0, 0.0, 400.0],
            loss=model["n"] * (objective - penalty_term),
            n_saturated=0,
        )

    @pytest.mark.parametrize("objective", [0.6, 0.5 + 1e-10])
    def test_rise_halves_the_step(self, objective: float) -> None:
        """
        Test that a step that raised the objective is halved from its start.

        A rise of 1e-10 is far above the rounding in the loss, so it is an overshoot
        too; that bounds the rounding allowance from above. The start and its
        objective stay those of the last accepted point.
        """
        model = self._model()
        _newton_update(model, self._sums(model, objective), 1e-8)
        assert list(model["step"]) == [0.0, 0.1]
        assert [model["intercept"], *model["coefs"]] == pytest.approx([0.0, 1.1])
        assert list(model["start"]) == [0.0, 1.2]
        assert model["objective"] == 0.5
        assert not model["converged"]

    def test_repeated_rises_keep_halving(self) -> None:
        """
        Test that a halved step that still raises the objective is halved again.

        Each halving is compared with the last accepted point. Comparing with the
        rejected point instead let a step halve only once, and in heavy-tailed
        replicas fits then cycled or failed.
        """
        model = self._model()
        _newton_update(model, self._sums(model, 0.6), 1e-8)
        _newton_update(model, self._sums(model, 0.55), 1e-8)
        assert list(model["step"]) == [0.0, 0.05]
        assert [model["intercept"], *model["coefs"]] == pytest.approx([0.0, 1.15])
        assert list(model["start"]) == [0.0, 1.2]
        assert model["objective"] == 0.5

    def test_rise_within_rounding_is_kept(self) -> None:
        """
        Test that a rise within the rounding of the summed loss isn't an overshoot.

        Near the minimum the rounding can exceed a step's real change. Treating it as
        an overshoot halved the step over and over on real data, so fits that full
        Newton steps converged on never converged. The kept point's Newton step is
        taken from the gradient and Hessian of the penalized mean loss.
        """
        model = self._model()
        _newton_update(model, self._sums(model, 0.5 + 1e-15), 1e-8)
        assert list(model["start"]) == [0.0, 1.0]
        assert model["objective"] == pytest.approx(0.5 + 1e-15, abs=1e-16)
        # Gradient -g / n + penalty * coef and Hessian h / n + penalty, at [0, 1].
        assert list(model["step"]) == pytest.approx([-0.001 / 0.3, 0.102 / 0.5])

    def test_singular_hessian_gives_the_smallest_step(self) -> None:
        """
        Test that a singular Hessian, as from collinear features, still gives a step.

        The step is the least-squares solution with the smallest norm: no change in
        the direction the Hessian has no curvature in. Only an unpenalized Hessian
        can be singular.
        """
        model = {**self._model(), "penalty": 0.0}
        _newton_update(model, self._sums(model, 0.4, h=[300.0, 0.0, 0.0, 0.0]), 1e-8)
        assert list(model["step"]) == pytest.approx([-0.001 / 0.3, 0.0])
        # Singular from the first step, as with collinear features: not separation.
        assert not model["separated"]

    def test_full_step_below_tol_has_converged(self) -> None:
        """Test that a model whose new Newton step is below ``tol`` has converged."""
        model = {**self._model(), "penalty": 0.0}
        sums = hl.Struct(**{**self._sums(model, 0.4), "g": [1e-9, -1e-9]})
        _newton_update(model, sums, 1e-8)
        assert model["converged"]

    def test_halved_step_above_tol_has_not_converged(self) -> None:
        """
        Test that a halved step above ``tol`` leaves the model unconverged.

        The halving keeps the counts of the last accepted point.
        """
        model = {**self._model(), "n_saturated": 3}
        _newton_update(
            model, hl.Struct(**{**self._sums(model, 0.6), "n_saturated": 9}), 1e-8
        )
        assert not model["converged"]
        assert model["n_saturated"] == 3

    @pytest.mark.parametrize("penalty", [0.0, 0.1])
    def test_lost_rank_marks_separated(self, penalty: float) -> None:
        """
        Test that a Hessian that loses rank after the first step marks separation.

        Along a separating direction the curvature vanishes as the coefficient grows,
        so the least-squares step stops moving it and looked converged. A penalty
        keeps curvature along every coefficient, so a penalized model keeps its rank.
        """
        model = {**self._model(), "penalty": penalty}
        _newton_update(model, self._sums(model, 0.4), 1e-8)
        assert not model["separated"]
        _newton_update(model, self._sums(model, 0.3, h=[300.0, 0.0, 0.0, 0.0]), 1e-8)
        assert model["separated"] == (penalty == 0)

    def test_step_halved_below_tol_has_converged(self) -> None:
        """Test that a model whose halved step is below ``tol`` has converged."""
        model = {**self._model(), "step": np.array([0.0, 1.5e-8])}
        _newton_update(model, self._sums(model, 0.6), 1e-8)
        assert model["converged"]


class TestTuneRegionalMutationRatePenalty:
    """Test the tune_regional_mutation_rate_penalty function."""

    @staticmethod
    def _expected_loss(
        rows: List[Dict], fold: str, penalty: float, features: List[str]
    ) -> float:
        """
        Return a fold's summed log loss under a model fit on the other rows.

        The fold is scored on the model's own features, without any it left out.
        """
        rows = [r for r in rows if all(r[f] is not None for f in ["y", *features])]
        fit_rows = [r for r in rows if r["fold"] != fold]
        model = _models(
            fit_regional_mutation_rate_models(
                _prepared(fit_rows, {rows[0]["context"]: features}), penalty
            )
        )[rows[0]["context"]]
        y, x = _arrays([r for r in rows if r["fold"] == fold], model.features)
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
        assert all(r.converged for r in tuning.values())
        for (penalty, fold), r in tuning.items():
            assert r.log_loss == pytest.approx(
                self._expected_loss(rows, fold, penalty, ["a", "b"]), rel=1e-6
            )

    def test_contexts_and_expressions(self) -> None:
        """
        Test that each context is tuned on its own rows, with other expressions.

        A held-out fold with only one label is still scored, and the folds that
        train on its rows are scored on models that include them.
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
        assert tuning[("ACG", "only_background")].n_test == 100
        for (context, fold), r in tuning.items():
            context_rows = aaa if context == "AAA" else acg + only_background
            assert r.log_loss == pytest.approx(
                self._expected_loss(context_rows, fold, 0.3, ["a", "b"]), rel=1e-6
            )

    def test_fits_that_differ_across_folds(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        """
        Test folds whose fits keep different features, or have no model.

        ACG's ``c`` varies only in odd rows, so the fit scored on odd rows, which is
        fit on even rows, leaves it out, and the one scored on even rows keeps it.
        CCG has no mutation in even rows, so the fit on them has no model and CCG's
        odd rows aren't scored, at any penalty. Penalties given as ints, including
        0, become floats.
        """
        rng = np.random.default_rng(11)
        rows = {
            "AAA": _training_rows(seed=3),
            "ACG": [
                {**r, "c": float(rng.normal()) if r["fold"] == "odd" else 0.0}
                for r in _training_rows(seed=4, context="ACG")
            ],
            "CCG": [
                {**r, "y": 0 if r["fold"] == "even" else r["y"]}
                for r in _training_rows(seed=5, context="CCG")
            ],
        }
        features = ["a", "b", "c"]
        with caplog.at_level(logging.WARNING):
            tuning = {
                (r.penalty, r.context, r.fold): r
                for r in tune_regional_mutation_rate_penalty(
                    _prepared(
                        [r for context_rows in rows.values() for r in context_rows],
                        features,
                        fold=True,
                    ),
                    [0, 1],
                ).collect()
            }
        assert sorted(tuning) == [
            (p, c, f)
            for p in [0.0, 1.0]
            for c, f in [
                ("AAA", "even"),
                ("AAA", "odd"),
                ("ACG", "even"),
                ("ACG", "odd"),
                ("CCG", "even"),
            ]
        ]
        for (penalty, context, fold), r in tuning.items():
            assert r.log_loss == pytest.approx(
                self._expected_loss(rows[context], fold, penalty, features), rel=1e-6
            )
        assert any(
            "no tuning model" in m and m.endswith(": CCG (odd).")
            for m in _warnings(caplog)
        )

    @pytest.mark.parametrize("fold", ["even", None])
    def test_one_fold(self, fold: Optional[str]) -> None:
        """Test that nothing is scored when every row is in one fold, or in none."""
        rows = [{**r, "fold": fold} for r in _training_rows()]
        ht = tune_regional_mutation_rate_penalty(
            _prepared(rows, {"AAA": ["a", "b"]}, fold=True), [0.01, 1000.0]
        )
        assert ht.count() == 0
        assert list(ht.row) == [
            "penalty",
            "context",
            "fold",
            "n_test",
            "log_loss",
            "converged",
        ]

    def test_non_converged_fits_warn(self, caplog: pytest.LogCaptureFixture) -> None:
        """
        Test that tuning fits that don't converge are counted, and still scored.

        Their rows say they didn't converge.
        """
        with caplog.at_level(logging.WARNING):
            ht = tune_regional_mutation_rate_penalty(
                _prepared(_training_rows(seed=2), {"AAA": ["a", "b"]}, fold=True),
                [0.3],
                max_iter=1,
            )
        assert [r.converged for r in ht.collect()] == [False, False]
        assert any(
            m.startswith("2 of 2 tuning models") and "penalties [0.3]" in m
            for m in _warnings(caplog)
        )

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
        (0.1, "AAA", "odd", 1.0),
        (0.1, "AAA", "even", 4.0),
        (0.1, "ACG", "odd", 3.0),
        (0.3, "AAA", "odd", 2.0),
        (0.3, "AAA", "even", 2.0),
        (0.3, "ACG", "odd", 3.5),
        (1.0, "AAA", "odd", 2.5),
        (1.0, "AAA", "even", 1.0),
        (1.0, "ACG", "odd", 4.5),
    ]
    """
    Summed losses of 8.0, 7.5, and 8.0 for penalties 0.1, 0.3, and 1.0.

    The lowest loss in any one context or fold is at another penalty.
    """

    @classmethod
    def _tuning_ht(
        cls,
        penalties: Optional[List[float]] = None,
        rows: Optional[List[tuple]] = None,
    ) -> hl.Table:
        """Return a tuning Table with the rows (default ``ROWS``) of the penalties."""
        return hl.Table.parallelize(
            [
                {"penalty": p, "context": c, "fold": f, "n_test": 10, "log_loss": l}
                for p, c, f, l in (cls.ROWS if rows is None else rows)
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
        """Test that the penalty with the lowest loss, over contexts and folds, wins."""
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

    def test_lower_end_of_0_does_not_warn(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Test that a chosen penalty of 0 doesn't warn: nothing lies below it."""
        rows = [(0.0, "AAA", "odd", 1.0), (1.0, "AAA", "odd", 2.0)]
        with caplog.at_level(logging.WARNING):
            assert (
                choose_regional_mutation_rate_penalty(self._tuning_ht(rows=rows), 1.0)
                == 0.0
            )
        assert not _warnings(caplog)

    @pytest.mark.parametrize("case", ["missing_row", "n_test"])
    def test_penalties_scored_on_different_rows_raise(self, case: str) -> None:
        """
        Test that penalties scored on different contexts, folds, or rows raise.

        A penalty missing a row, e.g. after a caller dropped fits that didn't
        converge, had a smaller summed loss and could win.
        """
        if case == "missing_row":
            ht = self._tuning_ht(rows=[r for r in self.ROWS if r[:2] != (0.3, "ACG")])
        else:
            ht = self._tuning_ht()
            ht = ht.annotate(
                n_test=hl.if_else(
                    (ht.penalty == 0.3) & (ht.context == "ACG"),
                    hl.int64(11),
                    ht.n_test,
                )
            )
        with pytest.raises(ValueError, match="weren't all scored on the same"):
            choose_regional_mutation_rate_penalty(ht, 1.0)

    def test_non_finite_loss_raises(self) -> None:
        """Test that a NaN summed loss raises; it used to win as the first penalty."""
        rows = [(p, c, f, float("nan") if p == 0.1 else l) for p, c, f, l in self.ROWS]
        with pytest.raises(ValueError, match="isn't finite at \\[0.1\\]"):
            choose_regional_mutation_rate_penalty(self._tuning_ht(rows=rows), 1.0)

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
        assert _rates(features_ht, _models_ht(), ["ACG", "AAA", "CCC"], 2.0)[1] == (
            pytest.approx({"AAA": _expit(0.5) / 0.6, "ACG": 2.0, "CCC": 1.0})
        )
        rates = _rates(features_ht, _models_ht(), ["ACG", "AAA", "CCC"], None)[1]
        assert rates["ACG"] is None
        assert rates["AAA"] == pytest.approx(_expit(0.5) / 0.6)
        assert rates["CCC"] == 1.0

    def test_missing_feature_rate_defaults_to_none(self) -> None:
        """
        Test that by default a context missing a model feature has a missing rate.

        A 1.0 filled in before normalizing ended up at 1 / factor; left missing, it
        stays out of the factors and can be filled with 1.0 afterwards.
        """
        features_ht = _features_ht([(1, 1001, 1.0, None)])
        ht = features_ht.annotate(
            adj_r=regional_mutation_rate_expr(
                features_ht.row, _models_ht(), ["AAA", "ACG"]
            )
        )
        rates = ht.collect()[0].adj_r
        assert rates["ACG"] is None
        assert rates["AAA"] == pytest.approx(_expit(0.5) / 0.6)

    def test_contexts(self) -> None:
        """Test that a context listed twice gets one rate, and no contexts no rates."""
        features_ht = _features_ht([(1, 1001, 0.5, 2.0)])
        rates = _rates(features_ht, _models_ht(), ["ACG", "ACG"])
        assert rates == {1: pytest.approx({"ACG": _expit(-1.0) / 0.5})}
        assert _rates(features_ht, _models_ht(), []) == {1: None}

    def test_no_models(self) -> None:
        """Test that with no models every window's rates are missing."""
        rates = _rates(
            _features_ht([(1, 1001, 0.0, 2.0)]), _models_ht().filter(False), ["AAA"]
        )
        assert rates == {1: None}

    @pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
    def test_non_finite_features_are_missing(self, value: float) -> None:
        """
        Test that a NaN or infinite feature value counts as a missing feature.

        A NaN rate passed the normalization's filter and made its context's factor,
        and every rate in that context, NaN; an infinite value gave a rate of 1/ave or
        0.
        """
        features_ht = _features_ht([(1, 1001, value, 2.0)])
        assert _rates(features_ht, _models_ht(), ["AAA", "ACG"])[1] == pytest.approx(
            {"AAA": 1.0, "ACG": _expit(-1.0) / 0.5}
        )
        assert _rates(features_ht, _models_ht(), ["AAA", "ACG"], None)[1]["AAA"] is None

    @pytest.mark.parametrize("skip", [False, True])
    def test_non_converged_models(
        self, caplog: pytest.LogCaptureFixture, skip: bool
    ) -> None:
        """
        Test that a model that didn't converge is applied, or skipped, with a warning.

        It was applied silently. With ``skip_non_converged`` its context gets 1.0, as
        a context without a model does.
        """
        models_ht = _models_ht()
        models_ht = models_ht.annotate(converged=models_ht.context != "AAA")
        features_ht = _features_ht([(1, 1001, 0.5, 2.0)])
        with caplog.at_level(logging.WARNING):
            ht = features_ht.annotate(
                adj_r=regional_mutation_rate_expr(
                    features_ht.row, models_ht, ["AAA", "ACG"], skip_non_converged=skip
                )
            )
            rates = ht.collect()[0].adj_r
        assert rates == pytest.approx(
            {"AAA": 1.0 if skip else _expit(0.0) / 0.6, "ACG": _expit(-1.0) / 0.5}
        )
        assert any(
            m.startswith("The models of AAA didn't converge") for m in _warnings(caplog)
        )

    def test_matches_the_fitted_model(self) -> None:
        """
        Test that the rate is the fitted model's prediction divided by ``ave``.

        The two features have different means and sds, so the test covers the mean
        subtraction, each feature's sd and coefficient, and the feature order, which
        isn't the order of the feature Table's fields. A window at the training means
        has rate 1, and a window missing a feature has no rate.
        """
        models_ht = fit_regional_mutation_rate_models(
            _prepared(_training_rows(), {"AAA": ["b", "a"]})
        )
        model = _models(models_ht)["AAA"]
        means = dict(zip(model.features, model.mean))
        values = [
            (-1.0, 9.0),
            (0.5, 5.0),
            (2.0, 1.0),
            (means["a"], means["b"]),
            (None, 1.0),
        ]
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
        # The values are (a, b); the model's features are (b, a).
        expected = _predict(model, np.array(values[:4])[:, ::-1]) / model.ave
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

    def test_missing_rates_are_left_out(self) -> None:
        """
        Test that a missing rate in a window is left out of its context's factor.

        ``regional_mutation_rate_expr`` leaves such rates missing when its
        ``missing_feature_rate`` is None.
        """
        rate_ht = _rate_ht()
        rate_ht = rate_ht.annotate(
            adj_r=hl.if_else(
                rate_ht.interval.start.position == 100,
                hl.dict({"AAA": hl.missing(hl.tfloat64), "ACG": 1.0, "CCC": 1.0}),
                rate_ht.adj_r,
            )
        )
        factors = _normalization_factors(_sites_ht(), rate_ht)
        assert factors == pytest.approx({"AAA": 2.0, "ACG": 0.5})

    @pytest.mark.parametrize("value", [float("nan"), float("inf")])
    def test_non_finite_rates_are_left_out(self, value: float) -> None:
        """
        Test that a NaN or infinite rate is left out of its context's factor.

        A NaN rate counts as defined in Hail, so it made the factor, and every
        normalized rate in that context, NaN.
        """
        rate_ht = _rate_ht()
        rate_ht = rate_ht.annotate(
            adj_r=hl.if_else(
                rate_ht.interval.start.position == 100,
                hl.dict({"AAA": value, "ACG": 1.0, "CCC": 1.0}),
                rate_ht.adj_r,
            )
        )
        factors = _normalization_factors(_sites_ht(), rate_ht)
        assert factors == pytest.approx({"AAA": 2.0, "ACG": 0.5})


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

    def test_missing_rates_stay_missing(self) -> None:
        """Test that a missing rate stays missing, and a missing dict stays missing."""
        ht = hl.Table.parallelize(
            [{"i": 0, "adj_r": {"AAA": None, "ACG": 1.0}}, {"i": 1, "adj_r": None}],
            hl.tstruct(i=hl.tint32, adj_r=hl.tdict(hl.tstr, hl.tfloat64)),
            key="i",
        )
        ht = ht.annotate(
            adj_r=normalize_regional_mutation_rate_expr(
                ht.adj_r, {"AAA": 2.0, "ACG": 0.5}
            )
        )
        assert [r.adj_r for r in ht.collect()] == [{"AAA": None, "ACG": 2.0}, None]
