"""
Functions to model a regional mutation rate from de novo mutations.

A regional mutation rate scales a site's expected variant count by how much more or
less often sites like it mutate in its region than on average. As in Gnocchi (Chen
et al., Nature 2024), one logistic model per trinucleotide context separates sites
with a de novo mutation (DNM) from a sample of non-mutated sites, using regional
genomic features of the window around each site. A window's rate is the model's
prediction for the window divided by its prediction for the average window, and it
is then normalized so that its expected-weighted mean over a set of sites is 1 in
each context.

The models are fit by Newton's method in Hail: each step aggregates every context's
gradient, Hessian, and loss over the training Table, and only the small per-context
linear systems are solved on the driver. The objective is the mean log loss, plus an
optional ridge penalty on the coefficients, and a step that raises it is halved. The
penalty can be chosen by the held-out log loss of models fit on all but one fold of
the training rows (:func:`tune_regional_mutation_rate_penalty` and
:func:`choose_regional_mutation_rate_penalty`). Alternatively, as in Gnocchi, each
context's features can first be screened one at a time
(:func:`screen_regional_mutation_rate_features`), and only those associated with the
label kept.
"""

import logging
import math
import sys
from typing import Any, Dict, List, Optional, Tuple, Union

import hail as hl
import numpy as np

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_utils")
logger.setLevel(logging.INFO)

MODELS_HT_SCHEMA = hl.tstruct(
    context=hl.tstr,
    features=hl.tarray(hl.tstr),
    mean=hl.tarray(hl.tfloat64),
    sd=hl.tarray(hl.tfloat64),
    coefs=hl.tarray(hl.tfloat64),
    intercept=hl.tfloat64,
    ave=hl.tfloat64,
    n_mutated=hl.tint32,
    n_background=hl.tint32,
    converged=hl.tbool,
)
"""
Row type of a regional mutation rate models Table.

:func:`fit_regional_mutation_rate_models` returns a Table of this type, and
:func:`regional_mutation_rate_expr` applies one.
"""

_SATURATED_PROBABILITY = 10 * sys.float_info.epsilon
"""
A fitted probability within this of 0 or 1 is numerically 0 or 1, as in R's glm.fit.
"""

_MODEL_EXPR_TYPE = hl.tstruct(
    positions=hl.tarray(hl.tint32),
    mean=hl.tarray(hl.tfloat64),
    sd=hl.tarray(hl.tfloat64),
    coefs=hl.tarray(hl.tfloat64),
    intercept=hl.tfloat64,
)


def _models_by_group_expr(
    models: Dict[str, Dict[str, Any]],
) -> hl.expr.DictExpression:
    """
    Get a dict expression of group to model, for evaluating models on rows.

    :param models: Dictionary of group to model, each with ``positions`` (the
        positions of its features in a row's feature array), ``mean``, ``sd``,
        ``coefs``, and ``intercept``.
    :return: Dict expression of group to model struct.
    """
    return hl.literal(
        {
            g: hl.Struct(**{f: m[f] for f in _MODEL_EXPR_TYPE})
            for g, m in models.items()
        },
        hl.tdict(hl.tstr, _MODEL_EXPR_TYPE),
    )


def _standardized_expr(
    x_expr: hl.expr.ArrayNumericExpression, model_expr: hl.expr.StructExpression
) -> hl.expr.ArrayNumericExpression:
    """
    Get the standardized values ``(x - mean) / sd`` of a model's features.

    :param x_expr: Feature values, indexed by the model's ``positions``.
    :param model_expr: Model struct (see :func:`_models_by_group_expr`).
    :return: Standardized value of each of the model's features, in its order.
    """
    return hl.range(hl.len(model_expr.positions)).map(
        lambda i: (x_expr[model_expr.positions[i]] - model_expr.mean[i])
        / model_expr.sd[i]
    )


def _linear_predictor_expr(
    z_expr: hl.expr.ArrayNumericExpression, model_expr: hl.expr.StructExpression
) -> hl.expr.Float64Expression:
    """
    Get the linear predictor ``intercept + sum(coefs * z)`` of a model.

    It is missing if any of the standardized values is missing.

    :param z_expr: Standardized feature values (:func:`_standardized_expr`).
    :param model_expr: Model struct (see :func:`_models_by_group_expr`).
    :return: Linear predictor expression.
    """
    return model_expr.intercept + hl.sum(
        hl.zip(model_expr.coefs, z_expr).map(lambda t: t[0] * t[1]),
        filter_missing=False,
    )


def _log_loss_expr(
    label_expr: hl.expr.Float64Expression, eta_expr: hl.expr.Float64Expression
) -> hl.expr.Float64Expression:
    """
    Get a row's log loss ``log(1 + exp(eta)) - label * eta``.

    It is written as ``log(1 + exp(-|eta|)) + max(eta, 0) - label * eta`` to avoid
    overflow.

    :param label_expr: 0/1 label.
    :param eta_expr: Linear predictor (:func:`_linear_predictor_expr`).
    :return: Log loss expression.
    """
    return (
        hl.log(1 + hl.exp(-hl.abs(eta_expr)))
        + hl.max(eta_expr, 0)
        - label_expr * eta_expr
    )


def prepare_regional_mutation_rate_training_ht(
    ht: hl.Table,
    features: Union[List[str], Dict[str, List[str]]],
    label_expr: Optional[hl.expr.NumericExpression] = None,
    context_expr: Optional[hl.expr.StringExpression] = None,
    fold_expr: Optional[hl.expr.StringExpression] = None,
) -> hl.Table:
    """
    Prepare a training Table for fitting and tuning regional mutation rate models.

    Each row gets its 0/1 label, its context, its fold if ``fold_expr`` is given, and
    the values of its features in order: the same ``features`` for every context, or
    its context's features when ``features`` is a dict. ``features`` is recorded in
    the globals, so the fit and the tuning use the same ones. Rows with a missing
    value are kept; :func:`fit_regional_mutation_rate_models` and
    :func:`tune_regional_mutation_rate_penalty` leave them out.

    :param ht: Training Table with the feature fields.
    :param features: Feature fields to model, which must be numeric or boolean and
        are converted to float64: one list for every context, or a dict of context to
        its own list.
    :param label_expr: 0/1 label of each row (1 for a site with a DNM), a number or a
        boolean. Default is ``ht.y``.
    :param context_expr: Context of each row. Default is ``ht.context``.
    :param fold_expr: Fold of each row, a string such as ``"odd"`` or ``"even"``
        for the chromosome number, needed for
        :func:`tune_regional_mutation_rate_penalty`. Rows with a missing fold are
        only used for fitting. Default is None.
    :return: Table keyed as ``ht`` with only ``label`` (float64), ``context``,
        ``fold`` (if ``fold_expr`` is given), and ``x``, the float64 values of the
        row's features in order (missing, with a dict, for a context not in it). Its
        only global is ``features``.
    """
    # hl.float64 also parses strings, so a string field would only fail, or be
    # misread, when the models are fit.
    feature_fields = (
        [f for fs in features.values() for f in fs]
        if isinstance(features, dict)
        else features
    )
    not_numeric = [
        f
        for f in dict.fromkeys(feature_fields)
        if not isinstance(ht[f], hl.expr.NumericExpression)
    ]
    if not_numeric:
        raise ValueError(
            "Features must be numeric or boolean fields; these are not:"
            f" {', '.join(not_numeric)}."
        )
    if fold_expr is not None and fold_expr.dtype != hl.tstr:
        raise ValueError(
            f"fold_expr must be a string expression, not {fold_expr.dtype}."
        )

    context_expr = ht.context if context_expr is None else context_expr
    # x holds a row's features in one array, so they are all converted to float64.
    if isinstance(features, dict):
        # Each row gets its own context's features.
        x = hl.switch(context_expr)
        for context, context_features in features.items():
            x = x.when(context, hl.array([hl.float64(ht[f]) for f in context_features]))
        x = x.or_missing()
        features_expr = hl.literal(features, hl.tdict(hl.tstr, hl.tarray(hl.tstr)))
    else:
        x = hl.array([hl.float64(ht[f]) for f in features])
        features_expr = hl.literal(list(features), hl.tarray(hl.tstr))
    fields = {
        "label": hl.float64(ht.y if label_expr is None else label_expr),
        "context": context_expr,
        **({"fold": fold_expr} if fold_expr is not None else {}),
    }
    ht = ht.select(**fields, x=x)

    return ht.select_globals(features=features_expr)


def _complete_training_rows(
    ht: hl.Table, fields: List[str]
) -> Tuple[hl.Table, Dict[str, List[str]]]:
    """
    Check a prepared training Table and drop the rows with a missing label or feature.

    NaN counts as defined in Hail, and one NaN or infinite feature value would make
    its context's model NaN, so the remaining rows' features must all be finite, and
    their labels 0 or 1.

    :param ht: Table from :func:`prepare_regional_mutation_rate_training_ht`.
    :param fields: Row fields ``ht`` must have.
    :return: ``ht`` without incomplete rows, and the features of each context, from
        its ``features`` global. With one list for every context, the contexts are
        those of the remaining rows.
    """
    missing = [f for f in fields if f not in ht.row]
    if "features" not in ht.globals:
        missing.append("the features global")
    if missing:
        raise ValueError(
            f"The training Table has no {', '.join(missing)}. Prepare it with"
            " prepare_regional_mutation_rate_training_ht"
            + (", with a fold_expr." if "fold" in missing else ".")
        )
    ht = ht.filter(
        hl.is_defined(ht.label) & hl.is_defined(ht.x) & hl.all(hl.is_defined, ht.x)
    )

    rows = ht.aggregate(
        hl.struct(
            n_not_finite=hl.agg.count_where(~hl.all(hl.is_finite, ht.x)),
            # NaN, like any other label but 0 or 1, is counted here.
            n_bad_label=hl.agg.count_where((ht.label != 0) & (ht.label != 1)),
            contexts=hl.agg.collect_as_set(ht.context),
        )
    )
    if rows.n_not_finite:
        raise ValueError(
            f"{rows.n_not_finite} training rows have a NaN or infinite feature value."
            " Make those values missing, or remove the rows."
        )
    if rows.n_bad_label:
        raise ValueError(
            f"{rows.n_bad_label} training rows have a label other than 0 or 1."
        )

    features = hl.eval(ht.features)
    if isinstance(features, list):
        # With one list for every context, model each context the rows have.
        features = {c: features for c in sorted(rows.contexts - {None})}

    return ht, {c: list(f) for c, f in features.items()}


def _init_logistic_models(
    ht: hl.Table,
    group_features: Dict[str, List[str]],
    ridge_penalties: Dict[str, float],
    log: bool,
) -> Dict[str, Dict[str, Any]]:
    """
    Set up a model for each group with rows of both labels and a varying feature.

    One pass over ``ht`` gets each group's counts and its features' means, minimums,
    and maximums. Features with a single value carry no information and can't be
    standardized, so they are left out. A second pass gets the standard deviations
    (:func:`_feature_sds`). Each model starts at the intercept of the group's label
    rate, with every coefficient 0.

    :param ht: Complete training rows (:func:`_complete_training_rows`) with a
        ``group`` field.
    :param group_features: Features of each group, in the order of its rows' ``x``.
    :param ridge_penalties: Per-row ridge penalty of each group.
    :param log: Whether to log groups without a model and features left out.
    :return: Dictionary of group to model, for the groups that have one.
    """
    stats = ht.aggregate(
        hl.agg.group_by(
            ht.group,
            hl.struct(
                n=hl.agg.count(),
                n_mutated=hl.agg.sum(ht.label),
                x=hl.agg.array_agg(
                    lambda v: hl.struct(
                        mean=hl.agg.mean(v), min=hl.agg.min(v), max=hl.agg.max(v)
                    ),
                    ht.x,
                ),
            ),
        )
    )
    models = {}
    for group, features in group_features.items():
        s = stats.get(group)
        # Compare the min and max rather than test sd > 0: a constant feature's
        # computed mean can differ from its value by rounding, leaving a tiny sd.
        keep = [i for i, v in enumerate(s.x) if v.min < v.max] if s else []
        # The label is float64, so its sum is too.
        n_mutated = int(round(s.n_mutated)) if s else 0
        if not keep or not 0 < n_mutated < s.n:
            if log:
                logger.info(
                    "No model for %s: no rows of one label, or no feature varies.",
                    group,
                )
            continue
        if log and len(keep) < len(features):
            constant = [f for i, f in enumerate(features) if i not in keep]
            logger.info("%s: features left out (constant): %s", group, constant)
        models[group] = {
            "features": [features[i] for i in keep],
            "positions": keep,
            "mean": [s.x[i].mean for i in keep],
            "sd": [1.0] * len(keep),
            "n": s.n,
            "n_mutated": n_mutated,
            "intercept": float(np.log(n_mutated / (s.n - n_mutated))),
            "coefs": [0.0] * len(keep),
            "converged": False,
            "penalty": ridge_penalties[group],
            # The objective at the last accepted coefficients (none yet).
            "objective": np.inf,
        }
    if models:
        for group, sds in _feature_sds(ht, models).items():
            models[group]["sd"] = sds

    return models


def _feature_sds(
    ht: hl.Table, models: Dict[str, Dict[str, Any]]
) -> Dict[str, List[float]]:
    """
    Get the standard deviation (ddof=0) of each model's features, around their means.

    ``hl.agg.stats`` gets the sd from the sum of squares in the same pass as the
    mean, which loses precision when a feature varies little next to its mean: a
    nearly constant feature can get an sd of NaN. Summing the squared deviations from
    the means, in a second pass, doesn't.

    :param ht: Complete training rows (:func:`_complete_training_rows`) with a
        ``group`` field.
    :param models: Models by group, with every ``sd`` still 1.
    :return: Dictionary of group to the sds of its model's features.
    """
    model_expr = _models_by_group_expr(models).get(ht.group)
    # With every sd 1, the standardized values are the deviations from the means.
    deviations = _standardized_expr(ht.x, model_expr)
    sum_sq = ht.aggregate(
        hl.agg.filter(
            hl.is_defined(model_expr),
            hl.agg.group_by(ht.group, hl.agg.array_sum(deviations.map(lambda d: d**2))),
        )
    )

    return {g: [float(np.sqrt(s / models[g]["n"])) for s in sum_sq[g]] for g in models}


def _newton_sums(
    ht: hl.Table, models: Dict[str, Dict[str, Any]]
) -> Dict[str, hl.Struct]:
    """
    Sum each group's log loss, gradient, and Hessian terms at its current model.

    With ``v = [1, z]`` for a row's standardized features and ``p`` its predicted
    probability, ``g`` sums ``(label - p) * v``, ``h`` sums ``p * (1 - p) * v * v^T``
    (flattened), and ``loss`` sums the log loss, in one pass over ``ht``.
    ``n_saturated`` counts the rows whose ``p`` is numerically 0 or 1.

    :param ht: Complete training rows (:func:`_complete_training_rows`) with a
        ``group`` field.
    :param models: Models to evaluate, by group.
    :return: Dictionary of group to a struct with ``g``, ``h``, ``loss``, and
        ``n_saturated``.
    """
    model_expr = _models_by_group_expr(models).get(ht.group)
    z = _standardized_expr(ht.x, model_expr)
    eta = _linear_predictor_expr(z, model_expr)
    p = hl.expit(eta)
    # The leading 1 is the intercept's term.
    v = hl.array([1.0]).extend(z)

    # Skip the rows of groups without a model, or whose model has converged.
    return ht.aggregate(
        hl.agg.filter(
            hl.is_defined(model_expr),
            hl.agg.group_by(
                ht.group,
                hl.struct(
                    g=hl.agg.array_sum(v.map(lambda a: a * (ht.label - p))),
                    h=hl.agg.array_sum(
                        hl.flatten(
                            v.map(lambda a: v.map(lambda b: p * (1 - p) * a * b))
                        )
                    ),
                    loss=hl.agg.sum(_log_loss_expr(ht.label, eta)),
                    n_saturated=hl.agg.count_where(
                        (p < _SATURATED_PROBABILITY) | (p > 1 - _SATURATED_PROBABILITY)
                    ),
                ),
            ),
        )
    )


def _newton_update(model: Dict[str, Any], sums: hl.Struct, tol: float) -> None:
    """
    Take one step of Newton's method, with step halving, on a model's objective.

    The objective is mean log loss + ``penalty`` / 2 * sum(coef ** 2), with the
    intercept unpenalized. If the last step raised it by more than the rounding in
    the summed loss, the step overshot, so the model goes back and takes half of it.
    Otherwise the next Newton step is taken from the current coefficients. The step
    is solved by least squares, so collinear features, which leave an unpenalized
    Hessian singular, still give the smallest step that solves it. The model has
    converged when a step, full or halved, changes no coefficient by more than
    ``tol``; a halved step that small means no step lowers the objective by more
    than its rounding.

    The Hessian's rank at the first step is kept. If it later falls, a direction has
    lost its curvature: along it the fitted probabilities of some rows are being
    pushed toward 0 or 1, as when the features separate the labels. The least-squares
    step then stops moving in that direction, so the model is marked ``separated``.
    A penalty adds curvature along every coefficient, so it keeps the rank.

    :param model: Model to update in place.
    :param sums: The model's loss, gradient, and Hessian sums (:func:`_newton_sums`).
    :param tol: Convergence tolerance.
    :return: None.
    """
    k = len(model["coefs"]) + 1
    beta = np.array([model["intercept"], *model["coefs"]])
    # Penalize every coefficient but the intercept.
    penalty = model["penalty"] * np.r_[0.0, np.ones(k - 1)]
    objective = sums.loss / model["n"] + np.sum(penalty * beta**2) / 2
    # Near the minimum, rounding in the mean loss can be larger than a step's real
    # change, so a smaller rise isn't an overshoot. Summing n rows rounds by up to
    # about n * eps of the sum, and each row's loss by about eps on its own.
    rounding = sys.float_info.epsilon * (model["n"] * model["objective"] + 1)
    model["halving"] = objective > model["objective"] + rounding
    if model["halving"]:
        model["step"] = model["step"] / 2
        model["converged"] = bool(np.max(np.abs(model["step"])) < tol)
    else:
        grad = -np.array(sums.g) / model["n"] + penalty * beta
        hess = np.array(sums.h).reshape(k, k) / model["n"] + np.diag(penalty)
        model["objective"] = objective
        model["n_saturated"] = sums.n_saturated
        model["start"] = beta
        model["step"], _, rank, _ = np.linalg.lstsq(hess, grad, rcond=None)
        model.setdefault("first_rank", rank)
        model["separated"] = model.get("separated", False) or rank < model["first_rank"]
        model["converged"] = bool(np.max(np.abs(model["step"])) < tol)
    beta = model["start"] - model["step"]
    model["intercept"] = float(beta[0])
    model["coefs"] = beta[1:].tolist()


def _fit_logistic_models(
    ht: hl.Table,
    group_features: Dict[str, List[str]],
    ridge_penalties: Dict[str, float],
    max_iter: int,
    tol: float,
    log: bool = True,
) -> Dict[str, Dict[str, Any]]:
    """
    Fit a logistic model, with a ridge penalty, for each group of a prepared Table.

    A group is a context, or, when tuning, a context in one penalty and fold. All
    groups are fit together by Newton's method (:func:`_newton_update`), with one
    pass over ``ht`` per step. A model that converged with a halved step, or ran out
    of steps, keeps its last accepted coefficients, the last point evaluated.

    A model marked ``separated`` (its Hessian lost rank) is not converged: its
    features separate the labels, at least in part, so the log loss has no minimum
    and the coefficient along that direction is wherever the steps stopped. For an
    unpenalized model, fitted probabilities numerically 0 or 1 are only logged, as
    R's glm.fit warns: they can come from separation, but also from a heavy-tailed
    feature at a real minimum. A penalty gives the objective a minimum, so penalized
    models aren't checked for them.

    :param ht: Complete training rows (:func:`_complete_training_rows`) with a
        ``group`` field.
    :param group_features: Features of each group, in the order of its rows' ``x``.
    :param ridge_penalties: Per-row ridge penalty of each group, 0 for none.
    :param max_iter: Maximum number of passes over ``ht`` for Newton steps,
        including halved ones.
    :param tol: A fit has converged when a step, full or halved, changes no
        coefficient by more than this.
    :param log: Whether to log groups without a model, features left out, and fits
        that don't converge. Default is True.
    :return: Dictionary of group to model, for the groups that have one.
    """
    models = _init_logistic_models(ht, group_features, ridge_penalties, log)
    if not models:
        return models

    for _ in range(max_iter):
        # Converged models keep their coefficients and leave the aggregation.
        todo = {g: m for g, m in models.items() if not m["converged"]}
        if not todo:
            break
        sums = _newton_sums(ht, todo)
        for group, model in todo.items():
            _newton_update(model, sums[group], tol)

    for model in models.values():
        # A halved step, or a step taken just before max_iter ran out, leads to a
        # point never evaluated, so the model goes back to its last accepted
        # coefficients, the best it knows.
        if "start" in model and (model["halving"] or not model["converged"]):
            model["intercept"] = float(model["start"][0])
            model["coefs"] = model["start"][1:].tolist()
        if model.get("separated"):
            model["converged"] = False

    if log:
        for group, model in models.items():
            if model.get("separated"):
                logger.warning(
                    "The model fit for %s did not converge: its log loss lost its"
                    " curvature in one direction, so its features separate the"
                    " labels, at least in part.",
                    group,
                )
            elif not model["converged"]:
                logger.warning("The model fit for %s did not converge.", group)
            if model.get("n_saturated") and not model["penalty"]:
                logger.warning(
                    "%s: %d training rows have a fitted probability numerically 0 or"
                    " 1, so its features may separate the labels.",
                    group,
                    model["n_saturated"],
                )
        logger.info("%d of %d contexts modeled.", len(models), len(group_features))

    return models


def _check_penalties(penalties: List[float]) -> None:
    """
    Check that ridge penalties were given and that each is finite and at least 0.

    :param penalties: Per-row ridge penalties.
    :return: None.
    """
    if not penalties:
        raise ValueError("No ridge penalties were given.")
    bad = [p for p in penalties if not (np.isfinite(p) and p >= 0)]
    if bad:
        raise ValueError(f"Ridge penalties must be finite and at least 0, not {bad}.")


def fit_regional_mutation_rate_models(
    ht: hl.Table,
    ridge_penalty: float = 0.0,
    max_iter: int = 50,
    tol: float = 1e-8,
) -> hl.Table:
    """
    Fit a logistic regional mutation rate model for each context.

    For each context, rows missing the label or any of its features are dropped, and
    features with a single value in the remaining rows are left out, since they carry
    no information and can't be standardized. The features are standardized with
    their mean and standard deviation (ddof=0), and the fit minimizes mean log loss +
    ``ridge_penalty`` / 2 * sum(coef ** 2), with the intercept unpenalized. A site's
    linear predictor is ``intercept + sum(coefs * (x - mean) / sd)``.

    The fit is by Newton's method, and each step is a pass over ``ht``, so ``ht``
    should be read from disk or checkpointed. A step that raises the objective is
    halved and tried again. A fit whose last step was halved, because the halved
    step fell below ``tol`` or ``max_iter`` ran out, keeps its last accepted
    coefficients. Each step is solved by least squares, so collinear features still
    fit without a penalty: their coefficients aren't unique, but the model's
    predictions are. Without a penalty, features that separate the labels, even in
    part, leave the loss without a minimum: the curvature along the separating
    direction vanishes as its coefficient grows, and when the Hessian loses rank the
    fit is marked not converged. Fitted probabilities numerically 0 or 1 in an
    unpenalized fit are logged, as R's glm.fit warns, without changing
    ``converged``. Contexts with no model, features left out because they are
    constant, and fits that don't converge are logged. Feature values must be finite
    and labels 0 or 1: anything else raises an error.

    :param ht: Training Table from
        :func:`prepare_regional_mutation_rate_training_ht`.
    :param ridge_penalty: Per-row ridge penalty, finite and at least 0. Default is
        0.0, no penalty.
    :param max_iter: Maximum number of Newton steps, including halved ones. Default
        is 50.
    :param tol: The fit has converged when a step, full or halved, changes no
        coefficient by more than this. Default is 1e-8.
    :return: Table keyed by ``context`` with one row per modeled context: the kept
        ``features``, their training ``mean`` and ``sd``, ``coefs``, ``intercept``,
        ``ave`` (the predicted probability with every feature at its mean),
        ``n_mutated`` and ``n_background`` (rows of each label), and ``converged``
        (see :data:`MODELS_HT_SCHEMA`), and a ``ridge_penalty`` global (float64).
    """
    _check_penalties([ridge_penalty])
    ht, context_features = _complete_training_rows(ht, ["label", "context", "x"])
    models = _fit_logistic_models(
        ht.annotate(group=ht.context),
        context_features,
        dict.fromkeys(context_features, ridge_penalty),
        max_iter,
        tol,
    )

    return hl.Table.parallelize(
        [
            {
                "context": context,
                "features": m["features"],
                "mean": m["mean"],
                "sd": m["sd"],
                "coefs": m["coefs"],
                "intercept": m["intercept"],
                # At the feature means every standardized value is 0.
                "ave": float(1.0 / (1.0 + np.exp(-m["intercept"]))),
                "n_mutated": m["n_mutated"],
                "n_background": m["n"] - m["n_mutated"],
                "converged": m["converged"],
            } for context, m in models.items()
        ],
        MODELS_HT_SCHEMA,
        key="context",
    ).annotate_globals(ridge_penalty=float(ridge_penalty))


def screen_regional_mutation_rate_features(
    ht: hl.Table,
    max_iter: int = 50,
    tol: float = 1e-8,
) -> hl.Table:
    """
    Test each feature of each context on its own for an association with the label.

    This is the univariate screen of Gnocchi's feature selection. For each context and
    each of its features, an unpenalized logistic model of the label on that feature
    alone, standardized (ddof=0), is fit on the context's rows where the label and that
    feature are defined, so a row missing another feature still counts. The Wald test
    of the feature's coefficient gives its standard error, from the inverse of the
    Hessian of the summed log loss at the fit, and its two-sided p-value. All the fits
    are made together, with one pass over the rows per Newton step, so ``ht`` should
    be read from disk or checkpointed.

    A fit that doesn't converge, for instance because the feature separates the
    labels, has a missing ``se`` and ``p_value``, so a cutoff on the p-value never
    selects it; these fits are counted in a warning. A feature with a single value in a
    context's rows, and a context whose rows all have one label, have no row. Feature
    values must be finite and labels 0 or 1: anything else raises an error.

    :param ht: Training Table from
        :func:`prepare_regional_mutation_rate_training_ht`.
    :param max_iter: Maximum number of Newton steps per fit, including halved ones.
        Default is 50.
    :param tol: Convergence tolerance of each fit. Default is 1e-8.
    :return: Table keyed by ``context`` and ``feature`` (str) with the feature's
        coefficient per standard deviation (``coef``), its standard error (``se``) and
        two-sided Wald p-value (``p_value``), all float64, the number of rows fit
        (``n``) and of mutated rows among them (``n_mutated``), both int32, and whether
        the fit converged (``converged``, bool).
    """
    missing = [f for f in ["label", "context", "x"] if f not in ht.row]
    if "features" not in ht.globals:
        missing.append("the features global")
    if missing:
        raise ValueError(
            f"The training Table has no {', '.join(missing)}. Prepare it with"
            " prepare_regional_mutation_rate_training_ht."
        )
    features = hl.eval(ht.features)
    ht = ht.filter(hl.is_defined(ht.label) & hl.is_defined(ht.context))
    if isinstance(features, list):
        # With one list for every context, screen each context the rows have.
        contexts = ht.aggregate(hl.agg.collect_as_set(ht.context))
        features = {c: features for c in sorted(contexts)}
    features = {c: list(f) for c, f in features.items()}

    # One row per row and feature, keeping the rows where that feature is defined.
    names = hl.literal(features, hl.tdict(hl.tstr, hl.tarray(hl.tstr))).get(ht.context)
    ht = ht.annotate(
        _names=names,
        _i=hl.or_missing(
            hl.is_defined(names) & hl.is_defined(ht.x), hl.range(hl.len(ht.x))
        ),
    )
    ht = ht.explode("_i")
    ht = ht.filter(hl.is_defined(ht.x[ht._i]))
    ht = ht.select(
        "label",
        group=ht.context + "|" + ht._names[ht._i],
        x=hl.array([ht.x[ht._i]]),
    )
    rows = ht.aggregate(
        hl.struct(
            n_not_finite=hl.agg.count_where(~hl.is_finite(ht.x[0])),
            n_bad_label=hl.agg.count_where((ht.label != 0) & (ht.label != 1)),
        )
    )
    if rows.n_not_finite:
        raise ValueError(
            f"{rows.n_not_finite} training rows have a NaN or infinite feature value."
            " Make those values missing, or remove the rows."
        )
    if rows.n_bad_label:
        raise ValueError(
            f"{rows.n_bad_label} training rows have a label other than 0 or 1."
        )

    groups = {f"{c}|{f}": (c, f) for c, fs in features.items() for f in fs}
    models = _fit_logistic_models(
        ht,
        {g: [f] for g, (_, f) in groups.items()},
        dict.fromkeys(groups, 0.0),
        max_iter,
        tol,
        log=False,
    )
    not_converged = sorted(g for g, m in models.items() if not m["converged"])
    if not_converged:
        logger.warning(
            "%d of %d screen fits (one per context and feature) did not converge, so"
            " they have no p-value: %s.",
            len(not_converged),
            len(models),
            ", ".join(not_converged),
        )

    # The Hessian of the summed log loss at the fit, inverted, is the covariance of
    # the intercept and coefficient.
    sums = _newton_sums(ht, models) if models else {}
    results = []
    for group, m in models.items():
        se = p_value = None
        if m["converged"]:
            cov = np.linalg.inv(np.array(sums[group].h).reshape(2, 2))
            se = float(np.sqrt(cov[1, 1]))
            p_value = math.erfc(abs(m["coefs"][0]) / se / math.sqrt(2))
        context, feature = groups[group]
        results.append(
            {
                "context": context,
                "feature": feature,
                "coef": m["coefs"][0],
                "se": se,
                "p_value": p_value,
                "n": m["n"],
                "n_mutated": m["n_mutated"],
                "converged": m["converged"],
            }
        )

    return hl.Table.parallelize(
        results,
        hl.tstruct(
            context=hl.tstr,
            feature=hl.tstr,
            coef=hl.tfloat64,
            se=hl.tfloat64,
            p_value=hl.tfloat64,
            n=hl.tint32,
            n_mutated=hl.tint32,
            converged=hl.tbool,
        ),
        key=["context", "feature"],
    )


def _rows_by_fit(ht: hl.Table, fits: List[hl.Struct], held_out: bool) -> hl.Table:
    """
    Copy each training row into every fit it is a training, or held-out, row of.

    Fit ``i`` trains on the rows outside ``fits[i].fold`` (including rows with a
    missing fold) and is scored on the rows in it.

    :param ht: Complete training rows (:func:`_complete_training_rows`) with a
        ``fold``.
    :param fits: Penalty and fold of each fit.
    :param held_out: Whether to copy each row into the fits it is held out of,
        rather than those it trains.
    :return: ``ht`` with one row per row and fit, and a ``group`` field naming the
        fit and the row's context, ``"{i}|{context}"``.
    """
    fit_folds = hl.literal([fit.fold for fit in fits], hl.tarray(hl.tstr))
    # A row with a missing fold trains every fit and is never held out.
    fit_ids = hl.range(len(fits)).filter(
        lambda i: (
            (fit_folds[i] == ht.fold)
            if held_out
            else hl.or_else(fit_folds[i] != ht.fold, True)
        )
    )
    ht = ht.annotate(fit=fit_ids).explode("fit")

    return ht.annotate(group=hl.str(ht.fit) + "|" + ht.context)


def _held_out_log_loss(
    ht: hl.Table, models: Dict[str, Dict[str, Any]]
) -> Dict[str, hl.Struct]:
    """
    Sum each group's log loss under its model, in one pass over ``ht``.

    :param ht: Held-out rows with a ``group`` field (:func:`_rows_by_fit`).
    :param models: Models by group; groups without one are skipped.
    :return: Dictionary of group to a struct with the number of rows (``n_test``)
        and their summed log loss (``log_loss``).
    """
    model_expr = _models_by_group_expr(models).get(ht.group)
    eta = _linear_predictor_expr(_standardized_expr(ht.x, model_expr), model_expr)

    return ht.aggregate(
        hl.agg.filter(
            hl.is_defined(model_expr),
            hl.agg.group_by(
                ht.group,
                hl.struct(
                    n_test=hl.agg.count(),
                    log_loss=hl.agg.sum(_log_loss_expr(ht.label, eta)),
                ),
            ),
        )
    )


def tune_regional_mutation_rate_penalty(
    ht: hl.Table,
    penalties: List[float],
    max_iter: int = 50,
    tol: float = 1e-8,
) -> hl.Table:
    """
    Get the held-out log loss of regional mutation rate models for each ridge penalty.

    For each penalty, context, and fold, the context's model is fit on the rows
    outside the fold (:func:`fit_regional_mutation_rate_models`) and scored on the
    rows in it. Rows with a missing fold are only used for fitting. All the fits are
    made together, with one pass over every fit's rows per Newton step, so ``ht``
    should be read from disk or checkpointed. Fits that don't converge are counted,
    with their penalties, in a warning and still scored, and contexts and folds with
    no model are logged. Feature values must be finite: a NaN or infinite value
    raises an error.

    :param ht: Training Table from
        :func:`prepare_regional_mutation_rate_training_ht`, with a ``fold``.
    :param penalties: Per-row ridge penalties to try, each finite and at least 0.
        Repeats are tried once.
    :param max_iter: Maximum number of Newton steps per fit, including halved ones.
        Default is 50.
    :param tol: Convergence tolerance of each fit. Default is 1e-8.
    :return: Table keyed by ``penalty`` (float64), ``context`` (str), and ``fold``
        (str), with the number of held-out rows with no missing label or feature
        (``n_test``, int64), their summed log loss (``log_loss``, float64), and
        whether the fit converged (``converged``, bool). A
        context and fold with no model (its training rows have only one label, or no
        feature varies) or no held-out rows has no rows, at every penalty alike, so
        each penalty's loss is summed over the same contexts and folds. The Table is
        empty when no held-out rows could be scored.
    """
    _check_penalties(penalties)
    ht, context_features = _complete_training_rows(
        ht, ["label", "context", "fold", "x"]
    )
    folds = sorted(ht.aggregate(hl.agg.collect_as_set(ht.fold)) - {None})
    # One fit per penalty and fold, on the rows outside the fold. They are all fit,
    # and scored, together: each row is copied into every fit it belongs to.
    fits = [
        hl.Struct(penalty=p, fold=f)
        for p in sorted({float(p) for p in penalties})
        for f in folds
    ]
    # Each context in each fit is its own model, so its own group.
    groups = {
        f"{i}|{c}": hl.Struct(**fit, context=c, features=features)
        for i, fit in enumerate(fits)
        for c, features in context_features.items()
    }
    models = _fit_logistic_models(
        _rows_by_fit(ht, fits, held_out=False),
        {g: group.features for g, group in groups.items()},
        {g: group.penalty for g, group in groups.items()},
        max_iter,
        tol,
        log=False,
    )
    # Whether a context has a model in a fold depends only on its training rows, so
    # it is the same at every penalty.
    no_model = sorted(
        {(group.context, group.fold) for group in groups.values()}
        - {(groups[g].context, groups[g].fold) for g in models}
    )
    if no_model:
        logger.warning(
            "These contexts and folds have no tuning model, since their training rows"
            " have one label or no feature varies, and are left out: %s.",
            ", ".join(f"{c} ({f})" for c, f in no_model),
        )
    not_converged = [g for g, m in models.items() if not m["converged"]]
    if not_converged:
        logger.warning(
            "%d of %d tuning models (one per penalty, fold, and context) did not"
            " converge, at penalties %s; they are still scored.",
            len(not_converged),
            len(models),
            sorted({groups[g].penalty for g in not_converged}),
        )
    scores = _held_out_log_loss(_rows_by_fit(ht, fits, held_out=True), models)

    return hl.Table.parallelize(
        [
            {
                "penalty": groups[g].penalty,
                "context": groups[g].context,
                "fold": groups[g].fold,
                **score,
                "converged": models[g]["converged"],
            }
            for g, score in scores.items()
        ],
        hl.tstruct(
            penalty=hl.tfloat64,
            context=hl.tstr,
            fold=hl.tstr,
            n_test=hl.tint64,
            log_loss=hl.tfloat64,
            converged=hl.tbool,
        ),
        key=["penalty", "context", "fold"],
    )


def choose_regional_mutation_rate_penalty(
    ht: hl.Table, default_penalty: float
) -> float:
    """
    Choose the ridge penalty with the lowest held-out log loss.

    The log loss is summed over contexts and folds, so one penalty is chosen for
    every context; on a tie, the smaller penalty is chosen. Every penalty must have
    been scored on the same contexts, folds, and numbers of held-out rows, as
    :func:`tune_regional_mutation_rate_penalty` does, so the summed loss ranks the
    penalties as the mean loss per held-out row does. A warning is logged when the
    chosen penalty is at either end of those tried, since the best one may then lie
    outside them (except below a penalty of 0).

    Each tuning fit is on the rows outside its fold, about half of them with two
    folds, and the chosen penalty is returned as is. Since it is per row, it shrinks
    a fit on all the rows about as much as the tuning fits, but the best penalty for
    more rows is usually somewhat smaller.

    :param ht: Held-out log loss per penalty
        (:func:`tune_regional_mutation_rate_penalty`).
    :param default_penalty: Penalty returned, with a warning, when ``ht`` is empty.
    :return: Chosen per-row ridge penalty.
    """
    rows = ht.aggregate(
        hl.agg.group_by(
            ht.penalty,
            hl.struct(
                log_loss=hl.agg.sum(ht.log_loss),
                scored=hl.agg.collect_as_set(
                    hl.struct(context=ht.context, fold=ht.fold, n_test=ht.n_test)
                ),
            ),
        )
    )
    if not rows:
        logger.warning(
            "No held-out rows to score the ridge penalties on; using %s.",
            default_penalty,
        )
        return default_penalty

    tried = sorted(rows)
    if any(rows[p].scored != rows[tried[0]].scored for p in tried):
        raise ValueError(
            "The ridge penalties weren't all scored on the same contexts, folds, and"
            " numbers of held-out rows, so their summed losses can't be compared."
        )
    loss = {p: rows[p].log_loss for p in tried}
    not_finite = [p for p in tried if not np.isfinite(loss[p])]
    if not_finite:
        raise ValueError(f"The held-out log loss isn't finite at {not_finite}.")
    # min keeps the first of tied penalties, so the smallest.
    penalty = min(tried, key=lambda p: loss[p])
    # Nothing lies below a penalty of 0.
    if len(tried) > 1 and (penalty == tried[-1] or 0 < penalty == tried[0]):
        logger.warning(
            "The chosen ridge penalty %s is at the end of those tried (%s).",
            penalty,
            tried,
        )

    return penalty


def _model_rate_expr(
    features_expr: hl.expr.StructExpression, model: hl.Struct
) -> hl.expr.Float64Expression:
    """
    Get a window's regional mutation rate under one context's model.

    The rate is ``expit(intercept + sum(coefs * (x - mean) / sd)) / ave``, and is
    missing if any of the model's features is missing, NaN, or infinite. The model's
    values go into
    the expression as constants, so it compiles to plain arithmetic. Applying the
    models through literals and :func:`_linear_predictor_expr` instead was several
    times slower and, with 32 models, made a class too large for the JVM.

    :param features_expr: Struct with the model's feature fields.
    :param model: Row of a models Table (:data:`MODELS_HT_SCHEMA`).
    :return: Regional mutation rate expression.
    """
    eta = model.intercept + sum(
        coef * (hl.float64(features_expr[f]) - mean) / sd
        for f, mean, sd, coef in zip(model.features, model.mean, model.sd, model.coefs)
    )

    # A NaN or infinite feature makes eta NaN or infinite; treat it as missing.
    return hl.bind(
        lambda eta: hl.or_missing(hl.is_finite(eta), hl.expit(eta) / model.ave), eta
    )


def regional_mutation_rate_expr(
    features_expr: hl.expr.StructExpression,
    models_ht: hl.Table,
    contexts: List[str],
    missing_feature_rate: Optional[float] = None,
    skip_non_converged: bool = False,
) -> hl.expr.DictExpression:
    """
    Get a window's regional mutation rate in each context.

    A window's rate in a context is ``expit(intercept + sum(coefs * (x - mean) /
    sd)) / ave`` of that context's model, as in Gnocchi. Contexts without a model get
    a rate of 1.0, and contexts whose model needs a feature the window is missing get
    ``missing_feature_rate``, which leaves their rates missing by default; a NaN or
    infinite feature value counts as missing.
    Models of contexts not in ``contexts`` are ignored. The rates are missing when
    every modeled context in ``contexts`` is missing one of its model's features, so
    filtering to defined rates keeps the windows the models cover. With one feature
    list for every context, that is any window missing a feature every model uses.

    Models that didn't converge (``converged`` is False) are applied as fit, with a
    warning naming them, unless ``skip_non_converged`` is set, which gives their
    contexts a rate of 1.0, as if they had no model.

    .. note::

        Because the non-mutated sites are a sample, the model is case-control, and
        the ratio of mutation rates between a window and the average window is
        closer to the odds ratio than to this ratio of probabilities. This form
        multiplies the log-scale regional effect by about ``1 - ave``, which pulls
        it toward 1, and caps the rate at ``1 / ave``. That is negligible when DNMs
        are a few percent of the training rows, but not when they are a large
        share, as for CpG contexts.

    :param features_expr: Struct with the models' feature fields, e.g. the row of a
        feature Table with one row per window.
    :param models_ht: Models Table (:func:`fit_regional_mutation_rate_models`). It
        has one small row per context and is collected.
    :param contexts: Contexts to include in the rates.
    :param missing_feature_rate: Rate of a context whose model needs a feature the
        window is missing, or None to leave it missing. Default is None: a missing
        rate is left out of the normalization factors and stays missing when the
        rates are normalized, so it can be filled afterwards, e.g. with 1.0 for no
        regional effect. A value given here is normalized like the other rates.
    :param skip_non_converged: Whether to give contexts whose model didn't converge a
        rate of 1.0 instead of applying the model. Default is False.
    :return: Dict of each context in ``contexts`` to the window's rate, missing when
        no modeled context has a rate. With ``missing_feature_rate`` None, the
        default, the rates of contexts missing a model feature are missing within the
        dict.
    """
    models = [m for m in models_ht.collect() if m.context in contexts]
    non_converged = sorted(m.context for m in models if not m.converged)
    if non_converged:
        logger.warning(
            "The models of %s didn't converge; %s.",
            ", ".join(non_converged),
            (
                "their contexts get a rate of 1.0"
                if skip_non_converged
                else "they are applied as fit"
            ),
        )
        if skip_non_converged:
            models = [m for m in models if m.converged]
    if not models:
        logger.warning("None of the contexts has a model, so every rate is missing.")
        return hl.missing(hl.tdict(hl.tstr, hl.tfloat64))

    def _rates(*model_rates: hl.expr.Float64Expression) -> hl.expr.DictExpression:
        rates = dict(zip([m.context for m in models], model_rates))
        if missing_feature_rate is not None:
            rates = {c: hl.or_else(r, missing_feature_rate) for c, r in rates.items()}

        # A context without a model has no regional effect.
        return hl.or_missing(
            hl.any(*[hl.is_defined(r) for r in model_rates]),
            hl.dict({c: rates.get(c, 1.0) for c in contexts}),
        )

    # Bind the models' rates so each is computed once per window.
    return hl.bind(_rates, *[_model_rate_expr(features_expr, m) for m in models])


def regional_mutation_rate_normalization_factors_agg_expr(
    context_expr: hl.expr.StringExpression,
    rate_expr: hl.expr.Float64Expression,
    weight_expr: hl.expr.Float64Expression,
) -> hl.expr.DictExpression:
    """
    Get an aggregation of the weighted mean regional mutation rate per context.

    Each context's factor is ``sum(weight * rate) / sum(weight)`` over the rows with
    a finite rate and a positive weight; dividing the rates by it makes their
    weighted mean 1. Contexts with no such rows are left out.

    :param context_expr: Each row's context.
    :param rate_expr: Each row's regional mutation rate, e.g. from a rate Table
        looked up by locus; missing outside a window.
    :param weight_expr: Each row's weight, e.g. its expected variant count.
    :return: Aggregation expression for a dict of context to normalization factor.
    """
    return hl.agg.filter(
        hl.or_else(hl.is_finite(rate_expr) & (weight_expr > 0), False),
        hl.agg.group_by(
            context_expr,
            hl.agg.sum(weight_expr * rate_expr) / hl.agg.sum(weight_expr),
        ),
    )


def normalize_regional_mutation_rate_expr(
    rate_expr: hl.expr.DictExpression,
    normalization_factors: Union[Dict[str, float], hl.expr.DictExpression],
) -> hl.expr.DictExpression:
    """
    Divide each context's regional mutation rate by its normalization factor.

    Contexts without a factor are left unchanged.

    :param rate_expr: Dict of context to rate, e.g. from
        :func:`regional_mutation_rate_expr`.
    :param normalization_factors: Dict, or dict expression, of context to factor,
        e.g. from :func:`regional_mutation_rate_normalization_factors_agg_expr`.
    :return: Dict of context to normalized rate.
    """
    if isinstance(normalization_factors, dict):
        normalization_factors = hl.literal(
            normalization_factors, hl.tdict(hl.tstr, hl.tfloat64)
        )

    return hl.dict(
        rate_expr.items().map(
            lambda kv: (kv[0], kv[1] / normalization_factors.get(kv[0], 1.0))
        )
    )
