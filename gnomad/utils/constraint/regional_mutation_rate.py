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
linear systems are solved on the driver. A step that raises the loss is halved.
"""

import logging
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
) -> hl.Table:
    """
    Prepare a training Table for fitting regional mutation rate models.

    Each row gets its 0/1 label, its context, and the values of its features in
    order: the same ``features`` for every context, or its context's features when
    ``features`` is a dict. ``features`` is recorded in the globals, so the fit uses
    the same ones. Rows with a missing value are kept;
    :func:`fit_regional_mutation_rate_models` leaves them out.

    :param ht: Training Table with the feature fields.
    :param features: Feature fields to model, which must be numeric or boolean and
        are converted to float64: one list for every context, or a dict of context to
        its own list.
    :param label_expr: 0/1 label of each row (1 for a site with a DNM), a number or a
        boolean. Default is ``ht.y``.
    :param context_expr: Context of each row. Default is ``ht.context``.
    :return: Table keyed as ``ht`` with only ``label`` (float64), ``context``, and
        ``x``, the float64 values of the row's features in order (missing, with a
        dict, for a context not in it). Its only global is ``features``.
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
    ht = ht.select(
        label=hl.float64(ht.y if label_expr is None else label_expr),
        context=context_expr,
        x=x,
    )

    return ht.select_globals(features=features_expr)


def _complete_training_rows(
    ht: hl.Table, fields: List[str]
) -> Tuple[hl.Table, Dict[str, List[str]]]:
    """
    Check a prepared training Table and drop the rows with a missing label or feature.

    NaN counts as defined in Hail, and one NaN or infinite feature value would make
    its context's model NaN, so the remaining rows' features must all be finite.

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
            " prepare_regional_mutation_rate_training_ht."
        )
    ht = ht.filter(
        hl.is_defined(ht.label) & hl.is_defined(ht.x) & hl.all(hl.is_defined, ht.x)
    )

    rows = ht.aggregate(
        hl.struct(
            n_not_finite=hl.agg.count_where(~hl.all(hl.is_finite, ht.x)),
            contexts=hl.agg.collect_as_set(ht.context),
        )
    )
    if rows.n_not_finite:
        raise ValueError(
            f"{rows.n_not_finite} training rows have a NaN or infinite feature value."
            " Make those values missing, or remove the rows."
        )

    features = hl.eval(ht.features)
    if isinstance(features, list):
        # With one list for every context, model each context the rows have.
        features = {c: features for c in sorted(rows.contexts - {None})}

    return ht, {c: list(f) for c, f in features.items()}


def _init_logistic_models(
    ht: hl.Table,
    group_features: Dict[str, List[str]],
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
            # The mean log loss at the last accepted coefficients (none yet).
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
    Take one step of Newton's method, with step halving, on a model's mean log loss.

    If the last step raised the mean log loss by more than the rounding in its sum,
    the step overshot, so the model goes back and takes half of it. Otherwise the
    next Newton step is taken from the current coefficients. The step is solved by
    least squares, so collinear features, which leave the Hessian singular, still
    give the smallest step that solves it. The model has converged when a step, full
    or halved, changes no coefficient by more than ``tol``.

    :param model: Model to update in place.
    :param sums: The model's loss, gradient, and Hessian sums (:func:`_newton_sums`).
    :param tol: Convergence tolerance.
    :return: None.
    """
    k = len(model["coefs"]) + 1
    beta = np.array([model["intercept"], *model["coefs"]])
    objective = sums.loss / model["n"]
    # Near the minimum, rounding in the mean loss can be larger than a step's real
    # change, so a smaller rise isn't an overshoot. Summing n rows rounds by up to
    # about n * eps of the sum, and each row's loss by about eps on its own.
    rounding = sys.float_info.epsilon * (model["n"] * model["objective"] + 1)
    model["halving"] = objective > model["objective"] + rounding
    if model["halving"]:
        model["step"] = model["step"] / 2
        model["converged"] = bool(np.max(np.abs(model["step"])) < tol)
    else:
        grad = -np.array(sums.g) / model["n"]
        hess = np.array(sums.h).reshape(k, k) / model["n"]
        model["objective"] = objective
        model["n_saturated"] = sums.n_saturated
        model["start"] = beta
        model["step"] = np.linalg.lstsq(hess, grad, rcond=None)[0]
        model["converged"] = bool(np.max(np.abs(model["step"])) < tol)
    beta = model["start"] - model["step"]
    model["intercept"] = float(beta[0])
    model["coefs"] = beta[1:].tolist()


def _fit_logistic_models(
    ht: hl.Table,
    group_features: Dict[str, List[str]],
    max_iter: int,
    tol: float,
    log: bool = True,
) -> Dict[str, Dict[str, Any]]:
    """
    Fit a logistic model for each group of rows of a prepared Table.

    All groups are fit together by Newton's method (:func:`_newton_update`), with one
    pass over ``ht`` per step. A model whose last step was halved, because it
    converged that way or ran out of steps, keeps its last accepted coefficients. A
    model with a fitted probability numerically 0 or 1 at those coefficients is
    marked not converged: its features separate the labels, at least in part, so the
    log loss has no minimum and the coefficients grow until the probabilities round
    off. R's glm.fit warns about the same fitted probabilities.

    :param ht: Complete training rows (:func:`_complete_training_rows`) with a
        ``group`` field.
    :param group_features: Features of each group, in the order of its rows' ``x``.
    :param max_iter: Maximum number of passes over ``ht`` for Newton steps,
        including halved ones.
    :param tol: A fit has converged when a step, full or halved, changes no
        coefficient by more than this.
    :param log: Whether to log groups without a model, features left out, and fits
        that don't converge. Default is True.
    :return: Dictionary of group to model, for the groups that have one.
    """
    models = _init_logistic_models(ht, group_features, log)
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
        # A model whose last step was halved sits at a point it never evaluated, so
        # it goes back to its last accepted coefficients, the best it knows.
        if model.get("halving"):
            model["intercept"] = float(model["start"][0])
            model["coefs"] = model["start"][1:].tolist()
        if model.get("n_saturated"):
            model["converged"] = False

    if log:
        for group, model in models.items():
            if model.get("n_saturated"):
                logger.warning(
                    "The model fit for %s did not converge: %d training rows have a"
                    " fitted probability of 0 or 1, so its features separate the"
                    " labels.",
                    group,
                    model["n_saturated"],
                )
            elif not model["converged"]:
                logger.warning("The model fit for %s did not converge.", group)
        logger.info("%d of %d contexts modeled.", len(models), len(group_features))

    return models


def fit_regional_mutation_rate_models(
    ht: hl.Table,
    max_iter: int = 50,
    tol: float = 1e-8,
) -> hl.Table:
    """
    Fit a logistic regional mutation rate model for each context.

    For each context, rows missing the label or any of its features are dropped, and
    features with a single value in the remaining rows are left out, since they carry
    no information and can't be standardized. The features are standardized with
    their mean and standard deviation (ddof=0), and the fit minimizes the mean log
    loss. A site's linear predictor is ``intercept + sum(coefs * (x - mean) / sd)``.

    The fit is by Newton's method, and each step is a pass over ``ht``, so ``ht``
    should be read from disk or checkpointed. A step that raises the loss is halved
    and tried again. A fit whose last step was halved, because the halved step fell
    below ``tol`` or ``max_iter`` ran out, keeps its last accepted coefficients.
    Each step is solved by least squares, so collinear features still fit: their
    coefficients aren't unique, but the model's predictions are. Features that
    separate the labels, even in part, leave the loss without a minimum, and the
    coefficients grow until some fitted probabilities are numerically 0 or 1; such a
    fit is marked not converged. Contexts with no model, features left out because
    they are constant, and fits that don't converge are logged. Feature values must
    be finite: a NaN or infinite value raises an error.

    :param ht: Training Table from
        :func:`prepare_regional_mutation_rate_training_ht`.
    :param max_iter: Maximum number of Newton steps, including halved ones. Default
        is 50.
    :param tol: The fit has converged when a step, full or halved, changes no
        coefficient by more than this. Default is 1e-8.
    :return: Table keyed by ``context`` with one row per modeled context: the kept
        ``features``, their training ``mean`` and ``sd``, ``coefs``, ``intercept``,
        ``ave`` (the predicted probability with every feature at its mean),
        ``n_mutated`` and ``n_background`` (rows of each label), and ``converged``
        (see :data:`MODELS_HT_SCHEMA`).
    """
    ht, context_features = _complete_training_rows(ht, ["label", "context", "x"])
    models = _fit_logistic_models(
        ht.annotate(group=ht.context), context_features, max_iter, tol
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
    )


def _model_rate_expr(
    features_expr: hl.expr.StructExpression, model: hl.Struct
) -> hl.expr.Float64Expression:
    """
    Get a window's regional mutation rate under one context's model.

    The rate is ``expit(intercept + sum(coefs * (x - mean) / sd)) / ave``, and is
    missing if any of the model's features is missing. The model's values go into
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

    return hl.expit(eta) / model.ave


def regional_mutation_rate_expr(
    features_expr: hl.expr.StructExpression,
    models_ht: hl.Table,
    contexts: List[str],
    missing_feature_rate: Optional[float] = 1.0,
) -> hl.expr.DictExpression:
    """
    Get a window's regional mutation rate in each context.

    A window's rate in a context is ``expit(intercept + sum(coefs * (x - mean) /
    sd)) / ave`` of that context's model, as in Gnocchi. Contexts without a model get
    a rate of 1.0, and contexts whose model needs a feature the window is missing get
    ``missing_feature_rate``. Models of contexts not in ``contexts`` are ignored. The
    rates are missing when no modeled context in ``contexts`` has one, e.g. for a
    window missing every feature, so filtering to defined rates keeps the windows the
    models cover. A NaN or infinite feature value gives a NaN rate.

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
        window is missing, or None to leave it missing. Default is 1.0, the model's
        rate for the average window.
    :return: Dict of each context in ``contexts`` to the window's rate, missing when
        no modeled context has a rate. With ``missing_feature_rate`` None, the rates
        of contexts missing a model feature are missing within the dict.
    """
    models = [m for m in models_ht.collect() if m.context in contexts]
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
    a rate and a positive weight; dividing the rates by it makes their weighted mean
    1. Contexts with no such rows are left out.

    :param context_expr: Each row's context.
    :param rate_expr: Each row's regional mutation rate, e.g. from a rate Table
        looked up by locus; missing outside a window.
    :param weight_expr: Each row's weight, e.g. its expected variant count.
    :return: Aggregation expression for a dict of context to normalization factor.
    """
    return hl.agg.filter(
        hl.is_defined(rate_expr) & hl.or_else(weight_expr > 0, False),
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
