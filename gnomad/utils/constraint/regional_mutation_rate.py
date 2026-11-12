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
    :param fold_expr: Fold of each row, e.g. odd or even chromosomes, needed for
        :func:`tune_regional_mutation_rate_penalty`. Default is None.
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

    features = hl.eval(ht.features)
    if isinstance(features, list):
        # With one list for every context, model each context the rows have.
        contexts = ht.aggregate(hl.agg.collect_as_set(ht.context))
        features = {c: features for c in sorted(contexts - {None})}

    return ht, {c: list(f) for c, f in features.items()}


def _init_ridge_logistic_models(
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
            # The penalized loss at the last accepted coefficients (none yet).
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

    :param ht: Complete training rows (:func:`_complete_training_rows`) with a
        ``group`` field.
    :param models: Models to evaluate, by group.
    :return: Dictionary of group to a struct with ``g``, ``h``, and ``loss``.
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
    Otherwise the next Newton step is taken from the current coefficients. The model
    has converged when a step changes no coefficient by more than ``tol``.

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
    # Near the minimum, the rounding in the summed loss, up to about n * eps of it,
    # can be larger than a step's real change, so a smaller rise isn't an overshoot.
    rounding = model["n"] * sys.float_info.epsilon * model["objective"]
    if objective > model["objective"] + rounding:
        model["step"] = model["step"] / 2
        model["converged"] = bool(np.max(np.abs(model["step"])) < tol)
    else:
        grad = -np.array(sums.g) / model["n"] + penalty * beta
        hess = np.array(sums.h).reshape(k, k) / model["n"] + np.diag(penalty)
        model["objective"] = objective
        model["start"] = beta
        model["step"] = np.linalg.solve(hess, grad)
        model["converged"] = bool(np.max(np.abs(model["step"])) < tol)
    beta = model["start"] - model["step"]
    model["intercept"] = float(beta[0])
    model["coefs"] = beta[1:].tolist()


def _fit_ridge_logistic_models(
    ht: hl.Table,
    group_features: Dict[str, List[str]],
    ridge_penalties: Dict[str, float],
    max_iter: int,
    tol: float,
    log: bool = True,
) -> Dict[str, Dict[str, Any]]:
    """
    Fit a ridge-penalized logistic model for each group of rows of a prepared Table.

    A group is a context, or, when tuning, a context in one penalty and fold. All
    groups are fit together by Newton's method (:func:`_newton_update`), with one
    pass over ``ht`` per step.

    :param ht: Complete training rows (:func:`_complete_training_rows`) with a
        ``group`` field.
    :param group_features: Features of each group, in the order of its rows' ``x``.
    :param ridge_penalties: Per-row ridge penalty of each group, greater than 0.
    :param max_iter: Maximum number of Newton steps.
    :param tol: A fit has converged when a Newton step changes no coefficient by more
        than this.
    :param log: Whether to log groups without a model, features left out, and fits
        that don't converge. Default is True.
    :return: Dictionary of group to model, for the groups that have one.
    """
    models = _init_ridge_logistic_models(ht, group_features, ridge_penalties, log)
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

    if log:
        for group, model in models.items():
            if not model["converged"]:
                logger.warning("The model fit for %s did not converge.", group)
        logger.info("%d of %d contexts modeled.", len(models), len(group_features))

    return models


def _check_penalties(penalties: List[float]) -> None:
    """
    Check that every ridge penalty is greater than 0.

    Without a penalty the objective can have no minimum, e.g. when a feature
    separates the labels, and collinear features leave it without a unique one.

    :param penalties: Per-row ridge penalties.
    :return: None.
    """
    if not penalties or min(penalties) <= 0:
        raise ValueError(f"Ridge penalties must be greater than 0, not {penalties}.")


def fit_regional_mutation_rate_models(
    ht: hl.Table,
    ridge_penalty: float,
    max_iter: int = 50,
    tol: float = 1e-8,
) -> hl.Table:
    """
    Fit a ridge-penalized logistic regional mutation rate model for each context.

    For each context, rows missing the label or any of its features are dropped, and
    features with a single value in the remaining rows are left out, since they carry
    no information and can't be standardized. The features are standardized with
    their mean and standard deviation (ddof=0), and the fit minimizes mean log loss +
    ``ridge_penalty`` / 2 * sum(coef ** 2), with the intercept unpenalized. A site's
    linear predictor is ``intercept + sum(coefs * (x - mean) / sd)``.

    The fit is by Newton's method, and each step is a pass over ``ht``, so ``ht``
    should be read from disk or checkpointed. Contexts with no model, features left
    out because they are constant, and fits that don't converge are logged.

    :param ht: Training Table from
        :func:`prepare_regional_mutation_rate_training_ht`.
    :param ridge_penalty: Per-row ridge penalty, greater than 0.
    :param max_iter: Maximum number of Newton steps. Default is 50.
    :param tol: The fit has converged when a Newton step changes no coefficient by
        more than this. Default is 1e-8.
    :return: Table keyed by ``context`` with one row per modeled context: the kept
        ``features``, their training ``mean`` and ``sd``, ``coefs``, ``intercept``,
        ``ave`` (the predicted probability with every feature at its mean),
        ``n_mutated`` and ``n_background`` (rows of each label), and ``converged``
        (see :data:`MODELS_HT_SCHEMA`).
    """
    _check_penalties([ridge_penalty])
    ht, context_features = _complete_training_rows(ht, ["label", "context", "x"])
    models = _fit_ridge_logistic_models(
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
    should be read from disk or checkpointed. Fits that don't converge are counted
    in a warning and still scored.

    :param ht: Training Table from
        :func:`prepare_regional_mutation_rate_training_ht`, with a ``fold``.
    :param penalties: Per-row ridge penalties to try, each greater than 0. Repeats
        are tried once.
    :param max_iter: Maximum number of Newton steps per fit. Default is 50.
    :param tol: Convergence tolerance of each fit. Default is 1e-8.
    :return: Table keyed by ``penalty`` (float64), ``context`` (str), and ``fold``
        (str), with the number of held-out rows with no missing label or feature
        (``n_test``, int64) and their summed log loss (``log_loss``, float64). A
        context and fold with no model (its training rows have only one label, or no
        feature varies) has no rows, and the Table is empty when no held-out rows
        could be scored.
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
    models = _fit_ridge_logistic_models(
        _rows_by_fit(ht, fits, held_out=False),
        {g: group.features for g, group in groups.items()},
        {g: group.penalty for g, group in groups.items()},
        max_iter,
        tol,
        log=False,
    )
    not_converged = sum(not m["converged"] for m in models.values())
    if not_converged:
        logger.warning(
            "%d of %d tuning fits did not converge; they are still scored.",
            not_converged,
            len(models),
        )
    scores = _held_out_log_loss(_rows_by_fit(ht, fits, held_out=True), models)

    return hl.Table.parallelize(
        [
            {
                "penalty": groups[g].penalty,
                "context": groups[g].context,
                "fold": groups[g].fold,
                **score,
            }
            for g, score in scores.items()
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


def choose_regional_mutation_rate_penalty(
    ht: hl.Table, default_penalty: float
) -> float:
    """
    Choose the ridge penalty with the lowest held-out log loss.

    The log loss is summed over contexts and folds, so one penalty is chosen for
    every context; on a tie, the smaller penalty is chosen. A warning is logged when
    the chosen penalty is at either end of those tried, since the best one may then
    lie outside them.

    :param ht: Held-out log loss per penalty
        (:func:`tune_regional_mutation_rate_penalty`).
    :param default_penalty: Penalty returned, with a warning, when ``ht`` is empty.
    :return: Chosen per-row ridge penalty.
    """
    loss = ht.aggregate(hl.agg.group_by(ht.penalty, hl.agg.sum(ht.log_loss)))
    if not loss:
        logger.warning(
            "No held-out rows to score the ridge penalties on; using %s.",
            default_penalty,
        )
        return default_penalty

    tried = sorted(loss)
    # min keeps the first of tied penalties, so the smallest.
    penalty = min(tried, key=lambda p: loss[p])
    if len(tried) > 1 and penalty in (tried[0], tried[-1]):
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
    missing if any of the model's features is missing. The model's values go into
    the expression as constants, so it compiles to plain arithmetic, which is about
    4 times faster than applying the models with :func:`_linear_predictor_expr`.

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
    models cover.

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
    :return: Dict of each context in ``contexts`` to the window's rate, or missing.
    """
    models = [m for m in models_ht.collect() if m.context in contexts]
    if not models:
        logger.warning("None of the contexts has a model, so every rate is missing.")
        return hl.missing(hl.tdict(hl.tstr, hl.tfloat64))

    def _rates(model_rates: hl.expr.DictExpression) -> hl.expr.DictExpression:
        def _rate(context: hl.expr.StringExpression) -> hl.expr.Float64Expression:
            rate = model_rates.get(context)
            if missing_feature_rate is not None:
                rate = hl.or_else(rate, missing_feature_rate)
            # A context without a model has no regional effect.
            return hl.if_else(model_rates.contains(context), rate, 1.0)

        return hl.or_missing(
            hl.any(hl.is_defined, model_rates.values()),
            hl.dict(hl.literal(contexts).map(lambda c: (c, _rate(c)))),
        )

    # Bind the models' rates so each is computed once per window.
    return hl.bind(
        _rates,
        hl.dict([(m.context, _model_rate_expr(features_expr, m)) for m in models]),
    )


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
