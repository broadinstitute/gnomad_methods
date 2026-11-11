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
gradient and Hessian over the training Table, and only the small per-context linear
systems are solved on the driver.
"""

import logging
from typing import Any, Dict, List, Optional, Union

import hail as hl
import numpy as np

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_mutation_rate")
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
"""Row type of the models Table returned by :func:`fit_regional_mutation_rate_models`."""

_MODEL_EXPR_TYPE = hl.tstruct(
    positions=hl.tarray(hl.tint32),
    mean=hl.tarray(hl.tfloat64),
    sd=hl.tarray(hl.tfloat64),
    coefs=hl.tarray(hl.tfloat64),
    intercept=hl.tfloat64,
)


def _models_expr(models: Dict[str, Dict[str, Any]]) -> hl.expr.DictExpression:
    """
    Get a dict expression of context to model, for evaluating models on rows.

    :param models: Dictionary of context to model, each with ``positions`` (the
        positions of its features in a row's feature array), ``mean``, ``sd``,
        ``coefs``, and ``intercept``.
    :return: Dict expression of context to model struct.
    """
    return hl.literal(
        {
            c: hl.Struct(**{f: m[f] for f in _MODEL_EXPR_TYPE})
            for c, m in models.items()
        },
        hl.tdict(hl.tstr, _MODEL_EXPR_TYPE),
    )


def _linear_predictor_expr(
    x_expr: hl.expr.ArrayNumericExpression, model_expr: hl.expr.StructExpression
) -> hl.expr.Float64Expression:
    """
    Get the linear predictor ``intercept + sum(coefs * (x - mean) / sd)`` of a model.

    It is missing if any of the model's features is missing.

    :param x_expr: Feature values, indexed by the model's ``positions``.
    :param model_expr: Model struct (see :func:`_models_expr`).
    :return: Linear predictor expression.
    """
    return model_expr.intercept + hl.sum(
        hl.range(hl.len(model_expr.positions)).map(
            lambda i: model_expr.coefs[i]
            * (x_expr[model_expr.positions[i]] - model_expr.mean[i])
            / model_expr.sd[i]
        ),
        filter_missing=False,
    )


def _prepare_training_ht(
    ht: hl.Table,
    context_features: Dict[str, List[str]],
    label_expr: hl.expr.NumericExpression,
    context_expr: hl.expr.StringExpression,
    fold_expr: Optional[hl.expr.StringExpression] = None,
) -> hl.Table:
    """
    Select each training row's label, context, fold, and its context's features.

    Rows with a missing label, a context not in ``context_features``, or a missing
    value for one of the context's features are dropped.

    :param ht: Training Table.
    :param context_features: Features to model for each context.
    :param label_expr: 0/1 label of each row (1 for a site with a DNM).
    :param context_expr: Context of each row.
    :param fold_expr: Optional fold of each row. Default is None.
    :return: Table with ``label`` (float64), ``context``, ``fold`` (if given), and
        ``x``, the values of the row's context's features in order.
    """
    all_features = list(
        dict.fromkeys(f for features in context_features.values() for f in features)
    )
    positions = hl.literal(
        {c: [all_features.index(f) for f in fs] for c, fs in context_features.items()},
        hl.tdict(hl.tstr, hl.tarray(hl.tint32)),
    )
    x = hl.array([hl.float64(ht[f]) for f in all_features])
    fields = {"label": hl.float64(label_expr), "context": context_expr}
    if fold_expr is not None:
        fields["fold"] = fold_expr
    ht = ht.select(**fields, x=positions.get(context_expr).map(lambda i: x[i]))

    return ht.filter(
        hl.is_defined(ht.label)
        & hl.is_defined(ht.x)
        & hl.all(lambda v: hl.is_defined(v), ht.x)
    )


def _fit_models(
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
    groups are fit together, with one pass over ``ht`` per Newton step.

    :param ht: Table from :func:`_prepare_training_ht` with a ``group`` field.
    :param group_features: Features of each group, in the order of its rows' ``x``.
    :param ridge_penalties: Per-row ridge penalty of each group, 0 or more.
    :param max_iter: Maximum number of Newton steps.
    :param tol: A fit has converged when no coefficient changes by more than this.
    :param log: Whether to log groups without a model, features left out, and fits
        that don't converge. Default is True.
    :return: Dictionary of group to model, for the groups that have one.
    """
    # Features with a single value carry no information and can't be standardized.
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
        n_mutated = int(round(s.n_mutated)) if s else 0
        keep = [i for i, v in enumerate(s.x) if v.min < v.max] if s else []
        if not n_mutated or n_mutated == s.n or not keep:
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
            "n": s.n,
            "n_mutated": n_mutated,
            "intercept": float(np.log(n_mutated / (s.n - n_mutated))),
            "coefs": [0.0] * len(keep),
            "converged": False,
            "penalty": ridge_penalties[group],
        }

    if not models:
        return models

    # Standard deviation (ddof=0) in a second pass, around the means, for precision.
    means = hl.literal(
        {
            c: hl.Struct(positions=m["positions"], mean=m["mean"])
            for c, m in models.items()
        },
        hl.tdict(
            hl.tstr,
            hl.tstruct(positions=hl.tarray(hl.tint32), mean=hl.tarray(hl.tfloat64)),
        ),
    )
    m_expr = means.get(ht.group)
    ss = ht.aggregate(
        hl.agg.filter(
            hl.is_defined(m_expr),
            hl.agg.group_by(
                ht.group,
                hl.agg.array_sum(
                    hl.range(hl.len(m_expr.positions)).map(
                        lambda i: (ht.x[m_expr.positions[i]] - m_expr.mean[i]) ** 2
                    )
                ),
            ),
        )
    )
    for group, m in models.items():
        m["sd"] = [float(np.sqrt(v / m["n"])) for v in ss[group]]

    # Newton's method on mean log loss + ridge_penalty / 2 * sum(coef ** 2), with the
    # intercept unpenalized, on the standardized features.
    for _ in range(max_iter):
        todo = {c: m for c, m in models.items() if not m["converged"]}
        if not todo:
            break
        model_expr = _models_expr(todo).get(ht.group)
        z = hl.range(hl.len(model_expr.positions)).map(
            lambda i: (ht.x[model_expr.positions[i]] - model_expr.mean[i])
            / model_expr.sd[i]
        )
        p = hl.expit(_linear_predictor_expr(ht.x, model_expr))
        v = hl.array([1.0]).extend(z)
        sums = ht.aggregate(
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
                    ),
                ),
            )
        )
        for group, m in todo.items():
            k = len(m["coefs"]) + 1
            beta = np.array([m["intercept"], *m["coefs"]])
            penalty = m["penalty"] * np.r_[0.0, np.ones(k - 1)]
            grad = -np.array(sums[group].g) / m["n"] + penalty * beta
            hess = np.array(sums[group].h).reshape(k, k) / m["n"] + np.diag(penalty)
            step = np.linalg.lstsq(hess, grad, rcond=None)[0]
            beta = beta - step
            m["intercept"] = float(beta[0])
            m["coefs"] = beta[1:].tolist()
            m["converged"] = bool(np.max(np.abs(step)) < tol)

    if log:
        for group, m in models.items():
            if not m["converged"]:
                logger.warning("The model fit for %s did not converge.", group)
        logger.info("%d of %d contexts modeled.", len(models), len(group_features))

    return models


def fit_regional_mutation_rate_models(
    ht: hl.Table,
    context_features: Dict[str, List[str]],
    ridge_penalty: float,
    label_expr: Optional[hl.expr.NumericExpression] = None,
    context_expr: Optional[hl.expr.StringExpression] = None,
    max_iter: int = 50,
    tol: float = 1e-8,
) -> hl.Table:
    """
    Fit a ridge-penalized logistic regional mutation rate model for each context.

    For each context, rows missing any of its features are dropped, and features with
    a single value in the remaining rows are left out, since they carry no
    information and can't be standardized. The features are z-scored (ddof=0), and
    the fit minimizes mean log loss + ``ridge_penalty`` / 2 * sum(coef ** 2), with the
    intercept unpenalized; a penalty of 0 gives an unpenalized fit. A site's linear
    predictor is ``intercept + sum(coefs * (x - mean) / sd)``.

    Contexts with no model, features left out because they are constant, and fits
    that don't converge are logged.

    :param ht: Training Table with the feature fields.
    :param context_features: Feature fields to model for each context.
    :param ridge_penalty: Per-row ridge penalty, 0 or more.
    :param label_expr: 0/1 label of each row (1 for a site with a DNM). Default is
        ``ht.y``.
    :param context_expr: Context of each row. Default is ``ht.context``.
    :param max_iter: Maximum number of Newton steps. Default is 50.
    :param tol: The fit has converged when no coefficient changes by more than this
        in a step. Default is 1e-8.
    :return: Table keyed by ``context`` with one row per modeled context: the kept
        ``features``, their training ``mean`` and ``sd``, ``coefs``, ``intercept``,
        ``ave`` (the predicted probability with every feature at its mean),
        ``n_mutated`` and ``n_background`` (rows of each label), and ``converged``.
    """
    ht = _prepare_training_ht(
        ht,
        context_features,
        ht.y if label_expr is None else label_expr,
        ht.context if context_expr is None else context_expr,
    )
    models = _fit_models(
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
                **{f: m[f] for f in ["features", "mean", "sd", "coefs", "intercept"]},
                "ave": float(1.0 / (1.0 + np.exp(-m["intercept"]))),
                "n_mutated": m["n_mutated"],
                "n_background": m["n"] - m["n_mutated"],
                "converged": m["converged"],
            }
            for context, m in models.items()
        ],
        MODELS_HT_SCHEMA,
        key="context",
    )


def tune_regional_mutation_rate_penalty(
    ht: hl.Table,
    context_features: Dict[str, List[str]],
    penalties: List[float],
    label_expr: Optional[hl.expr.NumericExpression] = None,
    context_expr: Optional[hl.expr.StringExpression] = None,
    fold_expr: Optional[hl.expr.StringExpression] = None,
    max_iter: int = 50,
    tol: float = 1e-8,
) -> hl.Table:
    """
    Get the held-out log loss of regional mutation rate models for each ridge penalty.

    For each penalty, context, and fold, the context's model is fit on the rows
    outside the fold (:func:`fit_regional_mutation_rate_models`) and scored on the
    rows in it. Rows with a missing fold are only used for fitting.

    :param ht: Training Table with the feature fields.
    :param context_features: Feature fields to model for each context.
    :param penalties: Per-row ridge penalties to try.
    :param label_expr: 0/1 label of each row. Default is ``ht.y``.
    :param context_expr: Context of each row. Default is ``ht.context``.
    :param fold_expr: Fold of each row, e.g. odd or even chromosomes. Default is
        ``ht.fold``.
    :param max_iter: Maximum number of Newton steps per fit. Default is 50.
    :param tol: Convergence tolerance of each fit. Default is 1e-8.
    :return: Table keyed by ``penalty``, ``context``, and ``fold``, with the number of
        held-out rows with no missing feature (``n_test``) and their summed log loss
        (``log_loss``). It is empty when no held-out rows could be scored.
    """
    ht = _prepare_training_ht(
        ht,
        context_features,
        ht.y if label_expr is None else label_expr,
        ht.context if context_expr is None else context_expr,
        ht.fold if fold_expr is None else fold_expr,
    )
    folds = sorted(ht.aggregate(hl.agg.collect_as_set(ht.fold)) - {None})
    # Fit i is penalty fits[i][0] on the rows outside fold fits[i][1]. Every fit is
    # done at once, each row joining the fits it belongs to, and so is the scoring.
    fits = [(float(p), f) for p in penalties for f in folds]
    fit_folds = hl.literal([f for _, f in fits], hl.tarray(hl.tstr))

    def _groups_ht(in_fold: bool) -> hl.Table:
        """Return the rows of every fit's training (or held-out) set, by group."""
        fit_ids = hl.range(len(fits)).filter(
            lambda i: (
                (fit_folds[i] == ht.fold)
                if in_fold
                else hl.or_else(fit_folds[i] != ht.fold, True)
            )
        )
        groups_ht = ht.annotate(fit=fit_ids).explode("fit")
        return groups_ht.annotate(group=hl.str(groups_ht.fit) + "|" + groups_ht.context)

    groups = {
        f"{i}|{c}": (penalty, c, features)
        for i, (penalty, _) in enumerate(fits)
        for c, features in context_features.items()
    }
    models = _fit_models(
        _groups_ht(in_fold=False),
        {g: features for g, (_, _, features) in groups.items()},
        {g: penalty for g, (penalty, _, _) in groups.items()},
        max_iter,
        tol,
        log=False,
    )

    test_ht = _groups_ht(in_fold=True)
    model_expr = _models_expr(models).get(test_ht.group)
    eta = _linear_predictor_expr(test_ht.x, model_expr)
    # log(1 + exp(eta)) - label * eta, written to avoid overflow.
    loss = hl.log(1 + hl.exp(-hl.abs(eta))) + hl.max(eta, 0) - test_ht.label * eta
    scores = test_ht.aggregate(
        hl.agg.filter(
            hl.is_defined(model_expr),
            hl.agg.group_by(
                test_ht.group,
                hl.struct(n_test=hl.agg.count(), log_loss=hl.agg.sum(loss)),
            ),
        )
    )
    rows = [
        {
            "penalty": groups[g][0],
            "context": groups[g][1],
            "fold": fits[int(g.split("|")[0])][1],
            **score,
        }
        for g, score in scores.items()
    ]

    return hl.Table.parallelize(
        rows,
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
    tuning_ht: hl.Table, default_penalty: float
) -> float:
    """
    Choose the ridge penalty with the lowest held-out log loss.

    The log loss is summed over contexts and folds, so one penalty is chosen for
    every context; on a tie, the smaller penalty is chosen. A warning is logged when
    the chosen penalty is at either end of those tried, since the best one may then
    lie outside them.

    :param tuning_ht: Held-out log loss per penalty
        (:func:`tune_regional_mutation_rate_penalty`).
    :param default_penalty: Penalty returned, with a warning, when ``tuning_ht`` is
        empty.
    :return: Chosen per-row ridge penalty.
    """
    loss = tuning_ht.aggregate(
        hl.agg.group_by(tuning_ht.penalty, hl.agg.sum(tuning_ht.log_loss))
    )
    if not loss:
        logger.warning(
            "No held-out rows to score the ridge penalties on; using %s.",
            default_penalty,
        )
        return default_penalty

    tried = sorted(loss)
    penalty = min(tried, key=lambda p: loss[p])
    if len(tried) > 1 and penalty in (tried[0], tried[-1]):
        logger.warning(
            "The chosen ridge penalty %s is at the end of those tried (%s).",
            penalty,
            tried,
        )

    return penalty


def compute_regional_mutation_rate_ht(
    features_ht: hl.Table,
    models_ht: hl.Table,
    contexts: List[str],
) -> hl.Table:
    """
    Compute the regional mutation rate of every window for each context.

    A window's rate in a context is ``expit(intercept + sum(coefs * (x - mean) /
    sd)) / ave`` of that context's model, as in Gnocchi. Contexts without a model,
    and windows missing one of a context's features, get a rate of 1.0. Models of
    contexts not in ``contexts`` are ignored. Windows are kept only if at least one
    modeled context in ``contexts`` has a rate.

    .. note::

        Because the non-mutated sites are a sample, the model is case-control, and
        the ratio of mutation rates between a window and the average window is
        closer to the odds ratio than to this ratio of probabilities. This form
        multiplies the log-scale regional effect by about ``1 - ave``, which pulls
        it toward 1, and caps the rate at ``1 / ave``. That is negligible when DNMs
        are a few percent of the training rows, but not when they are a large
        share, as for CpG contexts.

    :param features_ht: Feature Table with one row per window, e.g. keyed by an
        interval.
    :param models_ht: Models Table (:func:`fit_regional_mutation_rate_models`).
    :param contexts: Contexts to include in each window's rates.
    :return: ``features_ht`` keyed as before, with only an ``adj_r`` dict field
        mapping each context to its rate.
    """
    # The models Table has one small row per context.
    models = [m for m in models_ht.collect() if m.context in contexts]
    if not models:
        logger.warning("No context has a model, so the rate Table is empty.")
        return features_ht.select(adj_r=hl.empty_dict(hl.tstr, hl.tfloat64)).filter(
            False
        )

    all_features = list(dict.fromkeys(f for m in models for f in m.features))
    x = hl.array([hl.float64(features_ht[f]) for f in all_features])
    models_expr = _models_expr(
        {
            m.context: {
                **m,
                "positions": [all_features.index(f) for f in m.features],
            }
            for m in models
        }
    )
    aves = hl.literal({m.context: m.ave for m in models})
    rates = hl.dict(
        hl.literal([m.context for m in models]).map(
            lambda c: (c, hl.expit(_linear_predictor_expr(x, models_expr[c])) / aves[c])
        )
    )
    ht = features_ht.select(_rate=rates)
    ht = ht.filter(hl.any(lambda r: hl.is_defined(r), ht._rate.values()))

    # TODO: Windows missing a context's features get 1.0 here, the model's rate for
    #  the average window. Normalization then divides it by the context's factor, so
    #  these windows end up at 1 / factor instead of 1, unlike sites without a
    #  window. Consider leaving them missing here and filling 1.0 after
    #  normalization.
    return ht.select(
        adj_r=hl.dict(
            hl.literal(contexts).map(lambda c: (c, hl.or_else(ht._rate.get(c), 1.0)))
        )
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

    :param rate_expr: Dict of context to rate, e.g. the ``adj_r`` field of
        :func:`compute_regional_mutation_rate_ht`.
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
