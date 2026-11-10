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

Unlike the rest of :mod:`gnomad.utils.constraint`, the models are fit on the driver,
with scikit-learn on a pandas DataFrame of the collected training sites.
"""

import logging
from typing import Any, Dict, List, Optional, Union

import hail as hl
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("regional_mutation_rate")
logger.setLevel(logging.INFO)


def fit_regional_mutation_rate_model(
    df: pd.DataFrame,
    features: List[str],
    ridge_penalty: float,
    label_col: str = "y",
    max_iter: int = 10000,
) -> Optional[Dict[str, Any]]:
    """
    Fit a ridge-penalized logistic regional mutation rate model.

    Rows with any missing feature are dropped, and features with a single value in
    the remaining rows are left out, since they carry no information and can't be
    standardized. The features are z-scored, and the fit minimizes mean log loss +
    ``ridge_penalty`` / 2 * sum(coef ** 2), with the intercept unpenalized; a penalty
    of 0 gives an unpenalized fit. A site's linear predictor is
    ``intercept + sum(coefs * (x - mean) / sd)``.

    :param df: Training rows, with a 0/1 ``label_col`` (1 for a site with a DNM) and
        the feature columns.
    :param features: Feature columns to model.
    :param ridge_penalty: Per-row ridge penalty, 0 or more.
    :param label_col: Column with the 0/1 labels. Default is "y".
    :param max_iter: Maximum number of solver iterations. Default is 10000.
    :return: Dictionary with the kept ``features``, their training ``mean`` and
        ``sd`` (ddof=0), ``coefs``, ``intercept``, ``ave`` (the predicted
        probability with every feature at its mean), ``n_mutated`` and
        ``n_background`` (rows of each label), and ``converged`` (whether the solver
        stopped before ``max_iter`` iterations); or None if either label has no rows
        or no feature varies.
    """
    df = df[[label_col, *features]].dropna()
    labels = df[label_col].to_numpy(dtype=np.int64)
    n_mutated = int(labels.sum())
    n_background = len(labels) - n_mutated
    # Count distinct values rather than testing sd > 0: rounding can leave the sd of
    # a constant column at about 1e-17, which would make z explode at prediction.
    features = [f for f in features if df[f].nunique() > 1]
    if not n_mutated or not n_background or not features:
        return None

    x = df[features].to_numpy(dtype=np.float64)
    mean = x.mean(axis=0)
    sd = x.std(axis=0)
    # sklearn minimizes C * sum(log loss) + sum(coef ** 2) / 2, so C = 1 / (n *
    # penalty) gives the per-row penalty, and C = inf no penalty.
    c = np.inf if ridge_penalty == 0 else 1.0 / (len(labels) * ridge_penalty)
    fit = LogisticRegression(C=c, max_iter=max_iter)
    fit.fit((x - mean) / sd, labels)
    intercept = float(fit.intercept_[0])

    return {
        "features": features,
        "mean": mean.tolist(),
        "sd": sd.tolist(),
        "coefs": fit.coef_[0].tolist(),
        "intercept": intercept,
        "ave": float(1.0 / (1.0 + np.exp(-intercept))),
        "n_mutated": n_mutated,
        "n_background": n_background,
        "converged": bool(fit.n_iter_[0] < max_iter),
    }


def fit_regional_mutation_rate_models(
    df: pd.DataFrame,
    context_features: Dict[str, List[str]],
    ridge_penalty: float,
    context_col: str = "context",
    label_col: str = "y",
) -> List[Dict[str, Any]]:
    """
    Fit a regional mutation rate model for every context.

    Contexts with no model, features left out because they are constant, and fits
    that don't converge are logged.

    :param df: Training rows with ``context_col``, ``label_col``, and the feature
        columns.
    :param context_features: Features to model for each context.
    :param ridge_penalty: Per-row ridge penalty.
    :param context_col: Column with each row's context. Default is "context".
    :param label_col: Column with the 0/1 labels. Default is "y".
    :return: Model dictionaries (:func:`fit_regional_mutation_rate_model`) with a
        ``context`` field, for the contexts that have a model.
    """
    models = []
    for context, features in context_features.items():
        model = fit_regional_mutation_rate_model(
            df[df[context_col] == context], features, ridge_penalty, label_col
        )
        if model is None:
            logger.info(
                "No model for %s: no rows of one label, or no feature varies.",
                context,
            )
            continue
        constant = [f for f in features if f not in model["features"]]
        if constant:
            logger.info("%s: features left out (constant): %s", context, constant)
        if not model["converged"]:
            logger.warning("The model fit for %s did not converge.", context)
        models.append({"context": context, **model})
    logger.info("%d of %d contexts modeled.", len(models), len(context_features))

    return models


def predict_regional_mutation_rate_probability(
    model: Dict[str, Any], df: pd.DataFrame
) -> np.ndarray:
    """
    Get the predicted DNM probability of each row from a regional mutation rate model.

    :param model: Model dictionary (:func:`fit_regional_mutation_rate_model`).
    :param df: Rows with the model's feature columns.
    :return: Predicted probabilities.
    """
    z = (df[model["features"]].to_numpy(dtype=np.float64) - model["mean"]) / model["sd"]
    eta = model["intercept"] + z @ np.asarray(model["coefs"])

    return 1.0 / (1.0 + np.exp(-eta))


def tune_regional_mutation_rate_penalty(
    df: pd.DataFrame,
    context_features: Dict[str, List[str]],
    penalties: List[float],
    fold_col: str = "fold",
    context_col: str = "context",
    label_col: str = "y",
) -> pd.DataFrame:
    """
    Get the held-out log loss of regional mutation rate models for each ridge penalty.

    For each penalty, context, and fold, the context's model is fit on the rows
    outside the fold and scored on the rows in it. Rows with a missing fold are only
    used for fitting.

    :param df: Training rows with ``fold_col``, ``context_col``, ``label_col``, and
        the feature columns.
    :param context_features: Features to model for each context.
    :param penalties: Per-row ridge penalties to try.
    :param fold_col: Column with each row's fold, e.g. odd or even chromosomes.
        Default is "fold".
    :param context_col: Column with each row's context. Default is "context".
    :param label_col: Column with the 0/1 labels. Default is "y".
    :return: DataFrame with ``penalty``, ``context``, ``fold``, the number of held-out
        rows (``n_test``), and their summed log loss (``log_loss``). It is empty when
        no held-out rows could be scored.
    """
    folds = sorted(df[fold_col].dropna().unique())
    rows = []
    for penalty in penalties:
        for context, features in context_features.items():
            context_df = df[df[context_col] == context]
            for fold in folds:
                # A missing fold is not in any fold, whatever the column's dtype.
                in_fold = context_df[fold_col].eq(fold).fillna(False).to_numpy(bool)
                model = fit_regional_mutation_rate_model(
                    context_df[~in_fold], features, penalty, label_col
                )
                test_df = context_df[in_fold][[label_col, *features]].dropna()
                if model is None or test_df.empty:
                    continue
                rows.append(
                    {
                        "penalty": float(penalty),
                        "context": context,
                        "fold": fold,
                        "n_test": len(test_df),
                        "log_loss": float(
                            log_loss(
                                test_df[label_col],
                                predict_regional_mutation_rate_probability(
                                    model, test_df
                                ),
                                labels=[0, 1],
                                normalize=False,
                            )
                        ),
                    }
                )

    return pd.DataFrame(
        rows, columns=["penalty", "context", "fold", "n_test", "log_loss"]
    )


def choose_regional_mutation_rate_penalty(
    tuning_df: pd.DataFrame, default_penalty: float
) -> float:
    """
    Choose the ridge penalty with the lowest held-out log loss.

    The log loss is summed over contexts and folds, so one penalty is chosen for
    every context. A warning is logged when the chosen penalty is at either end of
    those tried, since the best one may then lie outside them.

    :param tuning_df: Held-out log loss per penalty
        (:func:`tune_regional_mutation_rate_penalty`).
    :param default_penalty: Penalty returned, with a warning, when ``tuning_df`` is
        empty.
    :return: Chosen per-row ridge penalty.
    """
    if tuning_df.empty:
        logger.warning(
            "No held-out rows to score the ridge penalties on; using %s.",
            default_penalty,
        )
        return default_penalty

    penalty = float(tuning_df.groupby("penalty")["log_loss"].sum().idxmin())
    tried = sorted(tuning_df["penalty"].unique())
    if len(tried) > 1 and penalty in (tried[0], tried[-1]):
        logger.warning(
            "The chosen ridge penalty %s is at the end of those tried (%s).",
            penalty,
            tried,
        )

    return penalty


def regional_mutation_rate_models_to_ht(models: List[Dict[str, Any]]) -> hl.Table:
    """
    Convert regional mutation rate model dictionaries to a Table keyed by context.

    :param models: Model dictionaries (:func:`fit_regional_mutation_rate_models`).
    :return: Table keyed by ``context`` with one row per model.
    """
    schema = hl.tstruct(
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

    return hl.Table.parallelize(models, schema=schema).key_by("context")


def regional_mutation_rate_expr(
    features_expr: hl.expr.StructExpression,
    model: Union[Dict[str, Any], hl.Struct],
) -> hl.expr.Float64Expression:
    """
    Get the regional mutation rate of a window from one context's model.

    The rate is ``expit(intercept + sum(coefs * (x - mean) / sd)) / ave``, as in
    Gnocchi, and is missing when any of the model's features is missing.

    .. note::

        Because the non-mutated sites are a sample, the model is case-control, and
        the ratio of mutation rates between a window and the average window is
        closer to the odds ratio than to this ratio of probabilities. This form
        multiplies the log-scale regional effect by about ``1 - ave``, which pulls
        it toward 1, and caps the rate at ``1 / ave``. That is negligible when DNMs are a
        few percent of the training rows, but not when they are a large share, as
        for CpG contexts.

    :param features_expr: Struct with the model's feature fields, e.g. a feature
        Table's row.
    :param model: Model dictionary (:func:`fit_regional_mutation_rate_model`) or a
        row of :func:`regional_mutation_rate_models_to_ht`.
    :return: Regional mutation rate expression.
    """
    eta = hl.float64(model["intercept"])
    for feature, mean, sd, coef in zip(
        model["features"], model["mean"], model["sd"], model["coefs"]
    ):
        eta = eta + coef * ((features_expr[feature] - mean) / sd)

    return hl.expit(eta) / model["ave"]


def compute_regional_mutation_rate_ht(
    ht: hl.Table,
    models: List[Union[Dict[str, Any], hl.Struct]],
    contexts: List[str],
) -> hl.Table:
    """
    Compute the regional mutation rate of every window for each context.

    Contexts without a model, and windows missing one of a context's features, get a
    rate of 1.0. Models of contexts not in ``contexts`` are ignored. Windows are kept
    only if at least one modeled context in ``contexts`` has a rate.

    :param ht: Feature Table with one row per window, e.g. keyed by an interval.
    :param models: Model dictionaries or rows of
        :func:`regional_mutation_rate_models_to_ht`, each with a ``context``.
    :param contexts: Contexts to include in each window's rates.
    :return: ``ht`` keyed as before, with only an ``adj_r`` dict field mapping each
        context to its rate.
    """
    rate_exprs = {
        m["context"]: regional_mutation_rate_expr(ht.row, m)
        for m in models
        if m["context"] in contexts
    }
    if not rate_exprs:
        logger.warning("No context has a model, so the rate Table is empty.")
        return ht.select(adj_r=hl.empty_dict(hl.tstr, hl.tfloat64)).filter(False)

    ht = ht.select(_rate=hl.struct(**rate_exprs))
    ht = ht.filter(hl.any([hl.is_defined(ht._rate[c]) for c in rate_exprs]))

    # TODO: Windows missing a context's features get 1.0 here, the model's rate for
    #  the average window. Normalization then divides it by the context's factor, so
    #  these windows end up at 1 / factor instead of 1, unlike sites without a
    #  window. Consider leaving them missing here and filling 1.0 after
    #  normalization.
    return ht.select(
        adj_r=hl.dict(
            [
                (
                    c,
                    (
                        hl.or_else(ht._rate[c], 1.0)
                        if c in rate_exprs
                        else hl.float64(1.0)
                    ),
                )
                for c in contexts
            ]
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
    normalization_factors: Dict[str, float],
) -> hl.expr.DictExpression:
    """
    Divide each context's regional mutation rate by its normalization factor.

    Contexts without a factor are left unchanged.

    :param rate_expr: Dict of context to rate, e.g. the ``adj_r`` field of
        :func:`compute_regional_mutation_rate_ht`.
    :param normalization_factors: Dict of context to factor, e.g. evaluated from
        :func:`regional_mutation_rate_normalization_factors_agg_expr`.
    :return: Dict of context to normalized rate.
    """
    factors_expr = hl.literal(normalization_factors, hl.tdict(hl.tstr, hl.tfloat64))

    return hl.dict(
        rate_expr.items().map(lambda kv: (kv[0], kv[1] / factors_expr.get(kv[0], 1.0)))
    )
