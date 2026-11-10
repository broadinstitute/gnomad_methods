"""Functions to build and apply the plateau and coverage models."""

import logging
from typing import Dict, Optional, Tuple, Union

import hail as hl

from gnomad.utils.constraint.constants import COVERAGE_CUTOFF
from gnomad.utils.constraint.context import annotate_mutation_type
from gnomad.utils.constraint.counts import _build_sum_agg_struct

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_utils")
logger.setLevel(logging.INFO)


def calibration_model_group_expr(
    exomes_coverage_expr: hl.expr.Int32Expression,
    cpg_expr: hl.expr.BooleanExpression,
    low_cov_cutoff: Optional[int] = None,
    high_cov_cutoff: int = COVERAGE_CUTOFF,
    upper_cov_cutoff: Optional[int] = None,
    skip_coverage_model: bool = False,
    additional_grouping_exprs: Optional[Dict[str, hl.expr.StringExpression]] = None,
    cpg_in_high_only: bool = False,
) -> hl.expr.StructExpression:
    """
    Get the calibration model grouping annotation for a variant.

    The calibration model expression is a struct with the following fields:

        - high_or_low_coverage: Whether the variant belongs to the high or low coverage
          calibration model. The variant is assigned to the high coverage model if the
          exome coverage is greater than or equal to 'high_cov_cutoff' and less than or
          equal to 'upper_cov_cutoff' (if provided). The variant is assigned to the low
          coverage model if `skip_coverage_model` is False and the exome coverage is
          greater than 'low_cov_cutoff' (if provided) and less than 'high_cov_cutoff'.
        - model_group: Struct of the fields the models are grouped by:

            - cpg: Whether the variant is a CpG (`cpg_expr`).
            - Any fields in `additional_grouping_exprs` (e.g. ``genomic_region``).

    The expression is missing for variants that belong to neither model.

    :param exomes_coverage_expr: Exome coverage expression.
    :param cpg_expr: CpG expression.
    :param low_cov_cutoff: Low coverage cutoff. Default is None.
    :param high_cov_cutoff: High coverage cutoff. Default is COVERAGE_CUTOFF.
    :param upper_cov_cutoff: Upper coverage cutoff. Default is None.
    :param skip_coverage_model: Whether to skip the coverage model. Default is False.
    :param additional_grouping_exprs: Optional Dictionary of additional expressions to
        group by. Default is None.
    :param cpg_in_high_only: Whether to only include `cpg_expr` in the model grouping
        for the high coverage model. Default is False.
    :return: StructExpression with ``high_or_low_coverage`` and ``model_group``
        fields.
    """
    high_cov_expr = exomes_coverage_expr >= high_cov_cutoff
    if upper_cov_cutoff is not None:
        high_cov_expr &= exomes_coverage_expr <= upper_cov_cutoff

    # Sites above `upper_cov_cutoff` fail the high coverage check, so the upper bound
    # on low coverage keeps them out of both models.
    low_cov_expr = hl.bool(False) if skip_coverage_model else hl.bool(True)
    low_cov_expr &= exomes_coverage_expr < high_cov_cutoff
    if low_cov_cutoff is not None:
        low_cov_expr &= exomes_coverage_expr > low_cov_cutoff

    # Define whether the variant should be included in the high or low coverage model.
    model_expr = (
        hl.case().when(high_cov_expr, "high").when(low_cov_expr, "low").or_missing()
    )

    cpg_expr = (
        hl.or_missing(model_expr == "high", cpg_expr) if cpg_in_high_only else cpg_expr
    )
    return hl.or_missing(
        hl.is_defined(model_expr),
        hl.struct(
            high_or_low_coverage=model_expr,
            model_group=hl.struct(
                cpg=cpg_expr,
                **(additional_grouping_exprs or {}),
            ),
        ),
    )


def build_models(
    ht: hl.Table,
    coverage_expr: hl.expr.Int32Expression,
    weighted: bool = False,
    keys: Tuple[str, ...] = (
        "context",
        "ref",
        "alt",
        "methylation_level",
    ),
    model_group_expr: Optional[hl.expr.StructExpression] = None,
    high_cov_definition: int = COVERAGE_CUTOFF,
    upper_cov_cutoff: Optional[int] = None,
    skip_coverage_model: bool = False,
    log10_coverage: bool = True,
    additional_grouping: Tuple[str, ...] = (),
) -> Tuple[Optional[hl.expr.ArrayNumericExpression], hl.expr.DictExpression]:
    """
    Build coverage and plateau models.

    This function builds models (plateau_models) using linear regression to calibrate
    mutation rate estimates against the proportion observed of each substitution,
    context, and methylation level in `ht`.

    Two plateau models are fit, one for CpG transitions, and one for the remainder of
    sites (transversions and non CpG transitions).

    The plateau models only consider high coverage sites, or sites above a median
    coverage of `high_cov_definition` and median coverage below `upper_cov_cutoff`.

    Plateau model: adjusts proportion of expected variation based on location in the
    genome and CpG status.

    The x and y of the plateau models:
        - x: `mu_snp` - mutation rate
        - y: proportion observed ('observed_variants' / 'possible_variants')

    This function also builds models (coverage models) to calibrate the proportion of
    expected variation at low coverage sites (sites below `high_cov_definition`).

    The coverage models are built by creating a scaling factor across all high coverage
    sites, applying this ratio to the low coverage sites, and running a linear
    regression.

    Coverage model: corrects proportion of expected variation at low coverage sites.
    Low coverage sites are defined as sites with median coverage < `high_cov_definition`.

    The x and y of the coverage model:

        - x: groupings of exome coverage at low coverage sites (log10 transformed if
          requested)
        - y: sum('observed_variants') / (``high_coverage_scale_factor`` *
          sum('possible_variants' * 'mu_snp')) at low coverage sites

    ``high_coverage_scale_factor`` = sum('observed_variants') /
    sum('possible_variants' * 'mu_snp') at high coverage sites

    .. note::

        This function expects that the input Table (`ht`) contains only high quality
        synonymous variants below 0.1% frequency.

    The following fields are expected in `ht`:

        - context: trinucleotide genomic context.
        - ref: the reference allele.
        - alt: the alternate allele.
        - methylation_level: methylation level.
        - mu_snp: mutation rate.
        - cpg: whether the variant is a CpG.
        - observed_variants: the number of observed variants in the dataset for each
          variant. Note that the term "variant" here refers to a specific substitution,
          context, methylation level, and coverage combination.
        - possible_variants: the number of possible variants in the dataset for each
          variant.

    :param ht: Input Table.
    :param coverage_expr: Expression that defines the coverage metric.
    :param weighted: Whether to weight the plateau models (a linear regression
        model) by 'possible_variants'. Default is False.
    :param keys: Annotations used to group observed and possible variant counts.
        Default is ("context", "ref", "alt", "methylation_level").
    :param model_group_expr: Expression with ``high_or_low_coverage`` and
        ``model_group`` annotations to group variants into high or low coverage models.
        If ``model_group`` includes a ``genomic_region`` field, the coverage model is
        built only on sites where it is "autosome_or_par"; otherwise all sites are
        used. If not provided, the ``calibration_model_group_expr`` function is used to
        define the grouping, which does not include ``genomic_region``.
    :param high_cov_definition: Coverage threshold for high/low classification.
        Sites with coverage above this cutoff are considered well covered. Default is
        ``COVERAGE_CUTOFF``.
    :param upper_cov_cutoff: Upper coverage cutoff. Sites with coverage above this
        cutoff are excluded from the high coverage Table. Default is None.
    :param skip_coverage_model: Whether to skip generating the coverage model. If set
        to True, None is returned instead of the coverage model. Default is False.
    :param log10_coverage: Whether to convert coverage sites with log10 when building
        the coverage model. Default is True.
    :param additional_grouping: Additional annotations to group by before
        counting the observed and possible variants. Default is ().
    :return: Tuple of (coverage model, plateau models). The coverage model is an
        ``[intercept, slope]`` ArrayNumericExpression, or None when
        ``skip_coverage_model`` is True.
    """
    if model_group_expr is None:
        # Define whether the variant should be included in the high or low coverage
        # model.
        model_group_expr = calibration_model_group_expr(
            coverage_expr,
            ht.cpg,
            low_cov_cutoff=0,
            high_cov_cutoff=high_cov_definition,
            upper_cov_cutoff=upper_cov_cutoff,
            skip_coverage_model=skip_coverage_model,
            cpg_in_high_only=True,
        )

    grouping = keys + additional_grouping
    mu_type_fields = ("cpg", "transition", "mutation_type", "mutation_type_model")
    # all() accepts a generator expression directly (no list needed) and
    # short-circuits on the first False.
    has_mu_type = all(x in ht.row for x in mu_type_fields)
    grouping += mu_type_fields if has_mu_type else ()

    grouping_exprs = {"build_model": model_group_expr}
    if not skip_coverage_model:
        grouping_exprs["exomes_coverage"] = hl.or_missing(
            model_group_expr.high_or_low_coverage == "low", coverage_expr
        )

    # Collapse per-site rows into one row per (context, ref, alt,
    # methylation_level, …) group with summed observed/possible counts.
    # model_group_expr and (optionally) exomes_coverage are carried along so
    # that plateau and coverage models can be built on the grouped table.
    # mu_snp is identical for all rows within a group so grabbing one value suffices.
    ht = (
        ht.group_by(*grouping, **grouping_exprs)
        .aggregate(
            mu_snp=hl.agg.take(ht.mu_snp, 1)[0],
            **_build_sum_agg_struct(
                fields_to_sum=["observed_variants", "possible_variants"], t=ht
            ),
        )
        .key_by(*keys)
    )

    if not has_mu_type:
        ht = annotate_mutation_type(ht)

    # Build plateau models.
    is_high_expr = ht.build_model.high_or_low_coverage == "high"
    agg_expr = {
        "plateau": hl.agg.filter(
            is_high_expr,
            build_plateau_models(
                ht.mu_snp,
                ht.observed_variants,
                ht.possible_variants,
                model_group_expr=ht.build_model.model_group,
                weighted=weighted,
            ),
        )
    }

    if not skip_coverage_model:
        # The coverage model is only built using the full dataset observed variants
        # so use the first element if the observed_variants is an array.
        obs_is_array = isinstance(ht.observed_variants, hl.expr.ArrayExpression)
        obs_expr = ht.observed_variants[0] if obs_is_array else ht.observed_variants

        # Restrict the coverage model to autosome/PAR sites when the model grouping
        # includes a genomic region, otherwise use all sites.
        model_group = ht.build_model.model_group
        autosome_or_par_expr = (
            model_group.genomic_region == "autosome_or_par"
            if "genomic_region" in model_group
            else hl.bool(True)
        )

        # Create a metric that represents the relative mutability of the exome calculated
        # on high coverage sites and will be used as scaling factor when building the
        # coverage model.
        agg_expr["high_coverage_scale_factor"] = hl.agg.filter(
            is_high_expr & autosome_or_par_expr,
            hl.agg.sum(obs_expr) / hl.agg.sum(ht.possible_variants * ht.mu_snp),
        )

        # Get the observed variant count and mu_snp for low coverage sites.
        agg_expr["coverage"] = hl.agg.filter(
            (ht.build_model.high_or_low_coverage == "low") & autosome_or_par_expr,
            hl.agg.group_by(
                ht.exomes_coverage,
                hl.struct(
                    obs=hl.agg.sum(obs_expr),
                    mu_snp=hl.agg.sum(ht.possible_variants * ht.mu_snp),
                ),
            ),
        )

    models = ht.aggregate(hl.struct(**agg_expr), _localize=False)

    # Build coverage model.
    coverage_model = None
    if not skip_coverage_model:
        coverage_model = models.coverage.map_values(
            lambda x: x.annotate(
                low_coverage_oe=x.obs / (models.high_coverage_scale_factor * x.mu_snp)
            )
        )

        # TODO: consider weighting here as well.
        coverage_model = (
            coverage_model.items()
            .aggregate(
                lambda x: build_coverage_model(
                    x[1].low_coverage_oe, x[0], log10_coverage=log10_coverage
                )
            )
            .beta
        )

    return coverage_model, models.plateau


def build_plateau_models(
    mu_snp_expr: hl.expr.Float64Expression,
    observed_variants_expr: Union[hl.expr.ArrayExpression, hl.expr.Int64Expression],
    possible_variants_expr: Union[hl.expr.ArrayExpression, hl.expr.Int64Expression],
    weighted: bool = False,
    cpg_expr: Optional[hl.expr.BooleanExpression] = None,
    model_group_expr: Optional[hl.expr.StructExpression] = None,
) -> Union[hl.expr.DictExpression, hl.expr.ArrayExpression, hl.expr.StructExpression]:
    """
    Build plateau models to calibrate mutation rate against proportion observed.

    Fits a linear regression of ``observed_variants_expr / possible_variants_expr``
    on ``mu_snp_expr``. When either observed or possible expressions are arrays
    (e.g., one model per downsampling), the regression is applied element-wise
    via ``hl.agg.array_agg``.

    When ``model_group_expr`` or ``cpg_expr`` is provided, the result is a
    ``DictExpression`` keyed by the grouping struct. Otherwise the result is the
    regression beta directly.

    :param mu_snp_expr: Mutation rate expression.
    :param observed_variants_expr: Observed variant counts (scalar or array).
    :param possible_variants_expr: Possible variant counts (scalar or array).
    :param weighted: If True, use weighted least squares with
        ``possible_variants_expr`` as weights. Default is False.
    :param cpg_expr: Boolean expression indicating CpG sites. When provided,
        adds a ``cpg`` field to the grouping struct.
    :param model_group_expr: Struct expression to group by in the aggregation.
    :return: Regression betas, optionally grouped by ``model_group_expr``
        (and/or ``cpg_expr``).
    """
    obs_is_array = isinstance(observed_variants_expr, hl.expr.ArrayExpression)
    pos_is_array = isinstance(possible_variants_expr, hl.expr.ArrayExpression)

    def _linreg(
        o: hl.expr.NumericExpression,
        p: hl.expr.NumericExpression,
    ) -> hl.expr.StructExpression:
        """
        Run linear regression of observed/possible on mutation rate.

        :param o: Observed variant count expression.
        :param p: Possible variant count expression.
        :return: Regression beta coefficients.
        """
        return hl.agg.linreg(
            o / p, [1, mu_snp_expr], weight=p if weighted else None
        ).beta

    if obs_is_array and pos_is_array:
        agg_expr = hl.agg.array_agg(
            lambda x: _linreg(*x),
            hl.zip(observed_variants_expr, possible_variants_expr),
        )
    elif obs_is_array:
        agg_expr = hl.agg.array_agg(
            lambda x: _linreg(x, possible_variants_expr), observed_variants_expr
        )
    elif pos_is_array:
        agg_expr = hl.agg.array_agg(
            lambda x: _linreg(observed_variants_expr, x), possible_variants_expr
        )
    else:
        agg_expr = _linreg(observed_variants_expr, possible_variants_expr)

    if model_group_expr is None and cpg_expr is None:
        return agg_expr

    model_group_expr = model_group_expr or hl.struct()
    if cpg_expr is not None:
        model_group_expr = model_group_expr.annotate(cpg=cpg_expr)

    return hl.agg.group_by(model_group_expr, agg_expr)


def build_coverage_model(
    low_coverage_oe_expr: hl.expr.Float64Expression,
    coverage_expr: hl.expr.Float64Expression,
    log10_coverage: bool = False,
) -> hl.expr.StructExpression:
    """
    Build coverage model.

    This function uses linear regression to build a model of coverage to correct
    proportion of expected variation at low coverage sites.

    The x and y of the coverage model:

        - x: `coverage_expr`
        - y: `low_coverage_oe_expr`

    :param low_coverage_oe_expr: The Float64Expression of observed:expected ratio
        for a given coverage level.
    :param coverage_expr: The Float64Expression of the coverage expression.
    :param log10_coverage: Whether to convert coverage sites by log10 when building the
        coverage model. Default is False.
    :return: StructExpression with intercept and slope of the model.
    """
    if log10_coverage:
        logger.info(
            "Converting coverage sites by log10 when building the coverage model."
        )
        coverage_expr = hl.log10(coverage_expr)

    return hl.agg.linreg(low_coverage_oe_expr, [1, coverage_expr])


def apply_plateau_models(
    mu_expr: hl.Float64Expression,
    plateau_models_expr: hl.ArrayExpression,
) -> Union[hl.ArrayExpression, hl.Float64Expression]:
    """
    Compute the predicted probability observed.

    The predicted probability observed is computed as the mutation rate adjusted by the
    plateau model.

    :param mu_expr: Mutation rate expression.
    :param plateau_models_expr: This can be either a single plateau model, where the
        first element is the intercept and the second element is the slope, or an array
        of plateau models.
    :return: Float64Expression of the predicted probability observed, or an
        ArrayExpression of them when ``plateau_models_expr`` is an array of models.
    """

    def _apply_model(plateau_model: hl.ArrayExpression) -> hl.Float64Expression:
        """
        Apply the plateau model to the mutation rate expression.

        :param plateau_model: ArrayExpression of the plateau model.
        :return: Float64Expression of the predicted probability observed.
        """
        slope = plateau_model[1]
        intercept = plateau_model[0]
        ppo_expr = mu_expr * slope + intercept

        return ppo_expr

    if plateau_models_expr.dtype.element_type == hl.tarray(hl.tfloat64):
        return plateau_models_expr.map(lambda x: _apply_model(x))

    return _apply_model(plateau_models_expr)


def coverage_correction_expr(
    coverage_expr: hl.Float64Expression,
    coverage_model: Union[Tuple[float, float], hl.expr.ArrayNumericExpression],
    low_coverage_expr: Optional[hl.BooleanExpression] = None,
    coverage_cutoff: Optional[int] = None,
    log10_coverage: bool = False,
) -> hl.Float64Expression:
    """
    Compute the coverage correction expression.

    .. note::

        One and only one of `low_coverage_expr` or `coverage_cutoff` must be specified.

    The coverage correction expression is computed as follows:

        - If the coverage is 0, the coverage correction is 0.
        - If the low coverage expression (`low_coverage_expr`) is True, or the
          coverage (`coverage_expr`) is below the coverage cutoff (`coverage_cutoff`),
          the coverage model is applied to the coverage.
        - Otherwise, the coverage correction is 1.

    :param coverage_expr: Float64Expression of the coverage.
    :param coverage_model: Coverage model as ``[intercept, slope]``, either a tuple of
        floats or an ArrayNumericExpression (as returned by ``build_models``).
    :param low_coverage_expr: Optional BooleanExpression indicating whether the site is
        a low coverage site, and the coverage model should be applied. Default is None.
    :param coverage_cutoff: Optional coverage cutoff. If specified, the coverage model
        is applied to sites with coverage below this cutoff. Default is None.
    :param log10_coverage: Whether to convert coverage sites by log10 when applying the
        coverage model. Default is False.
    :return: Float64Expression of the coverage correction.
    """
    if low_coverage_expr is None and coverage_cutoff is None:
        raise ValueError(
            "Either 'low_coverage_expr' or 'coverage_cutoff' must be specified!"
        )
    if low_coverage_expr is not None and coverage_cutoff is not None:
        raise ValueError(
            "Only one of 'low_coverage_expr' or 'coverage_cutoff' can be specified!"
        )

    if coverage_cutoff is not None:
        low_coverage_expr = coverage_expr < coverage_cutoff

    if log10_coverage:
        cov_corr_expr = hl.log10(coverage_expr)
    else:
        cov_corr_expr = coverage_expr

    return (
        hl.case()
        .when(coverage_expr == 0, 0)
        .when(
            low_coverage_expr,
            coverage_model[1] * cov_corr_expr + coverage_model[0],
        )
        .default(1)
    )


def apply_models(
    mu_expr: hl.expr.Float64Expression,
    plateau_models_expr: hl.expr.ArrayExpression,
    possible_variants_expr: hl.expr.Int64Expression,
    coverage_model: Optional[
        Union[Tuple[float, float], hl.expr.ArrayNumericExpression]
    ] = None,
    coverage_expr: Optional[hl.expr.Int32Expression] = None,
    cpg_expr: Optional[hl.expr.BooleanExpression] = None,
    model_group_expr: Optional[hl.expr.StructExpression] = None,
    high_cov_definition: int = COVERAGE_CUTOFF,
    log10_coverage: bool = True,
) -> hl.expr.StructExpression:
    """
    Apply calibration models to compute expected variant counts.

    Applies plateau and (optionally) coverage models to produce a struct with
    ``mu``, ``predicted_proportion_observed``, ``expected_variants``, and
    (when a coverage model is provided) ``coverage_correction``.

    .. note::

        When a coverage model is applied, a row whose ``model_group_expr`` is
        missing gets missing ``mu``, ``expected_variants``, and
        ``coverage_correction`` (unless its coverage is 0, which always gives a
        correction of 0). ``calibration_model_group_expr`` is missing for sites in
        neither model, such as sites above ``upper_cov_cutoff`` or at or below
        ``low_cov_cutoff``. ``hl.agg.sum`` skips missing values, so filter these rows
        out before applying the models if they should not silently drop out of
        aggregated expected counts.

    :param mu_expr: Mutation rate expression.
    :param plateau_models_expr: Single plateau model (array of [intercept,
        slope]) or an array of plateau models.
    :param possible_variants_expr: Possible variant counts to multiply the
        predicted proportion observed by.
    :param coverage_model: Coverage model as ``[intercept, slope]``, either a tuple of
        floats or an ArrayNumericExpression (as returned by ``build_models``).
        Default is None.
    :param coverage_expr: Int32Expression of the coverage. Required when
        ``coverage_model`` is provided or ``model_group_expr`` is None.
    :param cpg_expr: BooleanExpression indicating CpG sites. Required when
        ``model_group_expr`` is None.
    :param model_group_expr: Expression with ``high_or_low_coverage``
        annotation to group variants into high or low coverage models. If not
        provided, ``calibration_model_group_expr`` is used. Default is None.
    :param high_cov_definition: Coverage threshold for high/low classification.
        Default is ``COVERAGE_CUTOFF``.
    :param log10_coverage: Whether to log10-transform coverage when applying
        the coverage model. Default is True.
    :return: StructExpression with "mu", "predicted_proportion_observed",
        "expected_variants", and optionally "coverage_correction".
    """
    if coverage_model is not None and coverage_expr is None:
        raise ValueError(
            "If 'coverage_model' is specified, 'coverage_expr' must also be specified!"
        )
    if model_group_expr is None and (coverage_expr is None or cpg_expr is None):
        raise ValueError(
            "If 'model_group_expr' is not specified, 'coverage_expr' and 'cpg_expr' must"
            " be specified!"
        )

    if model_group_expr is None:
        # Get the annotations relevant for applying the calibration models.
        model_group_expr = calibration_model_group_expr(
            coverage_expr,
            cpg_expr,
            high_cov_cutoff=high_cov_definition,
            skip_coverage_model=coverage_model is None,
        )

    # Apply plateau models.
    ppo_expr = apply_plateau_models(mu_expr, plateau_models_expr)
    apply_expr = hl.struct(
        mu=mu_expr * possible_variants_expr,
        predicted_proportion_observed=ppo_expr,
        expected_variants=ppo_expr * possible_variants_expr,
    )

    # Get the coverage correction expression if a coverage model is provided.
    if coverage_model is not None:
        cov_corr_expr = coverage_correction_expr(
            coverage_expr,
            coverage_model,
            low_coverage_expr=model_group_expr.high_or_low_coverage == "low",
            log10_coverage=log10_coverage,
        )
        apply_expr = apply_expr.annotate(
            mu=apply_expr.mu * cov_corr_expr,
            expected_variants=apply_expr.expected_variants * cov_corr_expr,
            coverage_correction=cov_corr_expr,
        )

    return apply_expr
