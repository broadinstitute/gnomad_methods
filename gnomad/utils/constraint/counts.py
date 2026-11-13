"""Functions to count observed and possible variants and sum them by group."""

import logging
from typing import Dict, List, Optional, Tuple, Union

import hail as hl

from gnomad.utils.constraint.constants import DEFAULT_FIELDS_TO_SUM

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("constraint_utils")
logger.setLevel(logging.INFO)


def _resolve_row_annotation_expr(
    t: Optional[Union[hl.Table, hl.MatrixTable]] = None,
    annotation_name: Optional[str] = None,
    expr: Optional[hl.expr.Expression] = None,
    expr_param_name: Optional[str] = None,
) -> hl.expr.Expression:
    """
    Get a row annotation from a Table/MatrixTable, or return a provided expression.

    Provides a consistent pattern for functions that accept either an explicit
    Hail expression or fall back to a well-known field on a Table. This avoids
    duplicating "resolve expr or look it up on ht" logic across callers.

    If ``expr`` is provided it is returned directly. Otherwise, ``t`` and
    ``annotation_name`` must both be supplied and the named row field is looked up
    on ``t``.

    Example usage inside a public function::

        def my_function(
            ht: hl.Table,
            freq_expr: Optional[hl.expr.StructExpression] = None,
        ) -> ...:
            freq_expr = _resolve_row_annotation_expr(
                t=ht,
                annotation_name="freq",
                expr=freq_expr,
                expr_param_name="freq_expr",
            )
            # freq_expr is now guaranteed to be a valid expression —
            # either the caller's explicit value or ht.freq.

    :param t: Input Table or MatrixTable to look up ``annotation_name`` on. Required
        when ``expr`` is None.
    :param annotation_name: Name of the field to retrieve from ``t``. Required when
        ``expr`` is None.
    :param expr: Expression to return directly. When provided, ``t`` and
        ``annotation_name`` are ignored.
    :param expr_param_name: Human-readable name for ``expr``, used in log/error
        messages. Defaults to "expr".
    :return: The resolved Hail expression.
    """
    if expr is None and (t is None or annotation_name is None):
        raise ValueError("Either 't' and 'annotation_name' or 'expr' must be provided.")

    expr_param_name = expr_param_name or "expr"
    if expr is None and annotation_name in t.row:
        logger.info(
            "%s was not provided, using '%s'.", expr_param_name, annotation_name
        )
        expr = t[annotation_name]
    elif expr is None:
        raise ValueError(
            f"{expr_param_name} was not provided and '{annotation_name}' is "
            f"not present in the input Table or MatrixTable."
        )

    return expr


def variant_observed_expr(
    freq_expr: Optional[
        Union[hl.expr.StructExpression, hl.expr.ArrayExpression]
    ] = None,
    ht: Optional[hl.Table] = None,
    singleton: bool = False,
    max_af: Optional[float] = None,
    count_missing: bool = False,
) -> hl.expr.Int32Expression:
    """
    Return 1 if a variant meets frequency criteria, 0 otherwise.

    One of ``ht`` or ``freq_expr`` must be provided. If ``freq_expr`` is not provided,
    ``ht.freq`` is used. When ``freq_expr`` is an ArrayExpression, the first
    element is used.

    The returned count is 1 when:

        - ``singleton`` is True and ``freq_expr.AC == 1``.
        - ``max_af`` is not None and ``freq_expr.AC > 0`` and
          ``freq_expr.AF <= max_af``.
        - Neither ``singleton`` nor ``max_af`` is set but ``freq_expr`` is
          available: ``freq_expr.AC > 0``.
        - Neither ``singleton`` nor ``max_af`` is set and no ``freq_expr`` is
          available: unconditionally 1.

    :param ht: Input Table. Used to look up ``freq`` when ``freq_expr`` is None.
    :param freq_expr: StructExpression (or ArrayExpression of structs) with
        ``AC`` and ``AF`` fields. When an array, ``freq_expr[0]`` is used.
    :param singleton: Count only singletons (AC == 1). When True, ``max_af`` is
        ignored. Default is False.
    :param max_af: Maximum allele frequency threshold. Default is None (no
        cutoff).
    :param count_missing: Whether to count a variant whose count expression is
        missing (e.g. frequency is None). True counts it as 1, False as 0. Default
        is False.
    :return: Int32Expression equal to 0 or 1.
    """
    if ht is None and freq_expr is None:
        raise ValueError("Either ht or freq_expr must be provided.")

    if freq_expr is None and (max_af is not None or singleton or "freq" in ht.row):
        freq_expr = _resolve_row_annotation_expr(ht, "freq", freq_expr, "freq_expr")
    if isinstance(freq_expr, hl.expr.ArrayExpression):
        freq_expr = freq_expr[0]

    if singleton:
        count_var = freq_expr.AC == 1
    elif max_af is not None:
        count_var = (freq_expr.AC > 0) & (freq_expr.AF <= max_af)
    elif freq_expr is not None:
        count_var = freq_expr.AC > 0
    else:
        count_var = True

    return hl.int(hl.or_else(count_var, count_missing))


def variant_observed_and_possible_expr(
    freq_expr: hl.ArrayExpression,
    max_af: Optional[float] = None,
    use_possible_adj: bool = True,
) -> hl.expr.StructExpression:
    """
    Return per-variant observed and possible count expressions.

    For each element of the frequency array, ``observed_variants`` is 1 when the
    variant meets the frequency criteria (AC > 0, and AF <= ``max_af`` if set) and 0
    otherwise. ``possible_variants`` uses the same criteria, except that a missing
    frequency counts as possible: a variant absent from the dataset could still have
    been observed. A variant with a defined frequency that fails the criteria
    (AC == 0, or AF > ``max_af``) is therefore 0 in both, so it is excluded from the
    denominator as well as the numerator. With ``max_af`` set, this removes common
    variants from the possible counts.

    When ``freq_expr`` itself is missing (e.g. the variant is not in the dataset),
    ``observed_variants`` is missing and ``possible_variants`` is 1, or missing when
    ``use_possible_adj`` is False.

    When ``use_possible_adj`` is True (default), the possible count is a scalar
    derived from the first (adj) element of the frequency array. When False, it
    is an array with one entry per array element, matching the shape of
    ``observed_variants``.

    :param freq_expr: Array of frequency structs with ``AC`` and ``AF`` fields.
    :param max_af: Maximum allele frequency threshold. Variants above it are excluded
        from both counts. Default is None (no cutoff).
    :param use_possible_adj: If True, compute possible count from only the adj
        (first) frequency element. If False, compute per-downsampling. Default
        is True.
    :return: Struct with ``observed_variants`` (array of int) and
        ``possible_variants`` (int if ``use_possible_adj``, else array of int).
    """
    pos_expr = hl.array([freq_expr[0]]) if use_possible_adj else freq_expr
    pos_expr = pos_expr.map(
        lambda x: variant_observed_expr(freq_expr=x, max_af=max_af, count_missing=True)
    )
    pos_expr = pos_expr[0] if use_possible_adj else pos_expr
    return hl.struct(
        observed_variants=freq_expr.map(
            lambda x: variant_observed_expr(freq_expr=x, max_af=max_af)
        ),
        possible_variants=pos_expr,
    )


def weighted_sum_agg_expr(
    expr: Union[hl.expr.ArrayNumericExpression, hl.expr.NumericExpression],
    weight_expr: Union[hl.expr.ArrayNumericExpression, hl.expr.NumericExpression],
) -> Union[hl.expr.ArrayExpression, hl.expr.NumericExpression]:
    """
    Return the weighted aggregate sum of ``expr`` weighted by ``weight_expr``.

    Both parameters may be scalar or array numeric expressions:

        - Both arrays: elements are multiplied pairwise and summed per-element
          across rows (``hl.agg.array_sum``).
        - Both scalars: standard ``hl.agg.sum(expr * weight_expr)``.
        - Mixed: the scalar is broadcast across the array elements and summed
          per-element across rows.

    :param expr: Numeric expression (scalar or array) to be weighted.
    :param weight_expr: Numeric expression (scalar or array) to weight by.
    :return: Weighted aggregate sum expression.
    """
    expr_is_array = isinstance(expr, hl.expr.ArrayNumericExpression)
    weight_is_array = isinstance(weight_expr, hl.expr.ArrayNumericExpression)
    if expr_is_array and weight_is_array:
        return hl.agg.array_sum(hl.zip(expr, weight_expr).map(lambda x: x[0] * x[1]))
    elif not expr_is_array and not weight_is_array:
        return hl.agg.sum(expr * weight_expr)
    else:
        return hl.agg.array_sum(expr * weight_expr)


def counts_agg_expr(
    freq_expr: Optional[
        Union[hl.expr.ArrayExpression, hl.expr.StructExpression]
    ] = None,
    ht: Optional[hl.Table] = None,
    count_singletons: bool = False,
    max_af: Optional[float] = None,
    count_missing: bool = False,
) -> hl.expr.StructExpression:
    """
    Return an aggregation expression for variant and singleton counts.

    Aggregates per-variant counts (via ``variant_observed_expr``) across
    rows. Each variant contributes 0 or 1 to the count based on its frequency
    metadata (AC, AF). The result is a struct with ``variant_count`` and, when
    ``count_singletons`` is True, ``singleton_count``.

    One of ``freq_expr`` or ``ht`` must be provided. If ``freq_expr`` is not
    provided, ``ht.freq`` is used as the fallback.

    The shape of ``freq_expr`` controls whether the output counts are scalars or
    arrays:

        - **StructExpression** (single frequency entry with ``AC`` / ``AF``):
          returns scalar ``variant_count`` and ``singleton_count``.
        - **ArrayExpression** (array of frequency structs, e.g. one per
          downsampling): returns array-valued counts where each element corresponds
          to a position in the input array. Internally the array is mapped through
          ``variant_observed_expr`` and summed element-wise with
          ``hl.agg.array_sum``.

    When ``max_af`` is set, only variants with ``AF <= max_af`` and ``AC > 0``
    are counted. Singleton counting (``count_singletons=True``) counts only
    variants with ``AC == 1``, independent of the ``max_af`` filter.

    :param freq_expr: Frequency expression — an ArrayExpression of structs or a
        single StructExpression with ``AC`` and ``AF`` fields. If None, falls
        back to ``ht.freq``.
    :param ht: Input Table. Used to look up ``freq`` when ``freq_expr`` is None.
    :param count_singletons: If True, include a ``singleton_count`` field in
        the result. Default is False.
    :param max_af: Maximum allele frequency threshold. Variants with
        ``AF > max_af`` are excluded from ``variant_count``. Does not affect
        ``singleton_count``. Default is None (no cutoff).
    :param count_missing: Whether to count a variant whose frequency is missing.
        True counts it as 1, False as 0. Default is False.
    :return: Aggregation StructExpression with ``variant_count`` (and
        optionally ``singleton_count``). Values are scalars when ``freq_expr`` is
        a StructExpression, or arrays when it is an ArrayExpression.
    """
    if ht is None and freq_expr is None:
        raise ValueError("Either ht or freq_expr must be provided.")

    freq_expr = _resolve_row_annotation_expr(ht, "freq", freq_expr, "freq_expr")

    is_struct = isinstance(freq_expr, hl.expr.StructExpression)
    if is_struct:
        freq_expr = hl.array([freq_expr])

    params = {"variant_count": {"singleton": False, "max_af": max_af}}
    if count_singletons:
        params["singleton_count"] = {"singleton": True}

    count_expr = {
        k: freq_expr.map(
            lambda x: variant_observed_expr(
                freq_expr=x, **p, count_missing=count_missing
            )
        )
        for k, p in params.items()
    }

    agg_expr = {k: hl.agg.array_sum(v) for k, v in count_expr.items()}
    agg_expr = {k: agg_expr[k][0] for k in agg_expr} if is_struct else agg_expr

    return hl.struct(**agg_expr)


def count_observed_and_possible_by_group(
    ht: hl.Table,
    possible_expr: hl.expr.Int32Expression,
    observed_expr: hl.expr.ArrayExpression,
    additional_grouping: Union[List[str], Tuple[str, ...]] = ("methylation_level",),
    partition_hint: int = 100,
    weight_exprs: Optional[
        Union[
            List[str],
            Dict[str, Union[hl.expr.ArrayExpression, hl.expr.NumericExpression]],
        ]
    ] = None,
    additional_agg_sum_exprs: Optional[
        Union[
            List[str],
            Dict[str, Union[hl.expr.ArrayExpression, hl.expr.NumericExpression]],
        ]
    ] = None,
) -> hl.Table:
    """
    Aggregate observed and possible variant counts by substitution context group.

    Groups rows by ``context``, ``ref``, ``alt``, and any fields in
    ``additional_grouping``, then sums observed and possible counts within each
    group.

    :param ht: Input Table with ``context``, ``ref``, ``alt`` fields and any
        fields named in ``additional_grouping``.
    :param possible_expr: Per-variant possible count (scalar).
    :param observed_expr: Per-variant observed count (array, one element per
        array element of the frequency array).
    :param additional_grouping: Field names to append to the base
        ``(context, ref, alt)`` grouping. Default is ``("methylation_level",)``.
    :param partition_hint: Target number of partitions for the ``group_by``.
        Default is 100.
    :param weight_exprs: Weighted sums of ``possible_expr`` to include. Pass
        field names (looked up on ``ht``) or a dict mapping output names to
        weight expressions. Each produces
        ``weighted_sum_agg_expr(possible_expr, weight)``.
    :param additional_agg_sum_exprs: Extra fields to sum alongside
        observed/possible. Pass field names (looked up on ``ht``) or a dict
        mapping output names to expressions. Arrays use
        ``hl.agg.array_sum``; scalars use ``hl.agg.sum``.
    :return: Grouped Table with ``observed_variants``, ``possible_variants``,
        and any weighted/additional sum fields.
    """
    if isinstance(weight_exprs, list):
        weight_exprs = {k: ht[k] for k in weight_exprs}
    if isinstance(additional_agg_sum_exprs, list):
        additional_agg_sum_exprs = {k: ht[k] for k in additional_agg_sum_exprs}

    weight_exprs = weight_exprs or {}
    additional_agg_sum_exprs = additional_agg_sum_exprs or {}

    # Build the grouping struct for the variant count aggregation.
    grouping = hl.struct(context=ht.context, ref=ht.ref, alt=ht.alt)
    grouping = grouping.annotate(
        **{g: ht[g] for g in additional_grouping if g not in grouping}
    )
    logger.info(
        "The following annotations will be used to group the input Table rows when"
        " counting variants: %s.",
        ", ".join(grouping.keys()),
    )

    agg_expr = {
        "observed_variants": hl.agg.array_sum(observed_expr),
        "possible_variants": hl.agg.sum(possible_expr),
    }

    # Update the possible variant count aggregation expression to include weighted sums
    # of possible variant counts.
    agg_expr.update(
        {k: weighted_sum_agg_expr(possible_expr, v) for k, v in weight_exprs.items()}
    )

    # Get sum aggregation expressions for requested fields.
    agg_expr.update(_build_sum_agg_struct(exprs_to_sum=additional_agg_sum_exprs))

    # Apply each variant count aggregation in `agg_expr` to get counts for all
    # combinations of `grouping`.
    ht = ht.group_by(**grouping).partition_hint(partition_hint).aggregate(**agg_expr)

    return ht


def _build_sum_agg_struct(
    fields_to_sum: Optional[List[str]] = None,
    exprs_to_sum: Optional[
        Union[hl.expr.StructExpression, Dict[str, hl.expr.NumericExpression]]
    ] = None,
    t: Optional[Union[hl.Table, hl.expr.StructExpression]] = None,
) -> hl.expr.StructExpression:
    """
    Return an aggregation expression to sum fields or expressions in a Table/StructExpression.

    The aggregation expression is a struct with the sum or array_sum of the fields or
    expressions provided in `fields_to_sum` or `exprs_to_sum`.

    :param fields_to_sum: List of fields to sum. Must provide ``t`` if using this.
        Default is None.
    :param exprs_to_sum: Dictionary of expressions to sum. Default is None.
    :param t: Optional Table or StructExpression to get `fields_to_sum` from. Default
        is None.
    :return: Aggregation expression to sum fields or expressions in the Table.
    """
    if fields_to_sum is None and exprs_to_sum is None:
        raise ValueError("Either 'fields_to_sum' or 'exprs_to_sum' must be provided.")
    if fields_to_sum is not None and t is None:
        raise ValueError("t must be provided if 'fields_to_sum' is provided.")

    exprs_to_sum = exprs_to_sum or {}
    exprs_to_sum = hl.struct(**exprs_to_sum, **{f: t[f] for f in fields_to_sum or []})

    return hl.struct(
        **{
            k: (
                hl.agg.array_sum(v)
                if isinstance(v, hl.ArrayExpression)
                else hl.agg.sum(v)
            )
            for k, v in exprs_to_sum.items()
        }
    )


def aggregate_constraint_metrics_expr(
    t: Union[hl.Table, hl.StructExpression],
    fields_to_sum: Union[List[str], Tuple[str, ...]] = DEFAULT_FIELDS_TO_SUM,
    additional_exprs_to_sum: Optional[Dict[str, hl.expr.Expression]] = None,
) -> hl.expr.StructExpression:
    """
    Get an aggregation expression for the sum of specified fields and expressions.

    An aggregate sum or array sum is created for each field in ``fields_to_sum``.

    :param t: Input Table or StructExpression.
    :param fields_to_sum: Fields in ``t`` to get an aggregate sum expression for.
        Default is ``DEFAULT_FIELDS_TO_SUM``.
    :param additional_exprs_to_sum: Dictionary of additional expressions to get
        an aggregate sum expression for. Field names are the keys and expressions are
        the values. Default is None.
    :return: StructExpression with the sum of each field and expression.
    """
    return _build_sum_agg_struct(
        fields_to_sum=list(fields_to_sum),
        exprs_to_sum=additional_exprs_to_sum,
        t=t,
    )
