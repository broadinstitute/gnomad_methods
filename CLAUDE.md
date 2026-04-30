# gnomad_methods Project Reference

## Project Overview

Shared Hail utility library for gnomAD pipelines. Provides reusable functions for constraint analysis, variant QC, sample QC, VEP processing, resource management, and general genomics operations. Used as a dependency by `gnomad_constraint`, `gnomad_qc`, `gnomad_mnv`, and other gnomAD repos.

## Package Structure

| Directory | Purpose |
|-----------|---------|
| `gnomad/utils/constraint.py` | Constraint pipeline utilities (mutation rate, model building, o/e, pLI, z-scores) |
| `gnomad/utils/vep.py` | VEP annotation processing (consequence filtering, LOFTEE, MANE Select) |
| `gnomad/utils/annotations.py` | General variant annotations (age, frequency, quality) |
| `gnomad/utils/file_utils.py` | File I/O helpers, struct printing, array conversion |
| `gnomad/utils/filtering.py` | Frequency and variant filtering utilities |
| `gnomad/utils/gen_stats.py` | Statistical helper functions |
| `gnomad/utils/reference_genome.py` | Reference genome utilities |
| `gnomad/utils/release.py` | Release table formatting |
| `gnomad/utils/sparse_mt.py` | Sparse MatrixTable operations |
| `gnomad/resources/resource_utils.py` | Resource classes (`TableResource`, `VersionedTableResource`, etc.) |
| `gnomad/resources/grch38/gnomad.py` | GRCh38 resource paths and constants |
| `gnomad/resources/grch37/gnomad.py` | GRCh37 resource paths and constants |
| `gnomad/sample_qc/` | Sample QC (ancestry, relatedness, sex, platform) |
| `gnomad/variant_qc/` | Variant QC (random forest, evaluation, training) |
| `gnomad/assessment/` | Summary stats, validity checks |

## Code Style

### Formatting

Code is formatted with **black** (preview mode, line length 88), **isort** (profile `"black"`), and **autopep8** (aggressive=1, ignoring E201/E202/E203/E731). Linting uses **pylint** and **pydocstyle** (PEP 257 convention, ignoring D105/D107). Config is in `pyproject.toml`, `.pydocstylerc`, and `.pylintrc`.

Pre-commit hooks are installed. Run `pre-commit run --all-files` to check.

```bash
# Manual formatting
black gnomad/ tests/
isort --profile black --filter-files gnomad/ tests/
```

### Docstrings

Use **Sphinx-style** (`:param:`, `:return:`) docstrings:

- **Summary line**: Concise one-line description, followed by blank line.
- **Body**: Extended description if needed. Use `.. note::` for caveats.
- **Params**: `:param name: Description.` — include default values and when params are conditionally required.
- **Return**: `:return: Description.` — describe the structure, not just the type.
- **Code references**: Use double backticks (````field_name````).
- **Constants**: Document with a docstring on the line immediately after the assignment.

```python
COVERAGE_CUTOFF = 30
"""Minimum median exome coverage differentiating high from low coverage sites."""


def my_function(
    ht: hl.Table,
    max_af: Optional[float] = None,
) -> hl.Table:
    """Short summary of what the function does.

    Extended description with more detail about behavior, edge cases, or
    design decisions.

    .. note::

        Any important caveats go here.

    :param ht: Input Table with ``freq`` and ``context`` fields.
    :param max_af: Maximum allele frequency threshold. Default is None (no
        cutoff).
    :return: Table with ``expected_variants`` and ``observed_variants`` fields
        added.
    """
```

### Type Annotations

- **All functions** must have type annotations on parameters and return values.
- Use `typing.List`, `typing.Optional`, `typing.Union`, `typing.Tuple`, `typing.Dict` for generic types.
- For Hail expressions, use the `hl.expr.*` prefix: `hl.expr.StructExpression`, `hl.expr.BooleanExpression`, `hl.expr.Float64Expression`, `hl.expr.Int32Expression`, `hl.expr.ArrayExpression`, `hl.expr.NumericExpression`, etc.
- For Hail table/matrix types: `hl.Table`, `hl.MatrixTable`.
- Use `Tuple[str, ...]` for variable-length tuples (not `Tuple[str]` which means exactly one element).
- Never use mutable defaults (`List`, `Dict`) in function signatures — use `None` and assign inside the function body.

### Function Design

- **Expression-based where possible**: Prefer functions that take and return Hail expressions rather than Tables. This makes them composable and reusable in different contexts. Use `expr._indices.source` to recover the source Table when needed internally.
- **Single Table param named `ht`**: When a function takes one Table, name it `ht`. Use descriptive names only when taking multiple Tables (e.g., `mutation_ht`, `gencode_ht`).
- **Pure transformations**: Utils functions should be HTs in / HTs out. All file I/O (read, write, checkpoint) belongs in pipeline scripts, not utility functions. Exceptions exist for historical reasons but should not be added.
- **No lazy imports**: Always use top-level imports unless needed to resolve circular imports.
- **`Optional` for nullable params**: Always wrap with `Optional[...]` when the default is `None`.

### Testing

- **Policy**: Any new or modified function in a PR must have tests.
- **Format**: Class-based tests with `pytest`. One test class per function.
- **Test data**: Use `hl.Table.parallelize()` to create small inline test tables.
- **Fixtures**: Use `@pytest.fixture` for shared setup. Hail init is handled by a session-scoped fixture in `tests/conftest.py`.
- **Docstrings**: Both test classes and test methods should have docstrings.

```python
class TestMyFunction:
    """Test the my_function function."""

    def test_basic_case(self):
        """Test that basic input produces expected output."""
        ht = hl.Table.parallelize(
            [{"x": 1, "y": 2.0}],
            hl.tstruct(x=hl.tint32, y=hl.tfloat64),
        )
        result = ht.annotate(z=my_function(ht.x, ht.y)).collect()[0]
        assert result.z == 3.0
```

## Hail Best Practices

### When to checkpoint vs cache

- **`checkpoint(new_temp_file(...))`**: Use for intermediate results that feed into multiple downstream operations or that follow expensive computations (joins, aggregations). Materializes to disk and breaks the query plan, preventing Hail from re-executing the upstream DAG.
- **`.cache()`**: Use for small results that will be reused immediately. Keeps data in memory but doesn't break the query plan as reliably.
- **After checkpoint, `.count()` is free**: It reads materialized metadata rather than re-executing the query. Place `.count()` after a checkpoint when you need the count.

### Avoid `.count()` for logging

Never use `.count()` just to log how many rows a table has. On large tables this forces full materialization and can cause Spark shuffle failures. Only call `.count()` when the result is actually needed for computation.

### `naive_coalesce()` after aggressive filters

When filtering a large table down to a small subset (e.g., LOFTEE HC LoF from all variants), most partitions become empty. This causes shuffle skew in downstream `group_by` aggregations. Call `.naive_coalesce(N)` after the filter to rebalance.

### `_localize=False` for expression results

Use `ht.aggregate(expr, _localize=False)` when the result will be used as a Hail expression in downstream operations (e.g., annotating another table). This avoids collecting to Python and back.

### Field existence checks

Use `field_name in ht.row` to check if a field exists on a Table. Hail Tables do not have `.get()`.

### `hl.or_else`, `hl.or_missing`, `hl.is_defined`

- **`hl.or_else(expr, default)`**: Substitute `default` when `expr` is missing.
- **`hl.or_missing(condition, expr)`**: Return `expr` when `condition` is True, missing otherwise.
- **`hl.is_defined(expr)`**: Returns True/False (never missing) — no need to wrap in `hl.or_else`.
- **`divide_null(num, denom)`**: Safe division returning null when denominator is 0. Import from `hail.utils.misc`.

### `approx_quantiles` is approximate

`hl.agg.approx_quantiles` uses the t-digest algorithm and returns approximate percentiles. Document this with a `.. note::` block when using it in functions.

### Falsy value gotchas

When checking optional numeric parameters, always use `is not None` instead of truthiness checks. `if max_af:` will skip `max_af=0.0`, which is a valid value. Same applies to any numeric parameter where 0 is meaningful.

### Array schema uniformity

All elements of a Hail array field must have identical struct schemas. You cannot annotate only `array[0]` with extra fields while leaving `array[1+]` unchanged — Hail will reject the mixed schema. Promote such metadata to the parent struct level instead.

### Rank assignment with `order_by`

`ht.order_by(expr)` destroys the key. To rejoin ranked results:
1. `ht.add_index("_rank_idx")` before ordering.
2. `rank_ht.key_by("_rank_idx")` after ordering.
3. Use `hl.scan.count()` to assign 0-based ascending ranks: `ht.order_by(val).annotate(rank=hl.scan.count())`.

### Small table reconstruction

`hl.Table.parallelize(hl.eval(ht.my_array_global), schema=...)` reconstructs a small Hail Table from a global array without re-running any jobs.

## Key Resource Classes

```python
from gnomad.resources.resource_utils import (
    TableResource,           # .ht() to read
    MatrixTableResource,     # .mt() to read
    VersionedTableResource,  # .ht() with version param
    VariantDatasetResource,  # .vds() to read
)
```

Usage:
```python
resource = TableResource(path="gs://gnomad/v4.1/constraint/metrics.ht")
ht = resource.ht()  # Reads the table
```

## Key Constraint API

```python
from gnomad.utils.constraint import (
    # Mutation rate
    annotate_mutation_type,
    annotate_with_mu,
    calibration_model_group_expr,
    # Consequence grouping
    build_constraint_consequence_groups,
    # Counting
    single_variant_count_expr,
    single_variant_observed_and_possible_expr,
    count_observed_and_possible_by_group,
    get_counts_agg_expr,
    # Model building & application
    build_models,
    apply_models,
    aggregate_expected_variants_expr,
    # GERP
    calculate_gerp_cutoffs,
    # Constraint metrics
    compute_pli,
    oe_confidence_interval,
    calculate_raw_z_score,
    calculate_raw_z_score_sd,
    get_constraint_flags,
    # Ranking & binning
    rank_and_assign_bins,
    compute_percentile_thresholds,
    annotate_bins_by_threshold,
)
```

## Maintaining CLAUDE.md

When working on any gnomAD repo that has a CLAUDE.md file, proactively add useful information you discover during development — gotchas, non-obvious API behavior, schema quirks, resource path conventions, or anything else that would save a future developer (or Claude session) time. Keep additions concise and placed in the appropriate section.

## CI/CD

- **Pre-commit hooks**: black, autopep8, pydocstyle, isort (run automatically on commit)
- **CI** (`.github/workflows/ci.yml`): Runs on push to main and PRs — black check, isort check, pydocstyle, autopep8 check, pylint, pytest
- **Publishing** (`.github/workflows/publish.yml`): Triggered by version tags (`v[0-9]+.[0-9]+.[0-9]+`) — validates version, runs tests, publishes to PyPI
