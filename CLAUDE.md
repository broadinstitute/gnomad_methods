# gnomad_methods Project Reference

## Project Overview

Shared Hail utility library for gnomAD pipelines, published to PyPI as the
`gnomad` package. Provides reusable functions for variant QC, sample QC,
constraint analysis, VEP processing, resource management, and general genomics
operations. Used as a dependency by `gnomad_qc`, `gnomad_constraint`, and other
gnomAD repos — treat the public API as something downstream code depends on.

## Repo Layout

Stable subpackage-level map (do not list individual functions here — they
change; see "Finding code" below):

| Directory | Purpose |
|-----------|---------|
| `gnomad/utils/` | General-purpose utilities: annotations, filtering, VEP, constraint, sparse MTs, liftover, VCF export, etc. One module per topic. |
| `gnomad/resources/` | Resource classes (`resource_utils.py`), resource source config (`config.py`), and per-build resource definitions (`grch37/`, `grch38/`). |
| `gnomad/sample_qc/` | Sample QC: ancestry, relatedness, sex inference, platform inference, filtering. |
| `gnomad/variant_qc/` | Variant QC: random forest, training, evaluation, LD. |
| `gnomad/assessment/` | Release assessment: summary stats, validity checks. |
| `tests/` | Pytest suite, mirroring the `gnomad/` layout. |
| `docs/` | Sphinx docs; API reference is auto-generated from docstrings. |

## Finding Code

Do not assume a function exists or guess its module from memory — verify in
the repo. Function names and locations change between versions.

- **Search before writing**: before adding a new utility, grep the package for
  related keywords (`grep -rn "allele_frequency" gnomad/`). This library is
  large and the function you need often already exists, possibly under a
  different name than you'd guess.
- **Module docstrings and section headers**: each module starts with a
  docstring describing its scope — read it to confirm you're in the right
  place before diving into functions.
- **Tests show intended usage**: `tests/` mirrors the package layout; a test
  file is often the best usage example for a function.
- **Downstream usage**: sibling repos (e.g. `../gnomad_qc`) show how functions
  are used in real pipelines. Check there before changing a signature.
- **Generated docs**: https://broadinstitute.github.io/gnomad_methods/ — built
  from docstrings, so the docstring in the source is always authoritative.

## Development Commands

```bash
pip install -r requirements.txt -r requirements-dev.txt   # setup
python3 -m pre_commit install                             # install hooks

black gnomad tests                  # format
isort --profile black --filter-files gnomad tests
autopep8 --in-place gnomad          # comment formatting
pydocstyle gnomad tests             # docstring check
./lint                              # pylint (gnomad + tests)
python -m pytest                    # all tests
python -m pytest tests/utils/test_vep.py  # one file
./docs/build.sh                     # build docs (regenerates API reference)
```

Formatting config lives in `pyproject.toml`, `.pydocstylerc`, `.pylintrc`, and
`.pre-commit-config.yaml`. Black runs in preview mode with the default
88-character line length.

## Code Style

### Docstrings

Use **Sphinx-style** (`:param:`, `:return:`) docstrings. These are the public
documentation — the docs build renders them with Sphinx in `-W` mode, so
malformed RST in a docstring **fails CI**.

- **Summary line**: concise one-line description, then a blank line.
- **Body**: extended description if needed. Use `.. note::` for caveats.
- **Params**: `:param name: Description.` — include defaults and when params
  are conditionally required.
- **Return**: `:return: Description.` — describe the structure, not just the
  type.
- **Code references**: use double backticks (``` ``field_name`` ```).
- **Constants**: document with a docstring on the line immediately after the
  assignment.

```python
COVERAGE_CUTOFF = 30
"""Minimum median exome coverage differentiating high from low coverage sites."""


def my_function(
    ht: hl.Table,
    max_af: Optional[float] = None,
) -> hl.Table:
    """Short summary of what the function does.

    Extended description with detail about behavior or edge cases.

    .. note::

        Any important caveats go here.

    :param ht: Input Table with ``freq`` and ``context`` fields.
    :param max_af: Maximum allele frequency threshold. Default is None (no
        cutoff).
    :return: Table with ``expected_variants`` field added.
    """
```

### Type Annotations

- **All functions** must annotate parameters and return values.
- Use `typing.List`, `typing.Optional`, `typing.Union`, etc. for generics.
- For Hail expressions use the `hl.expr.*` prefix (`hl.expr.StructExpression`,
  `hl.expr.BooleanExpression`, ...); for tables use `hl.Table` /
  `hl.MatrixTable`.
- `Tuple[str, ...]` for variable-length tuples (`Tuple[str]` means exactly one
  element).
- Never use mutable defaults — use `None` and assign inside the body.
- Always wrap nullable params in `Optional[...]`.

### Function Design

- **Expression-based where possible**: prefer functions that take and return
  Hail expressions rather than Tables — they compose and reuse better. Use
  `expr._indices.source` internally to recover the source Table when needed.
- **Single Table param named `ht`**: use descriptive names only when a
  function takes multiple Tables (e.g. `mutation_ht`, `gencode_ht`).
- **Pure transformations**: utility functions should be HTs in / HTs out. File
  I/O (read, write, checkpoint) belongs in pipeline scripts, not utilities.
  Historical exceptions exist; don't add new ones.
- **No lazy imports**: top-level imports only, unless resolving a circular
  import.
- **Don't break downstream**: renaming or removing a public function is a
  breaking change for `gnomad_qc` and other consumers. CI installs your branch
  and pylints `gnomad_qc` against it — coordinate renames with downstream PRs.

## Testing

- **Policy**: any new or modified function in a PR must have tests.
- **Format**: class-based pytest, one test class per function, docstrings on
  both classes and methods.
- **Test data**: build small inline tables with `hl.Table.parallelize()`.
- **Fixtures**: `@pytest.fixture` for shared setup. Hail init is handled by a
  session-scoped autouse fixture in `tests/conftest.py`.
- **Local Spark is forced**: `tests/conftest.py` sets
  `HAIL_QUERY_BACKEND=spark` and pins the Spark driver to 127.0.0.1 *before*
  Hail is imported, so tests never submit to Hail Batch regardless of your
  `hailctl` config. Don't reorder or remove that setup.
- Coverage is sparse — many modules have no tests yet. Mirror the `gnomad/`
  layout when adding test files.

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

### checkpoint vs cache

- **`checkpoint(new_temp_file(...))`**: for intermediate results feeding
  multiple downstream operations or following expensive computations (joins,
  aggregations). Materializes to disk and breaks the query plan so Hail won't
  re-execute the upstream DAG.
- **`.cache()`**: for small results reused immediately; doesn't break the
  query plan as reliably.
- **After a checkpoint, `.count()` is free** — it reads materialized metadata.

### Avoid `.count()` for logging

Never call `.count()` just to log row counts — on large tables it forces full
materialization and can cause Spark shuffle failures. Only count when the
result is needed for computation.

### `naive_coalesce()` after aggressive filters

Filtering a large table to a small subset leaves most partitions empty,
causing shuffle skew in downstream `group_by` aggregations. Call
`.naive_coalesce(N)` after the filter to rebalance.

### `_localize=False` for expression results

Use `ht.aggregate(expr, _localize=False)` when the result feeds back into Hail
expressions (e.g. annotating another table) — avoids a round trip through
Python.

### Missingness helpers

- `hl.or_else(expr, default)`: substitute `default` when `expr` is missing.
- `hl.or_missing(condition, expr)`: `expr` when condition is True, else
  missing.
- `hl.is_defined(expr)`: returns True/False, never missing — no `hl.or_else`
  wrapper needed.
- `divide_null(num, denom)` (from `hail.utils.misc`): safe division, null when
  denominator is 0.

### Field existence checks

Use `field_name in ht.row` — Hail Tables have no `.get()`.

### Falsy value gotchas

Check optional numeric params with `is not None`, never truthiness:
`if max_af:` silently skips `max_af=0.0`.

### Array schema uniformity

All elements of a Hail array field must share an identical struct schema. You
can't annotate only `array[0]` with extra fields — Hail rejects the mixed
schema. Promote such metadata to the parent struct.

### Rank assignment with `order_by`

`ht.order_by(expr)` destroys the key. To rejoin ranked results:
`ht.add_index("_rank_idx")` before ordering, `key_by("_rank_idx")` after, and
use `hl.scan.count()` for 0-based ascending ranks.

### `approx_quantiles` is approximate

`hl.agg.approx_quantiles` uses t-digest and returns approximate percentiles —
document this with a `.. note::` when using it.

### Small table reconstruction

`hl.Table.parallelize(hl.eval(ht.my_array_global), schema=...)` rebuilds a
small Table from a global array without re-running jobs.

## Resources

`gnomad/resources/resource_utils.py` defines the resource class hierarchy
(`TableResource`, `MatrixTableResource`, `VersionedTableResource`, etc.).
Resources wrap a cloud path and expose a reader (`.ht()`, `.mt()`, `.vds()`);
versioned resources hold a dict of per-version resources with a default.
Build-specific resource definitions live in `gnomad/resources/grch37/` and
`grch38/`.

Public resources can be read from multiple cloud sources (gnomAD GCS buckets,
Google Cloud Public Datasets, AWS Open Data); the source is auto-detected or
overridden via the `GNOMAD_DEFAULT_PUBLIC_RESOURCE_SOURCE` env var — see
`gnomad/resources/config.py`.

## CI/CD & Releases

- **Pre-commit hooks**: black, autopep8, pydocstyle, isort (plus yaml/
  whitespace checks).
- **CI** (`.github/workflows/ci.yml`), on PRs and pushes to main:
  1. Lint job: black, isort, pydocstyle, autopep8, pylint, pytest.
  2. Docs job: builds Sphinx docs with `-W` (docstring RST errors fail the
     build); publishes to GitHub Pages on push to main.
  3. `gnomad_qc` job: installs this branch and runs `pylint --disable=R,C,W`
     over `gnomad_qc` — catches breaking API changes.
- **Releases**: bump `version` in `setup.py` (semver, based on changes since
  the last release), merge to main, then push a `v<X.Y.Z>` tag.
  `.github/workflows/publish.yml` validates the tag matches `setup.py` and
  publishes to PyPI. See CONTRIBUTING.md for details.

## Maintaining CLAUDE.md

When working in this repo, proactively add useful discoveries — gotchas,
non-obvious Hail behavior, schema quirks, conventions — to the appropriate
section. Keep it durable and lean:

- **No function inventories or import lists** — they go stale when functions
  are renamed or moved. Reference directories or modules at most, and prefer
  describing how to *find* things over enumerating them.
- Prefer facts that are expensive to rediscover (an afternoon lost to a Spark
  shuffle failure) over facts that are one grep away.
- Remove entries that are no longer true rather than accumulating corrections.
