# gnomad_methods Project Reference

## Project Overview

Shared Hail utility library for gnomAD pipelines, published to PyPI as the
`gnomad` package. Provides reusable functions for variant QC, sample QC,
constraint analysis, Ensembl VEP processing, resource management, and general
genomics operations. Used as a dependency by `gnomad_qc`, `gnomad_constraint`,
and other gnomAD repos — treat the public API as something downstream code
depends on.

**This is a library, not a pipeline.** It exposes APIs that other repos import;
it has no `__main__` entry point and is not run via Dataproc directly.
Public-API changes ripple through every consuming repo, so prefer additive
changes (new functions, new optional parameters with safe defaults) over
breaking renames.

**This is a public repo of genomic utility functions**, and the point of it is
reusability: future gnomAD team members and external users need to find and
understand code they can reuse. External groups have repeatedly asked the gnomAD
team for functionality that already existed here but that they couldn't find.
Optimize new code for discoverability — clear names, tight scopes, accurate
docstrings — not for cleverness.

## Repo Layout

Stable subpackage-level map (do not list individual functions here — they
change; see "Finding code" below):

| Directory | Purpose |
|-----------|---------|
| `gnomad/utils/` | General-purpose utilities: annotations, filtering, Ensembl VEP, constraint, sparse MTs, liftover, VCF export, etc. One module per topic. |
| `gnomad/resources/` | Resource classes (`resource_utils.py`), resource source config (`config.py`), and per-build resource definitions (`grch37/`, `grch38/`). |
| `gnomad/sample_qc/` | Sample QC: genetic ancestry, relatedness, sex inference, platform inference, filtering. |
| `gnomad/variant_qc/` | Variant QC: random forest, training, evaluation, LD. |
| `gnomad/assessment/` | Release assessment: summary stats, validity checks. |
| `tests/` | Pytest suite, mirroring the `gnomad/` layout. |
| `docs/` | Sphinx docs; API reference is auto-generated from docstrings. |

### What belongs here vs. in gnomad_qc

Generalized, reusable functions belong in gnomad_methods; gnomAD release
pipeline code belongs in `gnomad_qc`. Thin wrapper functions are acceptable
here.

The boundary is genuinely confusing in practice, and newcomers routinely search
both repos to find one thing. It is also not always clean: in places
gnomad_methods calls into `gnomad_qc`, where that `gnomad_qc` function is itself
mostly a wrapper around a gnomad_methods function. Trace the call chain before
assuming which repo owns a behavior.

## Finding Code

**Do not guess — read and confirm.** Do not assume a function exists, guess its
module from memory, or reason about its performance from its name. Function
names and locations change between versions.

This is the most common way work here goes wrong. A concrete example: asked to
make a gnomad_methods function more efficient, Claude proposed several changes
and asserted they would speed things up; none did, because it had not read the
function body or how the caller was using it. Read the actual implementation and
the actual call site before proposing a change, and don't claim a performance
win you haven't measured.

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

- **Summary line**: starts on the line *after* the opening `"""`, concise and
  one line, then a blank line. (Both styles exist in the repo; the
  summary-on-its-own-line form is overwhelmingly dominant — match it.)
- **Body**: extended description if needed. Use `.. note::` for caveats.
- **Params**: `:param name: Description.` — include defaults and when params
  are conditionally required.
- **Return**: `:return: Description.` — describe the structure, not just the
  type.
- **Do not document `:raises:`** — `:param:` and `:return:` only.
- **Code references**: use double backticks (``` ``field_name`` ```).
- **Constants**: document with a docstring on the line immediately after the
  assignment.

Be thorough but not wordy. Docstrings here are the public documentation, and an
over-long one is harder to follow than a tight one — favor coherence over
exhaustiveness.

```python
COVERAGE_CUTOFF = 30
"""Minimum median exome coverage differentiating high from low coverage sites."""


def my_function(
    ht: hl.Table,
    max_af: Optional[float] = None,
) -> hl.Table:
    """
    Short summary of what the function does.

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

- **All functions** must have type annotations on parameters and return values.
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
- **Tight scope, modest size**: a function should do one clearly nameable
  thing. Split large functions rather than adding parameters to them.
- **Readable over dense**: clear, easy-to-follow code beats compact code unless
  density is actually required for performance. Avoid long chains of Hail method
  calls — they are tempting but hard to read, hard to debug, and obscure where
  a shuffle or an evaluation is being triggered.
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
- **Each test must earn its place**: before adding a test, read the existing
  ones for that function and cover a case they don't. Near-duplicate tests cost
  runtime (Hail init is slow) and add no signal.
- **Test data**: build small simulated tables with `hl.Table.parallelize()` or
  `hl.utils.range_table()`. Avoid reading data from GCS in tests wherever
  possible.
- **Fixtures**: `@pytest.fixture` for shared setup. Hail init is handled by a
  session-scoped autouse fixture in `tests/conftest.py`.
- **Small tests do not prove scale**: passing on a toy table says nothing about
  how the function behaves on a real gnomAD dataset. Out-of-memory failures,
  shuffle failures, and partition skew only appear at scale. Do not describe a
  change as verified for production on the basis of unit tests alone.
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

### Lazy evaluation is the root of most surprises

Hail builds a query plan and executes nothing until something forces it. **Any
operation that converts Hail data into Python forces evaluation of the entire
upstream plan** — `.count()`, `.collect()`, `hl.eval()`, `.show()`,
`.aggregate()`, `.take()`. Before adding one of these, know what it will cause
to re-execute.

The most common mistake this produces: filter, then `.count()` for a log
message, then write. That runs the whole upstream computation twice. **Write or
checkpoint first, then log against the materialized result** — after a
checkpoint or write, `.count()` reads materialized metadata and is cheap.

### checkpoint vs cache

- **`checkpoint(new_temp_file(...))`**: for intermediate results feeding
  multiple downstream operations or following expensive computations (joins,
  aggregations). Materializes to disk and breaks the query plan so Hail won't
  re-execute the upstream DAG.
- **`.cache()`**: for small results reused immediately; doesn't break the
  query plan as reliably.
- **After a checkpoint, `.count()` is free** — it reads materialized metadata.

**Checkpoints are not free.** They cost compute and storage, and on gnomAD-scale
data that I/O is expensive. Callers of a library function generally don't know it
is checkpointing on their behalf, so a checkpoint buried in a utility can quietly
consume their storage. Look for places where a large function would genuinely
benefit — tables being joined, or several expensive computations that should be
forced once — and add checkpoints there, not everywhere. If a function needs a
`.count()` or another validity check that forces evaluation, place it after a
checkpoint when the upstream computation is expensive.

### Avoid `.count()` for logging

Never call `.count()` just to log row counts — on large tables it forces full
materialization and can cause Spark shuffle failures. Only count when the
result is needed for computation.

### Avoid shuffles

Shuffles are expensive, and shuffle failures are hard to diagnose because the
stack trace rarely points at the operation that caused them. Know which
operations shuffle — **re-keying a table shuffles it** — and avoid them where the
same result is reachable another way.

### Partitioning: `naive_coalesce` only, never `repartition`

Partitioning is an unsolved pain point in Hail. Rules of thumb:

- Large datasets need many partitions; small datasets do not.
- Joins benefit from more partitions. Joining a table with many partitions to
  one with few is pathologically slow — match them up first.
- Filtering a large table to a small subset leaves most partitions empty,
  causing shuffle skew in downstream `group_by` aggregations. Use
  `.naive_coalesce(N)` after the filter to rebalance.

**Never use `repartition()`** — neither to reduce nor to increase partition
count. To reduce, use `naive_coalesce()` (it avoids a shuffle). To increase,
repartition on read instead (e.g. the `_n_partitions` / `min_partitions`
argument on the read call), not on an in-memory dataset.

> **Caveat on `naive_coalesce` in pipelines**: in gnomAD pipelines
> `naive_coalesce` often runs **very** early, which means a large dataset gets
> collapsed into a small number of partitions up front and every downstream step
> pays for it. `naive_coalesce` is the right tool immediately after an aggressive
> filter; it is the wrong tool as a blanket early-pipeline call. Check where in
> the data flow the call actually lands before adding one.

### Aggregations are costly

Each aggregation pass over a large dataset costs real money. If a function needs
many aggregations, build a **single aggregator expression** (one `hl.struct` of
`hl.agg.*` expressions) and do one pass, rather than aggregating repeatedly.

### Always set the reference genome explicitly

Hail still defaults to **GRCh37**. gnomAD v3+ work is GRCh38. Set the reference
genome explicitly on `hl.init()` / locus construction / import calls rather than
relying on the default.

### `ht.aggregate(..., _localize=False)`

`_localize=False` keeps an aggregation result as a Hail expression instead of
returning it to Python, which avoids a round trip when the result feeds directly
into another expression.

Use it with caution: it is a private argument, it has been used in gnomAD code,
but the Hail team has advised against relying on it. **Confirm the current
recommendation with the Hail team before introducing new uses**, and prefer
supported alternatives where one exists.

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

## Execution and Infrastructure

This library has no entry point of its own; it runs inside whatever the consuming
code uses. In practice that is:

- `hailctl dataproc submit <cluster> script.py --param1 a --param2 b` — the most
  common path for gnomAD pipeline runs.
- Jupyter notebooks on a Dataproc cluster
  (`hailctl dataproc connect <cluster> nb`) — for exploratory work.
- Hail Batch / Query-on-Batch scripts — `python batch_script.py --param1 a`.
- Local Spark in local mode — for small tests.

gnomad_methods was designed by a group that primarily uses **Hail's Spark backend
on Google Cloud Dataproc**, and that assumption is baked into much of the code.
Functions should nonetheless work across backends and environments; don't add
code that only works on one.

### Cost and cluster conventions

- Run small tests **locally** with Hail's Spark backend in local mode. Data can
  be streamed directly from GCS via the GCS storage connector — no cluster
  needed.
- **Always read from the public gnomAD buckets where possible** to avoid
  requester-pays and egress charges.
- Use **us-central1-b**.
- Prefer **autoscaling** clusters.
- Prefer **preemptible** workers for initial job runs.

### Version and configuration hazards

- **Hail files are backwards-compatible, not forwards-compatible.** Data written
  by a newer Hail version cannot be read by an older one. Check the Hail version
  before assuming a path is readable, and don't casually bump the version used to
  write shared data.
- **Hail performance differs between versions.** There was a major performance
  regression from 0.2.130 to 0.2.131. Some work legitimately requires a newer
  version for new features — weigh the tradeoff rather than assuming newest is
  best.
- **GCP/GCS argument conventions change between versions, sometimes with large
  cost consequences.** A past change to how requester-pays buckets are specified
  caused egress fees to be charged on *all* writes, not just those actually
  touching a requester-pays bucket. Treat cloud configuration changes as
  cost-affecting until proven otherwise.
- An apparent bug in this library is sometimes a Hail version mismatch in the
  local environment. Check the installed version before debugging deeply.

## Never Do This

- **Never delete or overwrite data outside temp buckets.** When testing, write
  to a new temp path; do not overwrite an existing Table/MatrixTable.
- **Never remove code without asking first**, and never delete files that are
  referenced or built by the repo.
- **Never use `repartition()`** — see the partitioning section above.
- **Never claim a performance improvement you have not measured.**

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
