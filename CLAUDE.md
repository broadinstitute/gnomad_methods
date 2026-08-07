# gnomad_methods Project Reference

## Project Overview

Shared Hail utility library for gnomAD pipelines, published to PyPI as the
`gnomad` package. Provides reusable functions for variant QC, sample QC,
constraint analysis, Ensembl VEP processing, resource management, and general
genomics operations. Used as a dependency by `gnomad_qc`, `gnomad-constraint`,
and other gnomAD repos — treat the public API as something downstream code
depends on. See the README for the list of sibling repos.

**This is a library, not a pipeline.** It exposes APIs that other repos import;
it has no `__main__` entry point and is not run via Dataproc directly.
Public-API changes ripple through every consuming repo, so prefer additive
changes (new functions, new optional parameters with safe defaults) over
breaking renames.

**Check every change against `gnomad_qc`.** Before finishing a change to a
public function — signature, return schema, field names, defaults — grep
`gnomad_qc` for callers and report what breaks. **Flag incompatible call sites
for the user rather than fixing them silently**; the fix belongs in a
coordinated `gnomad_qc` PR, and the user decides what that looks like. Breaking
legacy code (e.g. v2 pipelines that are no longer run) is acceptable, so say
which version a broken caller belongs to.

CI catches only a subset of this: the `gnomad_qc` job installs this branch and
runs `pylint --disable=R,C,W` against `gnomad_qc` **main**, so it reports errors
only (missing names, bad attributes) and cannot see runtime breakage such as a
changed return schema or a reordered positional argument. A green CI run is not
evidence that downstream code still works.

**This is a public repo of genomic utility functions**, and the point of it is
reusability: future gnomAD team members and external users need to find and
understand code they can reuse. External groups have repeatedly asked the gnomAD
team for functionality that already existed here but that they couldn't find.
Optimize new code for discoverability — clear names, tight scopes, accurate
docstrings — not for cleverness.

## Repo Layout

**The repo layout lives in the README**, along with the list of related gnomAD
repos — read it there rather than duplicating it here. If a change adds,
removes, or moves a subpackage, update the layout table in the README as part of
the same pull request.

### What belongs here vs. in gnomad_qc

Generalized, reusable functions belong in gnomad_methods; gnomAD release
pipeline code belongs in `gnomad_qc`. Thin wrapper functions are acceptable
here.

The dependency runs one way only: `gnomad_qc` imports gnomad_methods, and
gnomad_methods never imports `gnomad_qc`. The boundary is still confusing in
practice, though, and newcomers routinely search both repos to find one thing —
partly because a `gnomad_qc` function is often a thin wrapper around a
gnomad_methods function of a similar name. Trace the call chain before assuming
which repo owns a behavior.

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
- **Downstream usage**: the sibling repos listed in the README (most usefully
  `gnomad_qc`, often checked out at `../gnomad_qc`) show how functions are used
  in real pipelines. Check there before changing a signature.
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

### Dependencies

Adding a new third-party library is a change to the package's install
requirements, and every consuming repo inherits it. Before introducing one,
check that it is compatible with the versions already pinned in
`requirements.txt` (and `requirements-dev.txt` for test-only libraries) — the
Hail, pandas, and numpy pins are the ones that usually conflict. If a new
library needs a pin loosened or bumped, **make that dependency update part of
the same pull request**, and say so in the PR description so reviewers know the
install surface changed. Prefer an existing dependency, or a few lines of code,
over a new one.

## Pull Requests

**Run the `/review` skill on every pull request before requesting review from a
team member** (Ben Weisburd's skill:
https://github.com/bw2/claude-code-review-skill). It reviews the diff under
multiple models and cross-validates the findings. Work through what it reports,
and note in the PR description that it was run and what you did or did not act
on. Human review time is the scarce resource here — don't spend it on things the
skill would have caught.

**If a change alters the structure of the repo, update the repo layout table in
the README in the same pull request.** Adding, removing, renaming, or moving a
subpackage all count. The README is the single source of truth for the layout —
it is what new contributors and external users read first, and a table that has
drifted from the tree is worse than no table. This applies to new dependencies
too: if `requirements.txt` changed, say so in the PR description.

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
- **Check optional numeric params with `is not None`, never truthiness.**
  `if max_af:` silently skips `max_af=0.0`, and `0.0` is a meaningful cutoff.
  Same for `0`, and for empty lists that a caller passed deliberately.

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
- **No file I/O in utilities**: reading, writing, and checkpointing belong in
  pipeline scripts, not here. A utility transforms what it is handed and returns
  the result. Historical exceptions exist; don't add new ones. Checkpointing in
  particular is not free and is not yours to spend: a caller has no idea a
  library function is checkpointing on their behalf, so a checkpoint buried in a
  utility quietly consumes their storage and compute. If a function genuinely
  needs one — a table about to be joined several ways, or an expensive
  computation that must be forced once — take a path or a flag from the caller
  rather than deciding for them.
- **Lazy imports are the exception, not the rule.** Use top-level imports.
  Three cases justify an import inside a function, and they are all already
  present in the repo: resolving a circular import; an optional or heavy
  dependency that shouldn't be required to import the package (`skl2onnx`,
  `ga4gh`); and picking a build-specific resource module at runtime (importing
  `grch37` vs `grch38` reference data based on a reference-genome argument).
  Anything else goes at the top of the file.
- **Private functions (`_name`)**: the leading underscore means "not part of the
  public API" — [docs/directives.py](docs/directives.py) skips private members,
  so a `_name` function does **not** appear in the generated API reference.
  Use it only for helpers that are implementation details of their own module,
  and default to public otherwise: this library exists to be reused, and a
  function that is private is a function nobody outside the repo can find. In
  `gnomad/resources/`, path-construction helpers and `_import_*` functions are
  private by convention — follow that when adding to those modules. Never make a
  function private that `gnomad_qc` or another repo already calls.
- **Don't break downstream**: renaming or removing a public function is a
  breaking change for `gnomad_qc` and other consumers. See the Project Overview
  for what CI does and does not catch — coordinate renames with downstream PRs.

## Testing

- **Policy**: any new or modified function in a PR must have tests.
- **Format**: class-based pytest, one test class per function, docstrings on
  both classes and methods, and type annotations on test methods (`-> None`)
  and fixtures like anywhere else.
- **Don't over-document tests**: the docstring says what case the test covers
  and, if it isn't obvious, why that case matters. Don't narrate the body —
  `# create a test table` above `hl.Table.parallelize(...)` is noise. If a
  literal in the test data is doing real work (a boundary value, a deliberate
  missing field), a short comment on *that* earns its place.
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

    def test_basic_case(self) -> None:
        """Test that basic input produces expected output."""
        ht = hl.Table.parallelize(
            [{"x": 1, "y": 2.0}],
            hl.tstruct(x=hl.tint32, y=hl.tfloat64),
        )
        result = ht.annotate(z=my_function(ht.x, ht.y)).collect()[0]
        assert result.z == 3.0
```

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

### Buckets and requester-pays

Released gnomAD data exists in more than one place, and which copy you read
determines who pays:

- `gs://gcp-public-data--gnomad` — the free public mirror, hosted by Google
  Cloud Public Datasets. **Read from here.** No requester-pays, no project
  billing for the read.
- `gs://gnomad-public-requester-pays` — where the gnomAD Production Team
  *writes* release data, which then syncs to the mirror above. It is public but
  **requester-pays**: reading it bills the requester's project. Don't hardcode
  paths into this bucket in examples, tests, or docstrings.
- `gs://gnomad` and other internal buckets — not public; internal work only.

**Prefer the resource classes over any literal path.** `gnomad/resources/` picks
the source for you (see `config.py`: it defaults to Google Cloud Public Datasets
and can be pointed elsewhere with `GNOMAD_DEFAULT_PUBLIC_RESOURCE_SOURCE`), so
`some_resource.ht()` reads the free copy without the caller having to know any
of the above. Note that the literal paths stored in `gnomad/resources/grch37/`
and `grch38/` are `gnomad-public-requester-pays` URLs — that is the write
location, and the source config rewrites them on read. Copying one of those
strings out of the source and using it directly bypasses the rewrite and incurs
charges.

### Cost and cluster conventions

- Run small tests **locally** with Hail's Spark backend in local mode. Data can
  be streamed directly from GCS via the GCS storage connector — no cluster
  needed.
- Use region **us-central1** (that's where the data lives — reading it from
  another region incurs egress charges). The zone within it is flexible:
  `us-central1-a`, `-b`, `-c`, and `-f` are all fine, and switching zones is a
  reasonable response to a capacity error.
- **Always autoscale, and always use preemptible workers**, unless you have hit
  a problem that requires otherwise. The usual reason to fall back to
  non-preemptible workers is a job that keeps dying in a shuffle: preemption
  during a large shuffle forces recomputation and can turn into a job that never
  finishes.

Treat any change to cloud configuration as cost-affecting until proven
otherwise, and check the installed Hail version before concluding that something
in this library is broken.

## Never Do This

- **Never delete or overwrite data outside temp buckets.** When testing, write
  to a new temp path; do not overwrite an existing Table/MatrixTable.
- **Never remove code without asking first**, and never delete files that are
  referenced or built by the repo.
- **Never add a bare `repartition(n)` to a large dataset** — it defaults to
  `shuffle=True` and forces a full shuffle. Set the partition count at read time
  (`_n_partitions` on the read call) or use `naive_coalesce()` to reduce it.
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
`gnomad/resources/config.py`, and "Buckets and requester-pays" above for which
copy of the data that resolves to and who pays for it.

## CI/CD

- **Pre-commit hooks**: black, autopep8, pydocstyle, isort (plus yaml/
  whitespace checks).
- **CI** (`.github/workflows/ci.yml`), on PRs and pushes to main:
  1. Lint job: black, isort, pydocstyle, autopep8, pylint, pytest.
  2. Docs job: builds Sphinx docs with `-W` (docstring RST errors fail the
     build); publishes to GitHub Pages on push to main.
  3. `gnomad_qc` job: installs this branch and runs `pylint --disable=R,C,W`
     over `gnomad_qc` — catches breaking API changes, with the limits described
     in the Project Overview.
- **Releases**: see CONTRIBUTING.md.

## Maintaining CLAUDE.md

Watch for useful discoveries while working here — gotchas, non-obvious Hail
behavior, schema quirks, conventions — and **propose them to the user rather
than adding them unilaterally**. Say what you learned and which section it
belongs in; let the user decide whether it is durable enough to write down. Keep
this file lean:

- **No function inventories or import lists** — they go stale when functions
  are renamed or moved. Reference directories or modules at most, and prefer
  describing how to *find* things over enumerating them.
- Prefer facts that are expensive to rediscover (an afternoon lost to a Spark
  shuffle failure) over facts that are one grep away.
- Remove entries that are no longer true rather than accumulating corrections.
