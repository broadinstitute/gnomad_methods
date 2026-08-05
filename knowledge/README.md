# gnomAD & Hail reference — institutional knowledge

Reference doc for anyone — AI or human — working with gnomAD data + Hail.
Contents range across:

- Classes of bug that **produce plausible-looking wrong output rather
  than crashes** (the original motivation — silent-fail traps).
- **Design considerations** for gnomAD-touching analyses — when to
  prefer one pattern over another, tradeoffs that aren't obvious from
  the code.
- **Useful-but-hidden facts** about the datasets, HT schemas, and
  cross-repo conventions — the kind of thing that's known to the
  methods group by folklore but not easy to find from source.

Prefer the listed patterns and considerations when writing new analysis
code, and use the symptoms + design notes as a checklist when
scoping or debugging.

---

## Which sub-doc do I need?

Load only the one you need — that's the point of the split.

| sub-doc | consult it when |
|---|---|
| [Missingness](missingness.md) | you touch a nullable field, a `hl.case` / `hl.if_else`, a `group_by` on a derived key, or a max/min over floats that could be NaN |
| [Hail idioms](hail-idioms.md) | you rebind `ht`, build a big literal or many-branch expression, join with `ht[key]`, sort with `order_by`, use scans, or reach for a Python UDF |
| [Lazy evaluation and materialization](lazy-eval-and-materialization.md) | you add a logging `.count()`, or you're choosing between `cache` / `persist` / `checkpoint` |
| [Partitioning, shuffles, and OOM](partitioning-shuffles-oom.md) | a job dies in a shuffle, you're about to `repartition`, or you want to know what shuffles in the first place |
| [Hidden and undocumented APIs](hidden-apis.md) | you meet an underscore-prefixed Hail parameter in gnomad_qc, or want the `hl.experimental` helpers we use |
| [gnomAD-specific](gnomad-specific.md) | you filter variants, index into `freq`, or cross HTs with different variant sets |
| [VDS](vds.md) | you subset or densify a VDS, or you're choosing between Table / MatrixTable / VDS |
| [Biology](biology.md) | reference builds, multiallelics, sex chromosomes, trios, AC/AF null semantics, or what QUAL and LOFTEE do and don't tell you |
| [Release conventions and data policy](dataset-conventions.md) | you wonder why the data is shaped the way it is — annotation versions, subsets, constraint cutoffs across releases, labels, requester-pays |
| [Codebase and analysis patterns](analysis-patterns.md) | you pick a starting table, a VEP transcript, a bucket scheme, or a per-stratum aggregation |
| [Dataproc / operational](dataproc-operational.md) | you run, size, package, or post-mortem a job on Dataproc or QoB, or you're wondering why storage costs what it does |
| [Testing](testing.md) | after any refactor of the join graph or aggregation semantics, and before you pay for a cluster run |
| [Working with an AI assistant](working-with-ai.md) | you're setting up `CLAUDE.md` / `CLAUDE.local.md`, or wondering what's burning tokens and cluster time |


Entries are tagged **`[Hail]`** (universal Hail behavior) and/or
**`[gnomAD data]`** (dataset- or schema-specific). Entries tagged
**`[internal only]`** describe data or artifacts available only to the
internal production team — unreleasable samples, the pre-release VDS,
internal QC tables. Everything untagged that way applies to anyone
working with released gnomAD data.

---

## For AI or human reading this file

**What this file is:** a shared reference maintained by the gnomAD
methods group. It captures data-model quirks, Hail idioms, codebase
conventions, and design-choice tradeoffs — the kind of institutional
knowledge that usually only gets learned by running into the problem
(or being told by someone who already did).

**When to consult it:** any time you're touching gnomAD annotation
tables, writing Hail aggregations, filtering by QC status, making
assumptions about coverage between HTs, or picking between two
reasonable-sounding analysis approaches. See the [When to consult this
](#when-to-consult-this) section at the bottom for the concrete trigger
list.

---

### How to contribute

**Rules — these come first** (and apply to AI assistants especially):

1. **Never guess. No claims without proof.** Verify every fact against
   the official Hail docs, the Hail source, or a tiny job you actually
   ran. Never reference a function that doesn't exist. If you can't
   verify it, don't add it.
2. **Propose, don't edit.** Surface candidate entries to the user and
   get a yes before writing. This doc bloats easily; only high-value or
   recurring issues belong.
3. **Be concise.** Tight entries, no padding. Prefer one line over a
   paragraph.
4. **Don't duplicate the official Hail docs.** Link them; capture only
   the non-obvious or anecdotal part.
5. **Tag every entry** (format below).
6. **Build the docs before opening the PR.** These pages are published
   to the documentation site, so CI renders them with
   `sphinx -W` — warnings are errors, and a malformed entry fails the
   build rather than merely looking odd on GitHub:

   ```sh
   pip install -r docs/requirements.docs.txt
   ./docs/build.sh
   ```

   Three things GitHub renders happily and Sphinx rejects:
   - **Adjacent `---` rules** with nothing between them.
   - **Skipping a heading level** — entries are `##` under the page's
     single `#` title. Don't reach for `###`.
   - **Non-Python placeholders in a ` ```python ` block** — `…` outside
     a string literal can't be lexed. Use `...`, or drop the fence's
     language tag.

**What clears the bar.** If during a session you learn something that
would help future readers — a bug you hit, a design choice you had to
reason through, a fact that took a while to piece together from code —
propose an entry. It should be:

1. **Non-obvious**: not readily discoverable from Hail / gnomAD docs or
   from a `grep` in the code. You either ran into it, or someone told
   you, or you had to reason your way to it.
2. **Generalizable**: applies beyond one specific analysis. If it's
   one-file-specific, it belongs in that file's docstring instead.
3. **Actionable**: a future reader can use it. Either "use this pattern
   instead of that one" or "know this fact before you make this
   choice" — not just background context.

**Item categories worth capturing:**
- **Silent-fail bugs** — the wrong pattern produces a plausible answer.
  Fixing these silently improves any analysis that would have hit them.
- **Design-choice tradeoffs** — e.g. "canonical vs MANE Select
  transcript picking, and when each is right"; "per-population EM vs
  full-cohort EM"; "when to drop OE-candidate × candidate pairs vs keep
  them." Documenting these lets a future analysis start with the right
  default.
- **Hidden-in-code facts** — things about HT schemas, path
  conventions, filter definitions, or upstream pipeline behavior that
  aren't apparent without reading the source. Reduces "why does this
  do that" debugging.
- **Cross-cutting conventions** — things like "adj is per-genotype and
  separate from variant-level filters"; "the v4 VDS is stored split";
  "test-mode output paths can differ between steps of the same
  pipeline, so a test run isn't necessarily self-consistent."

**Contribution format:** every item follows this shape:
- **Heading** (`##`) — one-line "the thing" summary, in whichever
  sub-document fits (see the routing table above). If nothing fits, a
  new sub-document is fine — add it to the routing table and to this
  README's `toctree` so the docs build can find it.
- **One-paragraph explanation** of the mechanism / rationale.
- **Code example** where relevant (fenced code block, small).
- **Symptom** for silent-fail items ("aggregation totals < input
  rowcount", etc.), OR **Tradeoff** for design items (what does each
  option cost / gain), OR **Where to find more** for hidden-fact items
  (file paths, upstream repo, etc.).
- **Recommendation / rule of thumb**.
- **Scope tag** — one of **`[Hail]`** (universal Hail behavior),
  **`[gnomAD data]`** (dataset/schema-specific), or both. Add
  **`[internal only]`** for anything that applies solely to the internal
  production team (unreleasable samples, internal bucket paths).
- **Hail version** — whenever the behavior is version-dependent, state
  the version it was verified on (`Verified on Hail 0.2.134`). Behavior
  noted for one version can change.

Add new items in whichever section fits best. If nothing fits, add a
new section — but tell the user first so we can decide whether that
section belongs.

Do **not** add:
- One-off analysis details (put in the analysis's docstring)
- Anything already well-covered in Hail docs
- Anything about specific code paths (belongs in code comments)
- Full explainers of biology / statistics / methods — this file is a
  quick-reference, not a textbook

**When to re-look at this doc:** on every Hail version bump, re-verify
the entries carrying a `Verified on Hail X` line — those are the ones
whose behavior can silently change under you. Also re-read it when a
gnomAD release changes schema (v4.1 → v5).

---

## When to consult this

Before writing new gnomAD-touching analysis code, especially:

- Anything that **filters variants** — see [gnomAD-specific](gnomad-specific.md#gnomad-specific)
  (filters, only_filters.ht, releasability).
- Anything with **`hl.case`, `hl.if_else`, set/dict ops on nullable
  fields** — see [Missingness](missingness.md#missingness).
- Anything that **group_by / agg_group_by on derived keys** — same.
- Anything that takes a **max / min over float data that can be NaN** —
  same section; `hl.max` and `hl.agg.max` disagree on NaN.
- Anything that builds a **many-branch expression or a per-stratum
  breakdown in one row op** (a `hl.case` with many `.when`s, or a dict
  comprehension of heavy sub-exprs) — the `ClassTooLargeException` item
  in [Hail idioms](hail-idioms.md#hail-idioms); prefer the loop-restrict-per-stratum
  pattern in [Codebase / analysis patterns](analysis-patterns.md#codebase-and-analysis-patterns).
- Anything that calls an **`hl.experimental` numeric kernel / iterative
  solver** (EM, statistical routines) — cast dtypes and assert value
  invariants first; see the un-validated-kernel item in
  [Hail idioms](hail-idioms.md#hail-idioms).
- Anything that adds a **logging `.count()` / `.show()`** to a pipeline
  — see [Lazy evaluation and materialization](lazy-eval-and-materialization.md#lazy-evaluation-and-materialization);
  it's free after a checkpoint and a full re-run before one.
- Anything that **indexes into `freq` arrays** — the freq[0] item.
- Anything that touches **trios, family data, or non-release samples**
  — releasability + trio gotchas.
- Anything that **picks a single VEP transcript** — the MANE / canonical
  item.
- Anything that **reads HTs whose filter provenance is unclear** — the
  bake-time-vs-read-time item.
- Anything that **subsets a VDS** (samples or sites) or densifies one —
  see [VDS](vds.md#vds); filtering one of the two MatrixTables alone corrupts
  the densified output silently.
- Anything involving **sex chromosomes, PAR, or multiallelics** — see
  [Biology](biology.md#biology).
- Anything that **densifies the v4 VDS and classifies genotypes**
  (het / hom-var counts, freq, phasing) — the high-AB correction +
  operation-order item under [gnomAD-specific](gnomad-specific.md#gnomad-specific).
- Anything you **run or debug on Dataproc** — log capture and
  post-mortems in [Dataproc / operational](dataproc-operational.md#dataproc--operational),
  especially when a job dies with no Python traceback.

And after any bug that looks like "the total is smaller than I
expected but there's no error message" — that's almost always a
`hl.case`/`hl.if_else`/`group_by` null drop.

```{toctree}
:hidden:
:maxdepth: 1

missingness.md
hail-idioms.md
lazy-eval-and-materialization.md
partitioning-shuffles-oom.md
hidden-apis.md
gnomad-specific.md
vds.md
biology.md
analysis-patterns.md
dataproc-operational.md
dataset-conventions.md
testing.md
working-with-ai.md
```
