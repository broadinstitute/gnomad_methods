# gnomAD & Hail reference — institutional knowledge

Reference doc for anyone — AI or human — working with gnomAD data + Hail,
maintained by the gnomAD methods group. It captures the kind of knowledge
that usually only gets learned by running into the problem, or by being
told by someone who already did:

- Classes of bug that **produce plausible-looking wrong output rather
  than crashes** (the original motivation — silent-fail traps).
- **Design considerations** for gnomAD-touching analyses — when to
  prefer one pattern over another, tradeoffs that aren't obvious from
  the code.
- **Useful-but-hidden facts** about the datasets, HT schemas, and
  cross-repo conventions — known to the methods group by folklore but
  not easy to find from source.

Consult it before you touch gnomAD annotation tables, write Hail
aggregations, filter by QC status, assume coverage carries over between
HTs, or choose between two reasonable-sounding approaches. Prefer the
listed patterns when writing new code, and use the symptoms and design
notes as a checklist when scoping or debugging.

---

## Which sub-doc do I need?

Load only the one you need — that's the point of the split.

**Hail behavior**

| sub-doc | consult it when |
|---|---|
| [Missingness](missingness.md) | you touch a nullable field, a `hl.case` / `hl.if_else`, a `group_by` on a derived key, or a max/min over floats that could be NaN — or a total comes out smaller than expected with no error, which is usually a null drop |
| [Hail idioms](hail-idioms.md) | you rebind `ht`, build a big literal or many-branch expression, join with `ht[key]`, sort with `order_by`, use scans, pass a Python function to `.map()` or need real Python per row, or call an `hl.experimental` numeric kernel |
| [Lazy evaluation and materialization](lazy-eval-and-materialization.md) | you add a logging `.count()`, or you're choosing between `cache` / `persist` / `checkpoint` |
| [Partitioning, shuffles, and OOM](partitioning-shuffles-oom.md) | a job dies in a shuffle, you're about to `repartition`, or you want to know what shuffles in the first place |
| [Hidden and undocumented APIs](hidden-apis.md) | you meet an underscore-prefixed Hail parameter in gnomad_qc, or want the `hl.experimental` helpers we use |
| [VDS](vds.md) | you subset or densify a VDS, or you're choosing between Table / MatrixTable / VDS |
| [Hail versions](hail-versions.md) | you pick, upgrade, or downgrade the Hail version for a pipeline, or an older Hail can't read a table |

**gnomAD data**

| sub-doc | consult it when |
|---|---|
| [gnomAD-specific](gnomad-specific.md) | you filter variants, index into `freq`, cross HTs with different variant sets, touch non-release samples, or densify the v4 VDS and classify genotypes |
| [Biology](biology.md) | reference builds, multiallelics, sex chromosomes, trios, AC/AF null semantics, or what QUAL and LOFTEE do and don't tell you |
| [Release conventions and data policy](dataset-conventions.md) | you wonder why the data is shaped the way it is — annotation versions, subsets, constraint cutoffs across releases, labels, requester-pays |
| [Codebase and analysis patterns](analysis-patterns.md) | you pick a starting table, a VEP transcript, a bucket scheme, or a per-stratum aggregation, or read an HT whose filter provenance is unclear |

**Running the work**

| sub-doc | consult it when |
|---|---|
| [Dataproc / operational](dataproc-operational.md) | you run, size, package, or post-mortem a job on Dataproc or QoB — especially one that died with no Python traceback — or you're wondering why storage costs what it does |
| [Testing](testing.md) | after any refactor of the join graph or aggregation semantics, and before you pay for a cluster run |
| [Working with an AI assistant](working-with-ai.md) | you're setting up `CLAUDE.md` / `CLAUDE.local.md`, or wondering what's burning tokens and cluster time |

Entries are tagged **`[Hail]`** (universal Hail behavior) and/or
**`[gnomAD data]`** (dataset- or schema-specific). Entries that only the
internal team can use also carry one or both of:

- **`[internal: production]`** — relies on QC outputs or production
  pipeline artifacts (variant-QC filter HTs, sample-QC metadata).
- **`[internal: research]`** — relies on internal data but not on QC
  outputs (the full VDS, non-releasable samples).

Everything without an internal tag applies to anyone working with
released gnomAD data. Internal code spans two repos that
depend on each other — `gnomad_qc` imports `gnomad_methods` — so check a
change to one against the other; see "Check every change against
`gnomad_qc`" in the repo's `CLAUDE.md`.

---

## How to contribute

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
  sub-document fits (see the routing tables above). If nothing fits, a
  new sub-document is fine — add it to the matching routing table and to this
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
  **`[internal: production]`** and/or **`[internal: research]`** for
  anything only the internal team can use (definitions above).
- **Hail version** — whenever the behavior is version-dependent, state
  the version it was verified on (`Verified on Hail 0.2.134`). Behavior
  noted for one version can change. A bug or regression tied to specific
  Hail versions goes in [Hail versions](hail-versions.md), with a link to
  its source.

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
whose behavior can silently change under you — and update
[Hail versions](hail-versions.md): read the new release's change log for
problems that affect gnomAD, and check whether `FileFormat.version` changed
in the source, since a format bump isn't always listed in the change log's
File Format notes. Also re-read it when a gnomAD release changes schema
(v4.1 → v5).

---

<!-- Docs-site navigation for the sub-docs (hidden on the site; GitHub shows it
     as a plain list). Keep it in sync with the routing tables above. -->

```{toctree}
:hidden:
:maxdepth: 1

missingness.md
hail-idioms.md
lazy-eval-and-materialization.md
partitioning-shuffles-oom.md
hidden-apis.md
vds.md
hail-versions.md
gnomad-specific.md
biology.md
dataset-conventions.md
analysis-patterns.md
dataproc-operational.md
testing.md
working-with-ai.md
```
