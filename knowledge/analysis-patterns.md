# Codebase and analysis patterns

## Row-set anchor determines what "could possibly appear"

**`[Hail]` + `[gnomAD data]`**

`hl.Table.filter()` removes rows; it can never add them. Similarly, an
"annotate" step's `.annotate(**...)` adds *columns* but not *rows* — the
row set is fully determined by the input HT's row set. Choice of the
starting HT (the "anchor") determines the outer bound on what downstream
code can see.

**Two ways this bites:**

1. **The anchor drops rows before any filter runs.** A sites table built
   by starting from the VEP table (`vep_ht.select("vep")`) and joining the
   filters onto it loses every variant that's in the filter table but
   missing from the VEP table — before any `.filter()` runs. Loosening the
   filter can't bring those variants back.

2. **A more permissive lookup table doesn't add rows.** Annotating a
   release-scoped table by looking up a more permissive one (say, one that
   also covers trio-only variants) still yields only the release-scoped
   rows. The trio-only rows were never in the anchor, so downstream
   left-joins see them as missing. The fix is to change the anchor — start
   from the table that contains the rows you want — not the lookup
   table.

**Symptom of getting this wrong:** you build a "permissive" annotation
path and downstream aggregations produce bit-identical output to the
release-scoped baseline. The permissive lookups happened, but the row
set never grew.

**Rule of thumb:** anchor on the widest-coverage table for the
semantics you want. Join everything else with nullability tolerance
downstream. When you have both a "release" and a "trio-inclusive"
consumer, you likely need two annotated HTs, one anchored on each row
set. Document the anchor choice in the function docstring so callers can
predict what's missing.

---

## Multi-source variants collapse to one "priority" tag

**`[gnomAD data]`**

Analyses that assign a single "priority" tag per variant (source-set →
priority-source) pick from a priority list. That priority list is a
design decision — not obvious from input tags.

Example: `SOURCE_PRIORITY = ["clinvar_plp", "hc_lof", "splice_path",
"clinvar_vus", "in_trans_oe_candidate", ..., "clinvar_blb"]`. A variant
tagged both `clinvar_blb` and `in_trans_oe_candidate` gets priority
= `in_trans_oe_candidate` (because `clinvar_blb` is at the *bottom* of
the list — a deliberate choice for that analysis).

**Symptom:** apparent tag counts don't sum to source-set counts.
**Fix:** document the priority order visibly in code and downstream
analysis narratives.

---

## VEP transcript picking: canonical vs MANE Select vs first vs explode-all

**`[gnomAD data]`**

A variant has many `transcript_consequences`. Different analyses pick
different representative transcripts:

- **Explode all** — one row per (variant, transcript). Right when you
  need per-transcript info.
- **Canonical** — `tc.canonical == 1`. Ensembl's designated canonical
  transcript per gene.
- **MANE Select** — `tc.mane_select` defined and non-empty. NCBI/EBI
  joint standard; generally preferred over canonical when available.
- **First** — `tcs[0]`. Order depends on VEP config; can be non-coding
  or a pseudogene transcript when the variant is intergenic.

The gotcha: picking `tcs[0]` when the first transcript is non-coding
gives a `gene_symbol = null` for a variant that clearly has a real
gene. Silently blanks out downstream per-gene tables. In one session
this caused MUC16/FBN3 to disappear from gene-symbol columns.

**Recommended:** filter to protein-coding Ensembl transcripts first,
then prefer canonical → MANE Select → first → `or_missing()`. Document
the pick order at the call site.

**The same trap applies to per-transcript pathogenicity scores.** REVEL
and similar predictors emit a score **per transcript**, not per variant.
Collapsing them with `max` is common and quietly picks whichever
transcript scored highest — which is frequently not the canonical or MANE
transcript, and not the one your consequence annotation came from
(*reported · project meeting*). If your consequence and your score are
chosen by different rules, they describe different transcripts. Pick the
transcript first, then take that transcript's score.

---

## `missing_*` / "other" default buckets: check what they actually catch

**`[Hail]` + `[gnomAD data]`**

A `hl.case().default("missing_af")` (or `"unknown"`, `"other"`, etc.) is
a *catch-all* for rows that didn't match any preceding `.when()`. It
does **not** mean "the AF value was null" — it means "the AF didn't fit
any of the buckets you defined."

A concrete failure:

```python
AF_BUCKETS = [
    ("Singleton", ...),
    ("(0, 1e-4)", 0.0, 1e-4),
    ("[1e-4, 1e-3)", 1e-4, 1e-3),
    ("[1e-3, 1e-2)", 1e-3, 1e-2),
    ("[1e-2, 5e-2)", 1e-2, 5e-2),  # ← ceiling at 5e-2
]
# Any variant with AF ≥ 5e-2 falls through to the default:
case.default("missing_af")
```

If the input includes higher-AF variants than the buckets were designed
for (e.g. a candidate list expanded up to AF 0.5), those common variants
land in `missing_af` — the same bucket as genuinely null-AF variants. The
label is misleading: in one analysis, 64% of the inputs landed there and
were read as a coverage gap, when they were common variants above the
ceiling of a rare-variant binning scheme.

**Rule of thumb:**
- Name catch-all default buckets by what they *catch* ("out_of_range",
  "af_gte_5e-2"), not by what you *hope* they catch ("missing_af").
- If a plausible input range extends past your bucket ceiling, add
  an explicit bucket for it rather than relying on the default.
- Before interpreting a large "missing" bucket as a coverage gap, sanity-
  check by looking at the actual input distribution. Ideally emit
  separate default and null buckets:
  ```python
  case = case.when(hl.is_missing(af_expr), "null_af")
  # ... range buckets ...
  return case.default("out_of_range")
  ```

**Symptom:** a large fraction of your rows end up in the default
bucket, and interpreting the fraction as a coverage / annotation-scope
issue leads you down a rabbit hole that's actually about bucket-range
design.

---

## Import paths from a resource module; don't type bucket paths in analysis code

**`[gnomAD data]`**

Dataset paths belong in a versioned resource module, not in the script
that uses them:

```python
from gnomad_qc.v3.resources import gnomad_v3_genotypes_vds   # good
vds_path = "gs://.../v3.1/raw/gnomad_v3.1.vds"               # bad
```

Two reasons, and the second is the one that bites. First, everyone
resolves the same version of the same artifact. Second, a hardcoded path
is invisible to the *next* release: when a dataset is re-cut, moved, or
superseded, the resource module gets updated and every consumer follows,
while pasted paths keep silently reading last year's data — a job that
still runs, on the wrong input.

This applies to this doc too: entries here name artifacts by what they
are and how to reach them through a resource function, not by bucket
path.

**Related, on cost rather than correctness:** before a large or
expensive run, it's worth asking someone who has run something similar
what cluster configuration they used. Compute for a genome-scale job is
easy to over- or under-provision by an order of magnitude, and both
directions cost real money.

---

## Aggregation-time vs bake-time filters — both need documenting

**`[gnomAD data]`**

If your on-disk HT has 1M rows but a downstream analysis shows 700k
rows, you need to know *where* the filter was applied — inside a
Dataproc job that wrote the HT, or inside the loader at read time.

Common example: mega-gene exclusion via `--exclude-gene-ids` gets baked
into some HTs at build time, but not others. Downstream loaders may
re-apply the exclusion at read time. `ht.count()` on the raw on-disk
HT gives one number; `_load_trio_ht(args).count()` gives another.
Neither is wrong.

**Rule of thumb:** in the docstring of any function that returns an
HT, list the filters applied. When passing HT paths as CLI overrides,
expect that the caller may not know what filters are baked in.

---

## Importing gnomad turns off your script's INFO logging

**`[gnomAD data]` (gnomad_methods code) · Verified on Python 3.11 with `gnomad.utils.filtering`**

27 gnomad_methods modules call `logging.basicConfig(...)` when imported, which
installs a handler on the root logger. A script that imports any of them and
then calls `logging.basicConfig(level=logging.INFO)` gets nothing:
`basicConfig` does nothing once the root logger has a handler, so the root
level stays at WARNING and the script's own INFO messages disappear.

```python
import gnomad.utils.filtering                      # installs a root handler
logging.basicConfig(level=logging.INFO)            # silently ignored
logging.basicConfig(level=logging.INFO, force=True)  # works (Python 3.8+)
logger.setLevel(logging.INFO)                      # also works: the gnomad_qc pattern
```

**What the library should do differently:** a library shouldn't configure the
root logger when imported — Python's logging HOWTO advises libraries to add no
handlers other than `NullHandler`. Removing these calls changes the logging
output of every downstream script that doesn't configure logging itself, so it
needs a coordinated change with gnomad_qc rather than a one-line fix.

---

## Per-stratum aggregates: loop-restrict-per-stratum, don't one-pass-inline-all

**`[Hail]`**

For per-population (or any per-stratum) counts / EM over a large keyed
dataset, prefer **restrict the cohort to each stratum, aggregate that
stratum on its own, then assemble the strata** over computing every
stratum in a single pass. Two payoffs:

1. **Codegen-safe** — sidesteps the `ClassTooLargeException` cliff (one
   heavy expression per stratum inlined into one row op; see
   [Hail idioms](hail-idioms.md#hail-idioms)).
2. **Cheaper shuffle** — each pass moves only that stratum's (smaller)
   data. Restricting a per-variant sample-set encoding to one group and
   re-indexing to that group's dense sample space shrinks the stored sets
   ~by the group's cohort fraction, so the co-partition / join shuffle is
   a fraction of the full-cohort size.

Derive the `"all"` stratum as the **element-wise sum over strata**, valid when the strata partition the cohort (every sample has a
group label); don't run a separate full-cohort pass for it.

**Tradeoff:** K passes instead of 1. But each pass is smaller *and it
actually completes at K≥3*, where the one-pass form both blows codegen and
shuffles the full-cohort sets. For a handful of strata this is a clear win;
for hundreds of tiny strata, reconsider (per-pass fixed overhead adds up).

---

[← back to the index](README.md)
