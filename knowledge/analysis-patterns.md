# Codebase and analysis patterns

## Row-set anchor determines what "could possibly appear"

**`[Hail]` + `[gnomAD data]`**

`hl.Table.filter()` removes rows; it can never add them. Similarly, an
"annotate" step's `.annotate(**...)` adds *columns* but not *rows* — the
row set is fully determined by the input HT's row set. Choice of the
starting HT (the "anchor") determines the outer bound on what downstream
code can see.

**Two concrete instances** where this trap fires in the variant-
cooccurrence pipeline:

1. **`assemble_sites_ht`** starts with `ht = vep_ht.select("vep")` —
   anchoring on `vep_ht`. If a variant is in the filter HT but missing
   from vep_ht, it's gone before any `.filter()` runs. Swapping the
   filter for something more permissive doesn't help.

2. **`annotate_phased_ht_with_sites`** takes a `phased_ht` and adds
   `v1_ann`/`v2_ann` per-endpoint via sites-HT lookup. The output's
   row set is `phased_ht`'s row set. If you swap in a more permissive
   sites HT but leave `phased_ht` pointing at a release-scoped table,
   trio-only pairs (absent from the release phased HT) still get no
   annotation row — they get null `v_ann` on the downstream left-join
   from consumers. The fix isn't in the sites HT, it's in what you pass
   as `phased_ht`. Anchor on the trio-comparison HT if you want trio-
   only pairs annotated.

**Symptom of getting this wrong:** you build a "permissive" annotation
path and downstream aggregations produce bit-identical output to the
release-scoped baseline. The permissive lookups happened, but the row
set never grew.

**Rule of thumb:** anchor on the widest-coverage table for the
semantics you want. Join everything else with nullability tolerance
downstream. When you have both a "release" and a "trio-inclusive"
consumer, you likely need two annotated HTs — one anchored on the
release pair space, one on the trio comparison. Document the anchor
choice in the function docstring so callers can predict what's missing.

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
a *catch-all* for pairs that didn't match any preceding `.when()`. It
does **not** mean "the AF value was null" — it means "the AF didn't fit
any of the buckets you defined."

Concrete failure mode from a session:

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

If the analysis populates the pair list from a dataset that includes
higher-AF variants (e.g., OE-candidate expansion up to AF 0.5), those
common variants land in `missing_af` — same bucket as genuinely null-AF
variants. The label is misleading. We saw a session where 438k of 682k
trio pairs (64%) landed in `missing_af` and were interpreted as "AC0
coverage gap" — actually they were common-AF variants overflowing Guo
2024's rare-focused bucket range.

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

**Symptom:** a large fraction of your pairs end up in the default
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

Derive the `"all"` stratum as the **element-wise sum over strata** (v2's
fold-sum), valid when the strata partition the cohort (every sample has a
group label); don't run a separate full-cohort pass for it.

**Tradeoff:** K passes instead of 1. But each pass is smaller *and it
actually completes at K≥3*, where the one-pass form both blows codegen and
shuffles the full-cohort sets. For a handful of strata this is a clear win;
for hundreds of tiny strata, reconsider (per-pass fixed overhead adds up).

**Where this lives:** the variant-cooccurrence per-population counts
(`compute_counts_by_pop`: `restrict_encoded_to_pops` per group → flat
count → assemble; `"all"` = sum over groups).

---

[← back to the index](README.md)
