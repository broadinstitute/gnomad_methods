# VDS

## Table vs MatrixTable vs VDS — what each one is for

**`[Hail]`**

- **Table** — rows and (optionally) a key. No sample dimension. Any
  per-sample data has to be exploded into rows or nested into an array
  field.
- **MatrixTable** — a rows × columns grid with an **entry** per
  (row, column). Dense: storage scales with `n_variants × n_samples`,
  including every hom-ref call.
- **VDS** — the sparse representation: two MatrixTables, a
  `reference_data` MT of **hom-ref blocks** (one entry covers a run of
  loci) and a `variant_data` MT holding only non-ref calls. Storage
  scales with the number of *calls*, not the grid.

The v4 exomes release is a VDS for exactly this reason — a dense MT of
that grid is not a thing you can hold. `hl.vds.to_dense_mt` expands it
back to the grid, so **always narrow first** (filter to your variant set
or intervals, and to your samples) and densify last. See the
[Hail VDS docs](https://hail.is/docs/0.2/vds/index.html) for the API.

**Rule of thumb:** if your analysis is per-variant, get to a Table as
early as you can. Stay in VDS/MT form only for the part that genuinely
needs per-sample entries.

---

## Column filters must be applied to BOTH reference_data and variant_data

**`[Hail]` · Verified on Hail 0.2.134**

A VDS's two MatrixTables must carry the **same columns in the same
order**. `to_dense_mt` pairs them **positionally**
(`hl.scan._densify(hl.len(var_cols), ref_entries)`), never by sample ID.
Filter columns on one and not the other and the densify **succeeds
silently, with every sample's reference data shifted onto a different
sample**.

Demonstrated with a per-sample tracer (sample *i* carries `DP = 10i` in
its reference blocks), dropping sample `0` from `variant_data` only:

| sample | expected DP | got |
|---|---|---|
| 2 | 20 | **10** |
| 3 | 30 | **20** |

No error, no warning — hom-ref depth, GQ, and every other
reference-block field lands on the wrong sample.

**Rule:** never filter one side alone. Either use
`hl.vds.filter_samples(vds, ht, keep=…)` (it filters both), or apply the
*same* predicate to both and rebuild:

```python
vmt = vds.variant_data.filter_cols(pred)
rmt = vds.reference_data.filter_cols(pred)
vds = hl.vds.VariantDataset(rmt, vmt)
```

**`validate()` is not a reliable backstop.** `vds.validate()` does check
that the column keys match — but on a VDS in the newer `LEN`
reference-block representation it crashes first with
`AttributeError: … no field … 'END'` (`variant_dataset.py:344` still
aggregates on `rd.END`), so on those datasets the check can't run at
all. Verify it yourself: `rd.count_cols() == vd.count_cols()`.

**Related:** `hl.vds.filter_samples(..., remove_dead_alleles=True)`
raises `AttributeError: … 'LA'` on a **split** VDS — recomputing local
alleles needs the `LA` field, which splitting removes.

---

## A VDS without `ref_block_max_length` gets slow interval filters — and Hail only warns

**`[Hail]` · Verified on Hail 0.2.134**

`hl.vds.filter_intervals` can only skip reference blocks that can't
reach the target interval if it knows the longest block in the dataset.
That bound lives in a global on `reference_data`
(`ref_block_max_length`), written by newer Hail. Read a VDS that predates
it and you get a **warning, not an error**:

> You are reading a VDS written with an older version of Hail. Hail now
> supports much faster interval filters on VDS, but you'll need to run
> either `hl.vds.truncate_reference_blocks(vds, ...)` and write a copy
> or patch the existing VDS in place with
> `hl.vds.store_ref_block_max_length(vds_path)`.

Nothing breaks; interval filtering just can't prune, so it reads far
more reference data than it needs. Hail emits enough warnings on a
normal job that this one is easy to scroll past — the symptom people
notice is "interval-filtered reads of this VDS are much slower than
they should be."

**Fix:** patch the dataset once with
`hl.vds.store_ref_block_max_length(vds_path)`. Check first with
`"ref_block_max_length" in vds.reference_data.globals` — cheaper than
re-reading the warning log.

---

## Split a VDS by columns, not `hl.vds.filter_samples`, if you'll rejoin/densify

**`[gnomAD data]` · `[internal only]` · Applies to:** v4 (the release VDS is split into UKB / non-UKB
cohorts for size; any analysis that subsets VDS samples with intent to
recombine).

`hl.vds.filter_samples()` doesn't only subset columns — after removing samples
it drops `variant_data` rows with no defined entry in the retained set (i.e.
sites with no alt-allele carrier in the subset), regardless of
`remove_dead_alleles`, in every Hail version through ≥0.2.128. That's fine for a
standalone subset. But if you split a VDS into cohorts intending to rejoin or
densify across the full variant set, the dropped sites mean the subset's
reference blocks (hom-ref calls) at those loci can't be placed on densify →
silently reduced AN for any variant private to the *other* subset. This bit the
v4.0 release freq; corrected in v4.1 (`fix_freq_an.py`, gnomad_production#1366).

Use gnomad_qc's `split_vds` pattern — filter columns on **both** datasets
(which keeps every variant site, unlike `filter_samples`) and reconstruct
the VDS:

```python
# BAD: drops sites with no alt carrier in the subset → AN loss on rejoin
ukb = hl.vds.filter_samples(vds, ukb_ht, keep=True)

# GOOD: gnomad_qc/v4/annotations/generate_freq.py::split_vds
vmt = vds.variant_data.filter_cols(pred)      # keeps ALL variant rows/sites
rmt = vds.reference_data.filter_cols(pred)    # same predicate — see entry above
rmt = rmt.filter_rows(hl.agg.count() > 0)     # prune now-empty ref blocks (safe)
ukb = hl.vds.VariantDataset(rmt, vmt)
```

**Symptom:** AN lower than expected for variants private to one cohort after
rejoining split VDSs; freq that won't reconcile with the unsplit callset; no
error.

**Canary:** the subset's `variant_data` row count should equal the input's (all
sites retained). If it shrank, filter_samples-style site dropping happened.

**Two things that made this hard to spot** when it happened, both worth
knowing before you trust a subset's variant list:

- `remove_dead_alleles=False` was expected to preserve every locus. It
  doesn't govern that at all — the row dropping happens regardless, and
  the subset behaves more like a project-VCF of its own samples than a
  view onto the callset.
- Comparing subsets by variant count is muddier than it looks with
  **multiallelics**. After densify-and-split, a site can carry call-stats
  information in a subset that has no carrier of that allele, because
  depth data alone is enough to define the site. So "same variants on
  both sides" is not a clean check, and the two subsets plus their union
  need not agree on row count.

**Where to find more:** `split_vds` in `gnomad_qc/v4/annotations/generate_freq.py`;
`fix_freq_an.py`; gnomad_production#1366.

---

[← back to the index](README.md)
