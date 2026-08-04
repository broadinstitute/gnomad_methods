# gnomAD-specific

## `filters.length() == 0` means "PASS AND AC>0", not just "PASS"

**`[gnomAD data]` · Applies to:** v4 (and v2, v3 by analogy — `filters` is a set of
*reasons for failure*, not a status string).

gnomAD's final_filter HT stores QC status as `Set[String]`. Empty set =
passes everything. But **`"AC0"` is one of the filter reasons** —
variants with `adj_freq.AC == 0` in the release cohort get
`filters = {"AC0"}` alongside sets for AS_VQSR, InbreedingCoeff,
Monoallelic, etc.

See [`gnomad_qc/v4/variant_qc/final_filter.py`](https://github.com/broadinstitute/gnomad_qc/blob/main/gnomad_qc/v4/variant_qc/final_filter.py):
```python
filters = {
    ...
    "AC0": adj_freq_expr.AC == 0,
    ...
}
```

**Consequence:** `filters.length() == 0` drops both fail-QC AND
AC0-in-release variants together. For trio-side analyses (or any
analysis that legitimately wants variants carried only by non-release
samples), this silently loses a large fraction of pairs.

**Correct patterns:**
- **Pure release-side analysis** (variant needs a release AC>0 to be
  meaningful): `filters.length() == 0` is correct.
- **Trio / non-release-inclusive analysis** ("PASS *or* AC0-only"):
  ```python
  filters.difference(hl.set(["AC0"])).length() == 0
  ```

**Symptom:** trio-adjacent aggregations with totals materially less than
input rowcount, or `missing_af` / `missing_source` buckets that grow
much larger than expected.

---

## `only_filters.ht` (all-variants) exists

**`[gnomAD data]` · `[internal only]`**

Alongside the release-scoped `final_filter.ht` there is an
**all-variants** sibling — same directory, `…final_filter.all_variants.only_filters.ht`
— reachable through `gnomad_qc.v4.resources.variant_qc` rather than by
typing a bucket path.

Schema: `(locus, alleles) → filters: Set[String]`. Covers **every v4
exome variant**, not just AC>0-in-release. Not obvious this exists
unless you already knew — the "normal" `final_filter.ht` sibling only
has the release-scoped subset.

Use this when you need PASS status for variants absent from release-
scoped tables (sites, freq, vep). Combine with the `.difference({"AC0"})`
pattern above to build "permissive" annotation paths for trio / non-
release / rare-tail analyses.

Genomes counterpart at the analogous path.

---

## `only_filters.ht` and release `final_filter.ht` are NOT interchangeable for AC0 semantics

**`[gnomAD data]` · `[internal only]`**

Both HTs have the same `filters: Set[String]` schema. But their **`AC0`
tag is defined against different cohorts**:

- **Release `final_filter.ht`** (from `gnomad_qc.v4.resources.variant_qc.final_filter(data_type='exomes')`):
  `AC0` fires when `adj_freq[release_cohort].AC == 0`. So "AC0" means
  "no release-sample carrier."
- **`all_variants.only_filters.ht`**: `AC0` appears to be defined against
  a wider joint-genotyping cohort, so it effectively never fires on real
  variants. Empirically, applying its filter check to a 682k-pair trio
  substrate classified 100% of pair endpoints as `filters={}` — no
  `AC0`-tagged variants at all.

**Rule of thumb:**
- Use `only_filters.ht` when you want *global* QC status (does this
  variant pass v4 joint-cohort QC).
- Use release `final_filter.ht` when you want *release-cohort* QC status
  (does this variant pass QC **and** have AC>0 in release samples).
- The two answer different questions. Picking the wrong one is silent —
  your PASS check will match a different variant population than you
  intend.

---

## `freq` is an array of strata — index it via `freq_index_dict`, not a hardcoded number

**`[gnomAD data]` · Verified against the released v4.1 exomes public HT**

Not `freq["all"].AF`. `freq` is an array of structs
(`AC` / `AF` / `AN` / `homozygote_count`), one per sample grouping —
**329 of them** in v4.1 exomes. The supported way to reach a stratum is
the `freq_index_dict` global, a dict from grouping label to array index:

```python
from gnomad.resources.grch38.gnomad import public_release
ht = public_release("exomes").ht()

ht.freq[ht.freq_index_dict["adj"]]              # whole callset, high-quality genotypes
ht.freq[ht.freq_index_dict["afr_XX_adj"]].AC    # one gen-anc × sex stratum
```

Label keys compose as group / sex-group / subset-group /
gen-anc-group / gen-anc-sex-group / subset-gen-anc-sex-group
(`adj`, `raw`, `XX_adj`, `non_ukb-raw`, `afr_adj`, `ami_XX_adj`,
`non_ukb_mid_XX_adj`, …). The same pattern applies to the filtering
allele frequency array via `faf_index_dict`.

For the record, `freq_index_dict["adj"] == 0` and `["raw"] == 1`, and
`freq_meta[0] == {"group": "adj"}` — so the familiar `freq[0]` does mean
"whole callset, adj." **Use the dict anyway**: hardcoding an index that
happens to be right today silently returns a *different stratum* if the
ordering ever changes, and reads as a magic number in review.

**Related globals:** `freq_meta` (ordered list of the grouping structs,
parallel to `freq`), `freq_meta_sample_count`, and `frequency_README` —
an explanation of this whole scheme carried on the HT itself.

**Where to find more:** gnomAD's
[v4 HTs help page](https://gnomad.broadinstitute.org/help/v4-hts)
documents the released HT schemas and this access pattern (source:
`browser/help/topics/v4-hts.md` in the gnomad-browser repo).

---

## Sample releasability skews annotation coverage

**`[gnomAD data]` · `[internal only]` · Applies to:** v4 (multi-cohort release model; less relevant in v2).

gnomAD joint-calls variants across ALL samples (releasable +
non-releasable) but the release freq HT only counts releasable samples.
Consequence:
- Variants carried only by non-releasable samples → AC=0 in release
  freq → dropped from sites HT (by the AC>0 filter) → absent from
  downstream release-scoped annotations.
- BUT those variants ARE in the "all variants" final_filter HT above,
  and usually in vep_ht (which is annotate-once at joint-call time).

Whenever your analysis spans releasable and non-releasable samples
(trios where some members are non-release, external cohorts joined to
gnomAD, family sub-sets), you'll see this gap. It commonly manifests as
a large fraction of trio-truth pairs with null AC/AF because at least
one endpoint is a "trio-only variant" absent from release annotations.

**Rule of thumb:** for any analysis touching non-release samples, use
the permissive filter path (`only_filters.ht` + AC0-allowed check)
rather than the release-scoped sites HT.

---

## "Absent from gnomAD" often means "not covered," not "not present"

**`[gnomAD data]` · reported · project meeting**

A variant missing from gnomAD is routinely read as evidence that nobody
in 800k people carries it. Often it means nobody **sequenced** it there.
Exome capture kits differ in what they cover, and coverage varies across
exons as a result, so "not in gnomAD" mixes two very different
statements:

- the site was callable in the cohort and no alt allele was seen, and
- the site wasn't well covered, so no call could be made either way.

The distinction matters most exactly where people lean on it hardest —
arguing a candidate variant is rare or novel because gnomAD lacks it.

**Rule of thumb:** before treating absence as evidence, check the
coverage at that locus (and `AN`, which is the same statement in
frequency terms — see [Biology](biology.md#biology)). Absence at a
well-covered site is informative; absence at a poorly-covered one is not.

---

## Different gnomAD HTs may have different variant sets

**`[gnomAD data]`**

An annotation join graph like `variant_filter_ht → sites_ht →
annotated_phase_ht` may not have identical row coverage between steps:
- `variant_filter_ht` may include tags for variants absent from
  `sites_ht`
- `sites_ht` may drop AC0 variants that `only_filters.ht` retains
- `vep_ht` and `freq_ht` are typically joint-call-scoped (broader) but
  `sites_ht` is release-scoped (narrower)

A downstream function that reads a `v_ann` struct may find:
- `v_ann.source` defined (from `variant_filter_ht`)
- `v_ann.ac` null (variant absent from release `sites_ht`)
- `v_ann.vep` defined or null depending on `vep_ht` scope

**Always check field-level nullability, not just struct-level.**
`hl.is_defined(v_ann)` does NOT guarantee `hl.is_defined(v_ann.ac)`.

---

## Identifying UKB samples in v4: `s.startswith("UKB")` vs `project_meta.ukb_sample`

**`[gnomAD data]` · `[internal only]` · Applies to:** v4 exomes. Any analysis that subsets or removes the UKB
cohort.

The authoritative UKB flag is `meta().project_meta.ukb_sample` (defined as
`project == "UKBB"` in
`gnomad_qc/v4/sample_qc/create_sample_qc_metadata_ht.py`). The sample-ID prefix
`s.startswith("UKB")` is a lighter, meta-free proxy used elsewhere in gnomad_qc
(e.g. `identify_trios.py`). Across all 955,213 v4 exomes meta samples the two
agree on every sample except **2** UKBB-project *control* samples that carry
control IDs rather than a UKB prefix: `CHMI_CHMI3_Nex1` (CHM1 mole) and
`Coriell_NA12878_NA12878` (NA12878), both non-release. Critically, **zero**
non-UKB samples carry the "UKB" prefix, so the prefix never leaks a non-UKB
sample into a UKB subset — it only leaves those 2 controls on the non-UKB side.

**Tradeoff:** the prefix is dependency-free and safe for bulk cohort splits (no
non-UKB leakage). Use `project_meta.ukb_sample` when you must remove/retain
*all* UKBB-project data exactly — e.g. UKB data-removal / governance — so those
2 controls go with the rest of UKB.

**Where to find more:** `create_sample_qc_metadata_ht.py`
(`ukb_sample = project == "UKBB"`); prefix use in `identify_trios.py`.

---

## Large globals structs ride along in every downstream IR — `select_globals()` early

**`[Hail]` + `[gnomAD data]`**

A table's globals are carried into every expression built from it, so a
large globals struct bloats every downstream IR whether or not you touch
it — and past a certain size it can trip Hail's optimizer outright, at
execute time, before any rows are processed. **If you only need the
rows, drop the globals right after reading:**

```python
ht = some_table.ht().select_globals()
```

The gnomAD instance: the v4 sample-QC `meta()` HT carries a
deeply-nested `global_annotation_descriptions` struct. Carrying it
through a nontrivial `select → annotate → filter → checkpoint` chain
raises `FatalError: AssertionError: type mismatch … bad ir from forward
lets`, with the dump showing the giant `project_meta` struct.

```python
ht = meta(data_type="exomes").ht().select_globals()
```

**Symptom:** "bad ir from forward lets" / a type-mismatch assertion on
an otherwise simple pipeline, raised at execute time before any rows are
processed. Any richly-annotated HT can do this; the meta HT is just the
one we hit.

---

## v4 genotype classification: match the high-AB het→hom-alt correction, and get the operation order right

**`[gnomAD data]` · Applies to:** v4 (exomes + genomes). Any analysis that densifies the
release VDS and classifies per-sample genotypes (het / hom-var /
hom-ref) for counting, freq, or phasing — if you want to reproduce the
released `freq`.

GATK <4.1.4.1 mis-called some true hom-alts as high-allele-balance
het-refs. gnomAD v4's released `freq` is `ab_adjusted_freq`:
`correct_for_high_ab_hets` reclassifies such calls het→hom-var. If you
classify genotypes yourself and skip this, your het / hom-var counts
won't match release — silently; the totals still look sane. Reclassify a
call het-ref→hom-var when:

```python
GT.is_het_ref() & adj & (AD[1]/DP > 0.9) & ~is_het_non_ref \
    & ~fixed_homalt_model & (af > 0.01)
```
(`fixed_homalt_model` from `meta().project_meta`; `af` = release
`freq[0].AF`; cutoffs are the `ab_cutoff=0.9` param on `densify_and_prep_vds_for_freq`, AF threshold `af_threshold=0.01` on `correct_for_high_ab_hets`.)

Two traps:

1. **The correction applies to RAW too, not just adj.** The released
   `freq[1]` (raw) is the raw-group aggregate of the *adj-gated*
   `high_ab_het`, so a corrected call lands in both raw and adj hom-var.
   Reclassify in raw as well — but keep it **gated on adj-pass** (a
   non-adj high-AB het is not in the correction set, so it stays het in
   raw). The gnomad_qc docstring says "raw is not adjusted," but the code
   + released data do adjust it — verify against released `freq[1].AC` if
   unsure.

2. **Capture `is_het_non_ref()` BEFORE the split.**
   `sparse_split_multi` downcodes a het-non-ref (e.g. `1/2`) to `0/1` at
   each biallelic row — indistinguishable from a true het-ref afterward.
   Test `is_het_non_ref()` post-split and it's always False, so the
   correction wrongly fires on split multiallelic calls. Capture it from
   the local `LGT` pre-split and carry it through as a passthrough entry.

**Order of operations** (matches `generate_freq.densify_and_prep_vds_for_freq`):
```
capture _het_non_ref (LGT, pre-split) → sparse_split_multi → to_dense_mt
→ adjusted_sex_ploidy_expr(GT) → adj = get_adj_expr(GT, …) → high-AB correction
```
Each dependency is load-bearing:
- **sex-ploidy before adj** — adj and the correction must see the
  haploid-adjusted GT; an XY chrX het dropped to missing by ploidy
  adjustment then can't be "corrected" into a hom-var.
- **correction after adj** — it's adj-gated.
- **`_het_non_ref` before the split** — see trap 2.

Reorder any of these and the classification silently diverges from
release.

**Where to find more:** `correct_for_high_ab_hets` +
`densify_and_prep_vds_for_freq` in
`gnomad_qc/v4/annotations/generate_freq.py`.

---

[← back to the index](README.md)
