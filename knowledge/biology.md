# Biology

Fewer of these came up in a specific recent session; contribute more as
you hit them.

## Multiallelic sites: split vs unsplit changes everything

**`[Hail]` + `[gnomAD data]`**

A single VCF row like `chr1:100 A→G,T` splits into two biallelic rows
`chr1:100 A→G` and `chr1:100 A→T`. Whether your HT is split determines:
- Whether `(locus, alleles)` is a unique key (split: yes; unsplit: no —
  same locus appears multiple times with a length-≥3 alleles array)
- Whether AC / AF are per-alt or aggregated
- Whether VEP consequences are per-alt or shared

Reading an unsplit HT and treating it as split (or vice versa) produces
plausible-looking wrong output. Check `hl.experimental.split_multi_hts`
or `hl.split_multi` provenance before joining. The v4 VDS is read split
by convention; v2 was mixed.

**Canary: an AF spike at exactly 0.5.** A pile-up of allele frequencies
landing on 50% is a classic signature of multiallelic mishandling — a QC
or splitting step that forces AF to 0.5 at multiallelic sites rather than
apportioning it (*reported · project meeting*). If a frequency histogram
has a spike at 0.5 that biology doesn't explain, check how multiallelics
were split before interpreting anything downstream.

---

## Reference genome coordinates: GRCh37 vs GRCh38

**`[Hail]` + `[gnomAD data]`**

gnomAD v2 is GRCh37; v3/v4 are GRCh38. Loci from different builds are
not comparable. Cross-referencing v2 and v4 requires liftOver (which
loses ~1–2% of variants and may flip alleles). Related: chromosome
naming (`1` vs `chr1`) and MT vs chrM. Hail's `hl.get_reference('GRCh38')`
uses `chrN` prefix.

**Hail defaults to GRCh37** (`[Hail]` · verified on 0.2.134). `hl.init()`
without `default_reference=` gives you GRCh37, so `hl.locus`,
`hl.parse_locus_interval`, and `hl.tlocus()` all build GRCh37 objects
unless told otherwise. Because GRCh38 in Hail uses `chr`-prefixed
contigs, the mismatch usually surfaces as a loud error (`invalid
interval expression: 'chr1:1000-2000'`, or a reference mismatch on a
join) — but it surfaces *at execute time*, which on Dataproc means
you've paid for the cluster to get there.

```python
hl.init(default_reference="GRCh38", tmp_dir="gs://…")                 # do this
hl.parse_locus_interval("chr1:1000-2000", reference_genome="GRCh38")  # or be explicit
```

---

## Reference false duplications make some medically relevant genes unreliable

**`[gnomAD data]` · Verified against the paper and the shipped browser flag**

Both GRCh37 and GRCh38 contain **false duplications** — regions
erroneously present twice in the reference. Reads that belong to one
real locus get split across the copies, so coverage and genotypes at the
real gene are wrong, and the effect is **reference-specific**: a gene
broken in GRCh38 can be fine in GRCh37 and vice versa. Wagner et al.
([*Nat Biotechnol* 2022](https://www.nature.com/articles/s41587-021-01158-1))
characterized 273 of the ~395 medically relevant genes that standard
benchmarks exclude for repetitiveness, and found that masking these false
duplications improved variant recall **from 8% to 100%** in affected
genes.

On GRCh38 the ones that matter for gnomAD are **CBS, KCNE1, and CRYAA**.
The browser flags exactly these three on GRCh38:

> Variant calls in this gene may be missing or unreliable due to false
> duplications in the GRCh38 reference.

**What this means for you:** absence or an odd frequency in one of these
genes is a reference artifact until proven otherwise. Because the defect
is build-specific, cross-checking the other build (v2, on GRCh37) is a
genuine diagnostic rather than a fallback. Long-read and assembly-based
callsets — the All of Us CMRG callset is the one gnomAD plans to
integrate — are the real fix.

---

## Sex chromosomes: X/Y ploidy, PAR regions

**`[gnomAD data]`**

- Males are hemizygous outside PAR on chrX (haploid) — AC/AN counting
  and AF calculation need special handling.
- Pseudoautosomal regions (PAR1, PAR2) are diploid in both sexes on
  chrX/chrY; the gnomAD freq HT typically treats PAR separately.
- `locus.in_autosome()` excludes chrX/Y/M (and chr sex-linked regions).
  Trio-analysis or PBT code that filters to autosomes drops these.

If you're doing anything sex-chromosome-inclusive, check that AC/AN
denominators respect ploidy and PAR boundaries.

**chrY ploidy can't be called from variants alone.** There is too little
variant data on Y to infer ploidy from calls, so it is derived from
reference-block coverage instead, normalized against an autosome
(gnomad_methods exposes the choice as
`impute_sex_ploidy(..., use_only_variants=...)`; the default path uses
reference blocks). X tolerates a variant-based approach in a way Y does
not (*reported · project meeting*) — worth knowing before you assume the
two chromosomes were inferred the same way.

---

## Trio / family / relatedness gotchas

**`[gnomAD data]` · `[internal only]`**

- **PBT ∩ release is a small subset.** Many trio members in gnomAD v4
  are non-releasable, so the intersection between PBT trios and release
  samples is much smaller than either alone (~2,318 of 730,947 release
  samples in v4). Any subtraction (release-minus-PBT) has this small-
  overlap consequence.
- **Trio-only variants exist.** Variants carried only by non-releasable
  trio members won't appear in release freq (AC=0). See the
  releasability item under gnomAD-specific.
- **Consanguinity / high inbreeding.** Some cohorts have elevated
  inbreeding coefficients. Filter thresholds tuned to outbred cohorts
  may not port cleanly.
- **De novo vs transmitted vs family-recurring.** Trio phasing (PBT /
  transmission-phase) works only for variants transmitted from a parent.
  De novos and mosaic sites don't get phased-by-transmission tags.
- **Relatedness inference depends on the cohort you run it on.**
  `pc_relate` run across a whole diverse callset behaved very differently
  from the same method run on a subset — at full scale it inferred
  individuals across diverse genetic ancestry groups to be related to one
  another, because the PC space it corrects against is not the same space
  (*reported · project meeting*). Relatedness results are not portable
  between a subset and the full callset; re-run rather than reuse, and be
  suspicious of a relatedness matrix that makes a whole ancestry group
  look like a family.

---

## AC=0 vs AF=0 vs AF=null are three distinct cases

**`[gnomAD data]`**

- **AC=0** (defined, zero) — variant was called at the locus, but no
  alt-allele carriers in the current cohort. AF = 0/AN is defined and
  zero.
- **AF=0** (defined, zero) — usually implies AC=0. Sometimes reflects
  rounding at very low AF (unlikely with float64).
- **AF=null** (missing) — usually because AN=0 (nobody callable at that
  locus in this cohort) or because the variant is absent from the HT
  entirely. `AF = AC / AN` is undefined at AN=0.

Downstream code often conflates these. Filters like `AF > 0` drop AC=0
AND missing-AF variants together; `hl.is_defined(AF)` distinguishes the
two-null case from the zero case.

**And AF moves when AN moves.** `AF = AC / AN`, so anything that guts AN
at a site — a coverage dropout, a quality problem affecting the
denominator more than the numerator — inflates AF without a single extra
carrier. A site whose AN has collapsed relative to its neighbours can
show a dramatically higher AF that is an artifact of callability, not
frequency (*reported · project meeting*).

**Rule of thumb:** read AN next to AF, always. Treat low-AN sites as
unstable rather than rare-and-interesting, and consider an AN floor
(e.g. a percentage of the cohort maximum) before ranking anything by AF.

---

## QUAL is frequency-coupled, not a site-quality score

**`[gnomAD data]` · reported · project meeting**

QUAL reads like "how good is this site," and it isn't. It is closer to
the probability that **at least one sample in the callset carries at
least one alt allele** — so it scales with how common the variant is, not
with how trustworthy the call is. A segmental duplication throwing
systematic error across thousands of samples can carry an extremely high
QUAL precisely *because* the error is everywhere.

**Use the QC filters and the allele-specific metrics instead** — the
`filters` set encodes the actual pass/fail decision, and the AS_* metrics
are the per-allele evidence. A high QUAL is not a reason to trust a site
that the filters flagged.

---

## LOFTEE evaluates one variant at a time

**`[gnomAD data]` · reported · project meeting**

LOFTEE's pLoF calls and flags are computed **per variant, in isolation**.
It has no view of a second variant that might change the interpretation:
a nearby frameshift restoring the reading frame, two pLoFs in phase on
the same haplotype versus on opposite haplotypes, or a rescuing variant
downstream. Doing that jointly is hard computationally and is outside
what the tool attempts.

**What this means for you:** if your analysis is co-occurrence-shaped —
compound heterozygosity, phasing, haplotype-level consequence — the
per-variant LoF annotation is an input, not an answer. Two `HC` pLoFs in
the same gene is not the same finding as two pLoFs *in trans*, and
LOFTEE cannot tell you which you have.

(One flag that is worth understanding for the same reason: pLoF variants
in the last exon are flagged as predicted to escape nonsense-mediated
decay — an escape prediction, again made one variant at a time.)

---

## PCR amplifies shorter fragments — so PCR+ samples look *less* expanded

**`[gnomAD data]` · reported · project meeting**

At repeat loci, the naive expectation is that PCR+ samples will show
longer repeat expansions than PCR-free ones. The observed direction is
the opposite: PCR preferentially amplifies **smaller** fragments, so
PCR+ samples appear *less* expanded.

Relevant any time you compare repeat genotypes across samples prepared
differently — the library prep is a confounder with a counterintuitive
sign, and "the PCR+ cohort has shorter repeats" is a technical artifact
before it is a biological finding.

---

## Synonymous O/E ≈ 1 is the constraint sanity check

**`[gnomAD data]` · reported · project meeting**

In a constraint analysis of essential genes the expected picture is:
missense depleted, pLoF close to absent, and **synonymous
observed/expected sitting near 1**. Selection doesn't act on synonymous
variation the way it acts on the other two, so that ratio is your
control.

If synonymous O/E drifts away from 1, suspect the pipeline before the
biology — coverage handling, the mutational model, or which bases were
included are the usual causes. It's the cheapest available check that a
constraint calculation is behaving.

---

[← back to the index](README.md)
