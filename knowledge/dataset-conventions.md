# Release conventions and data policy

Why the released data looks the way it does. These are project decisions
rather than Hail or biology facts — the kind of thing that's obvious to
whoever was in the room and invisible to everyone else.

Entries marked **reported · project meeting** come from the project's own
meeting record rather than from something re-derived here. They're
accurate as statements of what was decided; check the current release
notes before treating any of them as the present state of the data.

## Annotation versions are pinned per release, and deliberately not bumped

**`[gnomAD data]` · reported · project meeting**

A release's VEP and GENCODE versions are fixed for the life of that
release. v4 and v5 sit on **VEP 105 / GENCODE 39** — partly because
downstream products are coupled to it, notably pext, which depends on
GTEx being on the same GENCODE build.

Staying put has a known cost: a problem in `RNU4ATAC` (a coordinate shift
relative to ClinVar, fixed upstream in VEP 114) persists in the pinned
version. The remedy was to patch that one gene rather than move the whole
release, and that patch has shipped — the browser now annotates
`RNU4ATAC` with **VEP 115 / GENCODE v49** while everything else stays on
the pinned version, and flags the gene page to say so. So "which
annotation version is this?" can have a per-gene answer.

Moving isn't cheap either. VEP 115 more than doubles the transcript count
(protein-coding and lncRNA alike) and takes roughly 3.5 hours where 105
takes 45 minutes — and every transcript-level annotation downstream
shifts with it.

**What this means for you:** don't assume gnomAD's consequences match
whatever VEP you run locally. When your annotation disagrees with the
browser's, compare versions before you debug anything else. If you need a
newer annotation for a specific gene, that's a targeted re-run, not a
reason to expect the release to move.

---

## LOEUF cutoffs do not port across gnomAD versions

**`[gnomAD data]` · reported · project meeting**

The constraint distribution shifts as sample size grows — as N rises,
observed/expected converges toward a normal distribution centered near 1,
so a fixed threshold means something different in each release. v2's
widely-used LOEUF < 0.35 corresponds roughly to pLI 0.9; the equivalent
stringency in v4 is nearer **0.45**.

Reusing a v2 threshold on v4 data therefore silently changes what your
gene list *is*, with no error and a perfectly plausible-looking result.

The project's own position is stronger than "pick the right cutoff":
**don't encourage binary cutoffs on continuous metrics at all.** Use the
metric as a continuous ranking where you can, and if you must threshold,
state the version the threshold was calibrated on.

**Related:** derived per-gene metrics can move across releases for
entirely non-biological reasons. Per-gene mutation rates differed by up
to ~10× between v2 and v4 for some genes, traced to coverage-cutoff
changes in which bases were included in the model — not to any change in
mutation. Re-derive, don't carry over.

---

## Which subsets exist, and why others don't

**`[gnomAD data]` · reported · project meeting**

v4 ships **UKB** and **non-UKB**. Two categories of subset people ask
for do not exist, both deliberately:

- **non-TOPMed** was planned and then dropped. gnomAD holds only a small
  fraction of TOPMed (on the order of 30k of ~160k samples at the time),
  it was expensive to construct without the metadata, and user demand
  was low. It could return in a point release if demand appears.
- **Disease-specific subsets** (the v2-era non-neuro being the one still
  asked about) ended after v2, on two grounds: no single disease is
  over-enriched in the callset, and the label was inaccurate about what
  the subset actually contained.

**What this means for you:** if your analysis needs "gnomAD minus cohort
X," check whether the subset exists before designing around it. It
usually doesn't, and the frequencies you want may have to come from a
subset that does exist plus an explicit caveat.

---

## Age data is binned because the extremes are protected

**`[gnomAD data]` · reported · project meeting**

Released age information is bucketed: **under 20**, five-year bins from
20 to 80, and **over 80**. Ages below 18 and above 80 are protected
information, so finer resolution at the tails isn't available at any
access level. It isn't a display choice you can work around by
requesting the underlying data.

---

## "Remaining individuals," not "other" or "unassigned"

**`[gnomAD data]` · reported · project meeting**

The label for samples not assigned to a genetic ancestry group is
**"remaining individuals."** It was chosen over "other" and over
"unassigned" on genetic-counselor feedback that `un-` and `non-`
constructions should be avoided when describing people.

Worth knowing when you're writing anything user-facing off gnomAD data,
and when you're matching labels programmatically — the string is
`remaining` in the data.

---

## Superseded releases stay online, unsupported

**`[gnomAD data]` · reported · project meeting**

Old versions aren't taken down — v2 remains available, partly because
its GRCh37 coordinates still have real users. But the browser defaults
to the newest release, superseded datasets carry a "legacy dataset"
banner, and questions about them are not answered.

**What this means for you:** building new work on a superseded release
means no support if something looks wrong, and the newer release may
have fixed exactly the thing you're about to report. Check whether your
question is version-specific before assuming it's a bug.

---

## Joint exome + genome FAF isn't defined for every variant

**`[gnomAD data]` · reported · project meeting**

Where a variant's exome and genome filtering allele frequencies differ
enough that combining them would be misleading, no combined value is
produced. Treat a missing joint FAF as "these two callsets disagree
here," not as "no data" — and go look at the two separately.

---

## The large Hail tables are requester-pays

**`[gnomAD data]` · Verified against Hail 0.2.134 API**

Reading the full Hail Tables/VDS means paying for the egress yourself,
which needs a billing project configured or the read fails.

```python
hl.init(gcs_requester_pays_configuration="my-project")
# or, to scope it to particular buckets:
hl.init(gcs_requester_pays_configuration=("my-project", ["bucket-a", "bucket-b"]))
```

On Dataproc, allow the buckets at cluster start:

```bash
hailctl dataproc start CLUSTER --requester-pays-allow-buckets BUCKET
# --requester-pays-allow-all also exists; prefer naming buckets
```

Two things that make this bite harder than it should. Costs are
**region-dependent** — reads from a distant region are markedly more
expensive than in-region ones, which is the same trap as the
compute/data region mismatch in
[Dataproc / operational](dataproc-operational.md#dataproc--operational).
And some third-party tools give you **no way to pass a billing project
at all**, so a requester-pays bucket is simply unreadable from them;
you'll need to copy what you need out first.

---

[← back to the index](README.md)
