# Partitioning, shuffles, and OOM

## What an OOM actually looks like

**`[Hail]` · from a captured local-backend log, Hail 0.2.137**

Heap exhaustion in Hail is usually reported *inside the lowering or
encoding machinery*, not at the line of your code that caused it. A real
one, from a local-backend job that collected too much:

```
LoweringPipeline: ERROR: error while applying lowering 'EvalRelationalLets'
java.lang.OutOfMemoryError: Java heap space
    at java.util.Arrays.copyOf(Arrays.java:3745)
    at java.io.ByteArrayOutputStream.grow(ByteArrayOutputStream.java:120)
    at is.hail.io.StreamBlockOutputBuffer.writeBlock(OutputBuffers.scala:286)
    at is.hail.io.BlockingOutputBuffer.writeDouble(OutputBuffers.scala:229)
    at __C16collect_distributed_array_table_collect.__m87ENCODE_SFloat64$_TO_o_float64(...)
    at __C16collect_distributed_array_table_collect.__m84ENCODE_SBaseStructPointer_TO_r_struct_of_...
    at is.hail.backend.BackendUtils.$anonfun$runCDA$2(BackendUtils.scala:95)
    at is.hail.backend.local.LocalBackend$$anon$1.$anonfun$mapCollectPartitions$3(LocalBackend.scala:87)
```

**How to read it.** The generated class name is the diagnosis. Here
`collect_distributed_array_table_collect` and the `ENCODE_*` frames say
this died **encoding results to send back**, i.e. the answer being
collected was too big — not the computation. The long
`ENCODE_SBaseStructPointer_TO_r_struct_of_o_binaryANDo_binaryAND…` frame
spells out the struct being encoded, which usually tells you *which*
table it was.

**Read the generated frame first:**
- `…_table_collect` / `ENCODE_*` → too much coming back to the driver;
  see the driver-memory entry in
  [Dataproc / operational](dataproc-operational.md#dataproc--operational).
- `…native_writer` + `ClassTooLargeException` → codegen size, not
  memory; see [Hail idioms](hail-idioms.md#hail-idioms).
- Frames inside `RegionPool` / `Region` allocation → genuine per-partition
  data volume; shrink partitions or the row.

**On a cluster it looks different.** An executor that runs out of heap is
usually killed before it can print this — you see the executor loss
instead (`ExecutorLostFailure`, container killed), and the OOM itself is
only in that executor's log, not the driver log. That's why "no
`OutOfMemoryError` anywhere in the driver log" does **not** mean memory
was fine; see the fleet-vs-code entry below.

## What actually shuffles — and the `group_by` that doesn't

**`[Hail]` · Verified on Hail 0.2.134**

Operations that reorder or summarize across rows shuffle. A
non-exhaustive list:

- **Changing row order / key** — `key_by` / `key_rows_by`, liftover,
  `import_vcf` on an unordered VCF.
- **Joins** — foreign-key joins especially
  (`ht.annotate(x=ht2[ht.not_a_key])`); see the `ht[key_expr]` entry in
  [Hail idioms](hail-idioms.md#hail-idioms).
- **`repartition()`** — avoidable with `shuffle=False`.
- **Anything summarizing across rows** — `group_by` / `group_rows_by`,
  `annotate_cols` with an aggregation, `pca`, `pc_relate`.

The pair that looks identical and isn't:

```python
ht.aggregate(hl.agg.group_by(ht.gene, hl.agg.sum(ht.AC)))  # no shuffle -> dict
ht.group_by(ht.gene).aggregate(sum_ac=hl.agg.sum(ht.AC))   # shuffles   -> Table
```

The second re-keys the table by `gene` (its IR carries a `TableKeyBy`),
which is a shuffle. The first aggregates in a single pass and hands back
a plain `dict`.

**The catch the "avoid shuffles" advice usually omits:** that dict is
built in **driver memory**. It's the right tool for a few hundred
groups and the wrong one for a few million — at which point you want
the shuffle. See the driver-memory entry in
[Dataproc / operational](dataproc-operational.md#dataproc--operational).

---

## repartition vs naive_coalesce vs setting partitions at read time

**`[Hail]` · Verified on Hail 0.2.134**

| call | what it does | cost |
|---|---|---|
| `ht.repartition(n)` | full shuffle, equal-sized partitions (`shuffle=True` is the default) | a shuffle, at the next action |
| `ht.repartition(n, shuffle=False)` | merges existing partitions; `n` cannot exceed the current count | no shuffle |
| `ht.naive_coalesce(n)` | merges *adjacent* partitions, no rebalancing — can leave you badly skewed | no shuffle |
| `hl.read_table(path, _n_partitions=n)` | sets partitioning as the table is read | no shuffle on Spark; an extra stage on Query-on-Batch (see [hidden APIs](hidden-apis.md#hidden-and-undocumented-apis-we-rely-on)) |

**All of these are lazy** — the call returns instantly and the shuffle
happens at the next action, which is why a repartition usually gets
blamed on whatever line executed after it. (`n_partitions()` is itself
not free on a repartitioned plan.)

**Rule of thumb (anecdotal, from repeated production runs): avoid
repartitioning unless the dataset is small.** A repartition on a large
table buys balanced partitions at the price of a full shuffle — usually
the most fragile stage in the job. Prefer, in order: set the partition
count **at read time** (`_n_partitions`) or **at write time**; use
`naive_coalesce` when you only need *fewer* partitions after heavy
filtering; reach for `repartition` last.

**`naive_coalesce` also narrows the work before it.** It doesn't shuffle, so
it runs in the same stage as everything since the last materialization: that
upstream work runs on the coalesced number of tasks. Locally (0.2.134), an
expensive per-row annotation followed by `naive_coalesce(2)` ran as one 2-task
stage instead of 10 tasks. To shrink the partition count after a heavy filter
without also shrinking the parallelism of the expensive work, checkpoint the
filtered table first and coalesce after (or set the count at read time).

---

## A shuffle's output partition count is inherited from its input

**`[Hail]` · Verified on Hail 0.2.134**

This is the mechanism behind the rule above. `group_by(...).aggregate(...)`
does not choose a partition count — it carries through whatever the table
feeding it had. So a narrow input produces a narrow output no matter how many
rows come out.

Measured on one derive:

| source read at | rows out | partitions out | rows/partition |
|---|---|---|---|
| 99 (as written) | 973,085,319 | 99 | ~9.8 M |
| `_n_partitions=5000` | 973,085,319 | 4,493 | ~217 K |

Same content both times — only the layout changed.

**Why `repartition` afterwards is the wrong fix.** It pays for the problem
twice: the aggregate still runs on the narrow count, and then a *second* full
shuffle rebalances its output. The narrow stage is the slow one, and
repartitioning cannot make it any wider retroactively. On Spark, widening the
read costs nothing (`_n_partitions` re-splits off the table's own index) and the
aggregate never narrows in the first place. On Query-on-Batch the wider read
isn't free (see
[hidden APIs](hidden-apis.md#hidden-and-undocumented-apis-we-rely-on)).

**Corollary for the input you don't control.** The 99 above was not a choice
anyone made in this pipeline — it was the partitioning of an upstream published
table, inherited silently. When a shuffle stage is slower than its row count
justifies, check the partition count of what feeds it before tuning anything
downstream.

---

## Gene- or region-scale runs: `filter_intervals` keeps the table's partitions — re-read with `_intervals`

**`[Hail]` · Verified on Hail 0.2.134**

`hl.filter_intervals` (and `read_table(..., _intervals=..., _filter_intervals=True)`)
keeps the table's own partitioning and drops the partitions that don't overlap
— it never splits them. A few genes on a table written with large partitions
leave a handful of partitions, and the run happens on that many tasks however
big the cluster is: 5 genes gave 26 partitions, and the run was killed after
4 hours without the autoscaler adding workers.

Reading with `_intervals=` and *without* `_filter_intervals` makes the
intervals the partitions — one per interval, no shuffle — and still reads only
the rows inside them. So compute finer intervals over the rows you want and
re-read:

```python
genes_ht = hl.filter_intervals(hl.read_table(path), gene_intervals)
intervals = genes_ht._calculate_new_partitions(n_partitions)
ht = hl.filter_intervals(hl.read_table(path, _intervals=intervals), gene_intervals)
```

`_calculate_new_partitions` returns contiguous intervals over the filtered
table's key range, so they span the gaps between genes; the second
`filter_intervals` keeps only the gene rows. (A local test went from 3
partitions to 8 with the same rows.) `_intervals` can't be combined with
`_n_partitions`. gnomad_qc's `get_gnomad_v4_vds(chrom=..., n_partitions=...)`
uses this same pattern.

**The same mechanism co-partitions tables before a chain of joins.** Tables
keyed the same way but written separately have different partition bounds;
read them all with intervals computed once from the base table
(`base_ht._calculate_new_partitions(n)`), as gnomad_qc's
`create_release_sites_ht.py` does; gnomad_methods wraps this as
`gnomad.utils.file_utils.repartition_for_join`. When the base table is sparse —
a few genes — its intervals cover everything between the genes, so a single task
can end up reading a large stretch of the other table: filter the other table
to the genes first.

---

## Shuffle failures at scale are usually fleet-shaped, not code-shaped — shard the job

**`[gnomAD data]` · anecdotal, v4 exomes scale**

A `group_by` over the full genome that dies at scale is more often
losing shuffle blocks to preemption than hitting a real memory or skew
problem. Signature, from the driver log: `FetchFailed` /
`ExecutorLostFailure`, with large counts of `Container from a bad node`,
ending in a stage-retry abort — and no `OutOfMemoryError` anywhere.

What worked on a genome-wide v4 job with a large shuffle, after the full pass
failed repeatedly: **shard the job by chromosome** and run each shard on
its own smaller cluster with fewer preemptible workers. Tuning executor
memory did not help, because memory was never the problem.

**Distinguish before you tune** (see
[Dataproc / operational](dataproc-operational.md#dataproc--operational) for how to pull the
log):
- `FetchFailed` + `ExecutorLostFailure` + no OOM → lost shuffle blocks;
  shard the work, or use non-preemptible workers.
- One task at `(N-1)/N` for hours → data skew; shard the skewed key.
- Disk-full during shuffle spill → worker boot-disk sizing; see the
  `ht[key_expr]`-is-a-join entry in [Hail idioms](hail-idioms.md#hail-idioms).

**Hail's own shuffler is the other lever.** Wrapping the shuffling write in
`hl._with_flags(use_new_shuffle="1")` moves the blocks off Spark's
executor-local disk, which is what preemption takes with it. It is part of the recipe, not a
substitute for the fleet change — see
[Hidden and undocumented APIs](hidden-apis.md#backend-feature-flags-hl_set_flags--hl_get_flags--hl_with_flags).

**If you can't avoid a big shuffle, isolate it.** Long-standing advice
that predates the run above and still holds: run the shuffling step on
non-preemptible workers only (they cost more, so keep the core count
down — less is more), and do large shuffles **one at a time**, writing
an intermediate between them, rather than chaining several into one
job. Cheap machines can come back for the non-shuffling stages.

---

## How many partitions?

**`[Hail]` · anecdotal · gathered from production runs, ~2023**

There's no formula, but the failure modes are known:

- **Too many is not free.** A dataset written with ~50k partitions was
  slow to *read*; ~10k was materially better. Per-partition overhead is
  real.
- **Set it at read time to suit the pipeline** rather than
  repartitioning mid-job (see the entry above).
- **A join that grows each row wants more partitions than its inputs
  had.** Joining several large tables multiplies bytes per row while the
  partition count stays fixed, so partitions that were well-sized going
  in are oversized coming out. A release-scale join of six or seven
  tables was one place this bit. `repartition_for_join` builds in that
  headroom: by default it returns 1.1× the base table's partition count.
- Aim for roughly constant *work per partition*, not a constant
  partition count.


---

## Check partition counts and key ranges from a table's metadata, without Hail

**`[Hail]` · Verified on Tables written by 0.2.134 and 0.2.139**

A Table directory's `metadata.json.gz` holds the full schema (`table_type`) and
the row count of every partition (`components.partition_counts.counts`).
`rows/metadata.json.gz` holds each partition's key range (`_jRangeBounds`) and
part file (`_partFiles`). Reading them takes seconds and no cluster — enough to
tell whether a field exists, how many partitions will feed a shuffle, whether
they're skewed, and which partitions overlap your genes. (0.2.139 stopped
writing the per-partition index metadata, `index/*.idx/metadata.json.gz`, but
both of these files are still there.)

```bash
gcloud storage cat gs://<bucket>/<table>.ht/metadata.json.gz | gunzip | python -m json.tool | less
```

The same file records the Hail version that wrote the table; see
[Hail versions](hail-versions.md#hail-versions). MatrixTable directories are
laid out differently.

---

[← back to the index](README.md)
