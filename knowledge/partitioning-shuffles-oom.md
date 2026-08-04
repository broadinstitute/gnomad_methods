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
| `hl.read_table(path, _n_partitions=n)` | sets partitioning as the table is read | no shuffle |

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

---

## Shuffle failures at scale are usually fleet-shaped, not code-shaped — shard the job

**`[gnomAD data]` · anecdotal, v4 exomes scale**

A `group_by` over the full genome that dies at scale is more often
losing shuffle blocks to preemption than hitting a real memory or skew
problem. Signature, from the driver log: `FetchFailed` /
`ExecutorLostFailure`, with large counts of `Container from a bad node`,
ending in a stage-retry abort — and no `OutOfMemoryError` anywhere.

What worked on the v4 variant-pair-list step, after the genome-wide pass
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

**If you can't avoid a big shuffle, isolate it.** Long-standing advice
that predates the run above and still holds: run the shuffling step on
non-preemptible workers only (they cost more, so keep the core count
down — less is more), and do large shuffles **one at a time**, writing
an intermediate between them, rather than chaining several into one
job. Cheap machines can come back for the non-shuffling stages.

---

## How many partitions?

**anecdotal · gathered from production runs, ~2023**

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
  tables was one place this bit.
- Aim for roughly constant *work per partition*, not a constant
  partition count.

---

[← back to the index](README.md)
