# Hail idioms

## `.filter()` / `.select()` / `.annotate()` produce a new source object — don't mix with pre-op references

**`[Hail]`**

Hail tracks the source object of every expression. After
`sub = ht.filter(cond)` (or `.select` / `.annotate` / any re-binder),
expressions built off `ht` and expressions built off `sub` are considered
to come from *different sources*; mixing them raises
`Cannot combine expressions from different source objects`. The error
names the operation, not the two sources, so in a long
`.select(a=…, b=…, c=…, …)` where one of ten fields still points at
`ht` it takes a minute to find the stale reference.

```python
sub = ht.filter(cond)
sub.select(x=sub.foo)     # GOOD — every reference goes through the new binding
sub.select(x=ht.foo)      # BAD — ht.foo is from the pre-filter source
```

**Rule of thumb:** whenever you assign the result of `.filter` /
`.select` / `.annotate` to a new variable, refer to fields *only*
through that new variable. Two-liner filter-then-project is the common
site because it's tempting to reach back into the outer name.

**Same-name reassignment (`ht = ht.annotate(...)`) has the same failure
mode.** If you bind an expression to `ht` at the top of a function, then
later do `ht = ht.annotate(...)`, the earlier expression is now bound to
a stale source object. The bug is easy to miss because it looks like
"just modifying `ht` in place":

```python
# BAD -- in_training_genes is bound to ht #1, but ht is rebound below.
in_training_genes = ht.is_phaplo_gene | ht.is_ddg2p_nonlof_gene
if af_threshold_changed:
    ht = ht.annotate(is_common=recomputed)      # ht is now object #2
# ...200 lines later...
ht = ht.annotate(
    label=hl.if_else(in_training_genes & ht.is_blb_clnsig, 0, ...),
    # ^ crashes: in_training_genes from #1, ht.is_blb_clnsig from #2
)
```

Two fixes:

1. **Re-bind the expression right before you use it.** Cheap if it's
   built from a couple of fields:

   ```python
   ht = ht.annotate(is_common=recomputed)
   in_training_genes = ht.is_phaplo_gene | ht.is_ddg2p_nonlof_gene  # rebind
   ht = ht.annotate(label=hl.if_else(in_training_genes & ht.is_blb_clnsig, 0, ...))
   ```

2. **Wrap the expression as an inline closure that takes the current
   `ht` and resolves at call time.** Useful when the expression is
   reused at several points across multiple `ht.annotate(...)` calls,
   so rebinding at every use-site is noisy:

   ```python
   def _tg(h):
       return h.is_phaplo_gene | h.is_ddg2p_nonlof_gene

   ht = ht.annotate(is_common=recomputed)              # ht rebound
   ht = ht.annotate(label=... _tg(ht) & ht.is_blb_clnsig ...)   # OK
   ht = ht.annotate(is_reference=... _tg(ht) ...)              # OK, still resolves fresh
   ```

   The closure defers the field lookup until call time, so it always
   picks up the current binding of `ht` at the use site.

**Symptom:** `Cannot combine expressions from different source objects`
inside a function that has one or more `ht = ht.annotate(...)` /
`ht = ht.filter(...)` lines between where the offending expression was
bound and where it's used.

---

## `hl.literal(large_collection)` blows up JVM codegen and Spark RPC

**`[Hail]`**

Hail compiles literals into the query bytecode: a `hl.literal({...})`
with hundreds of thousands of elements produces a code object that
takes JVM codegen minutes of CPU, or OOMs before any Hail work starts.
The Dataproc analogue: a large literal used inside a *per-row*
expression is replicated into every task's serialized closure,
overflowing `spark.rpc.message.maxSize` (default ~128 MB, raised to
~835 MB in some gnomad_qc configs) on pipelines with tens of thousands
of partitions.

**Two failure modes, same root cause:**
- Local: `hl.literal(<~730k-element set>)` in a unit test hung the JVM
  for minutes before returning — no error, no progress bar, just CPU.
  Bigger sets OOM the driver.
- Dataproc: an `hl.literal(<set of sample indices>)` embedded in a
  per-row count expression got replicated into every task and killed the
  job at task-scheduling time with a serialized-task-too-big error.

**Fix pattern:** keep the large collection in the *data plane* (a
Table you join or semi-join against), not the *code plane*.

```python
# BAD: literal in a row-scoped expression → replicated to every task
mt = mt.annotate_rows(
    in_set=hl.literal(sample_index_set).contains(mt.sample_idx)
)

# GOOD: put it in a Table, join once, checkpoint if the set is big
set_ht = hl.Table.parallelize(
    [{"sample_idx": s} for s in sample_index_set],
    hl.tstruct(sample_idx=hl.tint32),
    key="sample_idx",
)
mt = mt.annotate_rows(in_set=hl.is_defined(set_ht[mt.sample_idx]))
```

If you *must* use `hl.literal` for a large set (small keep-set, one-off
transform), materialize it once inside a checkpointed transform and
hand the checkpointed HT to downstream steps — never inside a
`map`/`filter` lambda that runs per row.

**Symptom:** local Hail hangs before printing anything Hail-related
(JVM CPU spikes for minutes with no `Stage 0` line); Dataproc job dies
with `RpcSizeExceededException` or serialized-task-too-big errors,
often at task-scheduling time before any row work begins.

---

## `ClassTooLargeException` from inlining many sub-expressions into one row expression

**`[Hail]`**

Hail compiles a row expression into JVM bytecode. Building one expression
that inlines *K copies* of a heavy sub-expression — a Python loop that
stuffs a per-stratum breakdown (or a big `hl.case` / dict) into a single
`.annotate` / `.select` — compiles to one method/class that overflows the
JVM **64 KB class-size limit**. It's a scaling cliff, not a gradual
slowdown: the same code works at K≤2 and dies at K≥3 with
`is.hail.relocated.org.objectweb.asm.ClassTooLargeException: Class too
large: __C…native_writer`.

Concrete case: a per-group count built
```python
hl.dict([
    (g, hl.struct(raw=_count_from_sets(...), adj=_count_from_sets(...)))
    for g in grps
])
```
inside one `vp.select(gt_counts_by_grp=…)`. `_count_from_sets` is ~dozens
of set intersections; 2 groups compiled fine, 3 groups (`afr,nfe,eas`) blew
the class limit before any row work ran.

**Fix:** don't inline all K into one expression. Compute each in a
**separate pass / table**, checkpoint, then join-assemble the pieces:
```python
# BAD: K heavy sub-exprs inlined into one row expression → ClassTooLarge at K≥3
ht = ht.annotate(by_stratum=hl.dict([(s, heavy_expr(s)) for s in strata]))

# GOOD: one (checkpointed) table per stratum, then join into the dict
per = {s: compute(restrict_to(ht, s)).checkpoint(tmp(s)) for s in strata}
ht = ht.annotate(
    by_stratum=hl.dict([(s, per[s][ht.key].counts) for s in strata])
)
```
The K separate joins are separate IRs (fine); the assembled dict of K
struct-lookups is a small expression. See the per-stratum pattern under
[Codebase / analysis patterns](analysis-patterns.md#codebase-and-analysis-patterns).

**Symptom:** `ClassTooLargeException` naming a generated `__C…` class
(often `…native_writer`), raised at compile/execute time before rows are
processed; the same code path works with fewer strata / `hl.case`
branches and fails as you add more. Not memory, not data — pure codegen
size.

---

## `hl.experimental` numeric kernels don't validate inputs — check dtype and value invariants yourself

**`[Hail]` · `int32` requirement verified on Hail 0.2.134; the hang is reported from a production run**

Experimental / numeric Hail functions (EM solvers, statistical kernels)
are thin wrappers over a numeric routine with far fewer guardrails than
core Hail: they assume a specific input dtype, assume the input already
satisfies the routine's invariants, and don't cap iterations. Three
failure *classes*, none documented, all silent or misleading:

- **Strict, undocumented dtype requirements.** The kernel wants one
  specific numeric type and errors — or worse, misbehaves — on the
  wrong one. gnomAD aggregations default to `int64` / `float64`
  (`hl.agg.sum` / `count` promote), so you frequently need an explicit
  cast the signature doesn't advertise.
- **No input validation + uncapped iteration.** Iterative solvers
  assume well-formed input and have no max-iteration bound. Feed one
  input that violates an invariant it relies on (negative counts,
  totals exceeding the cohort) and it *oscillates instead of
  converging* — spinning forever on a single partition with no error
  and no progress, indistinguishable from a hung executor.
- **Degenerate input → NaN / Inf, not an error.** Division-based
  statistics silently return NaN on a zero denominator and propagate it
  into whatever consumes the result.

`hl.experimental.haplotype_freq_em` exhibits all three at once: it
requires `array<int32>` (cast `gt_counts.map(hl.int32)`); its iterations
are uncapped, so a set-encoding bug that yields negative cells or
`sum(gt_counts) > n_samples` makes it loop forever; and its
`p_chet = (h1·h2) / ((h0·h3) + (h1·h2))` is NaN when there are zero
double-het pairs, so guard with `hl.if_else(hl.is_nan(p),
hl.missing(...), p)`.

**Symptom:** a step using an experimental kernel either errors on a type
mismatch, or stalls for hours at `(N-1)/N` on one partition with no
traceback (looks exactly like the data-skew case under
[Dataproc](dataproc-operational.md#dataproc--operational)), or emits NaNs that quietly poison a
downstream aggregation.

**Rule of thumb:** treat any `hl.experimental` / numeric kernel as
un-validated. Before calling: cast inputs to the exact dtype it wants;
assert the value invariants it depends on (non-negativity,
`sum <= n_samples` — canary #1); and bound partition count so one
malformed row caps worst-case wall time instead of hanging the job.
Guard NaN on the output of any ratio.

---

## Aggregating across samples? Localize entries, then hoist each field into its own array

**`[Hail]` + `[gnomAD data]` · Verified in current `gnomad_methods` code**

To aggregate across samples per variant without paying the MatrixTable
machinery, the frequency code turns the MT into a Table whose rows carry
an array of every sample's entry:

```python
ht = mt.localize_entries("entries", "cols")
```

That much is a well-known trick. The part that isn't: **do not aggregate
by reaching into that array of structs.** Project each field you need
into its **own flat array** first, once, and aggregate over those:

```python
# gnomad_methods gnomad/utils/annotations.py :: agg_by_strata
ht = ht.select(
    *select_fields,
    **{ann: ht.entries.map(lambda e: e[ann]) for ann in select_expr.keys()},
)
```

The comment on that line in the source is explicit that this exists for
memory:

> Pull out each annotation that will be used in the array aggregation
> below as its own ArrayExpression. This is important to prevent memory
> issues when performing the below array aggregations.

**Why it matters:** an array-of-structs keeps the whole entry struct
live for every sample while you aggregate, and each field access inside
the aggregation walks that struct again. Hoisting gives you arrays of
exactly the values being aggregated — one narrow array per field instead
of one wide array of everything. On a release-scale callset that is the
difference between a job that fits and one that blows up. Gating the
`select_entries` down to only the fields the aggregation needs, *before*
localizing, is the same idea applied a step earlier.

**Rule of thumb:** the shape of the expression, not just the volume of
data, decides whether an aggregation fits in memory. If you're iterating
a per-row array of structs, project first.

---

## `ht[key_expr]` is a join, not an indexed point lookup

**`[Hail]`**

`other[ht.key]` / `ht.annotate(x=other[ht.key].x)` compiles to a join against
`other`, not an indexed random-access lookup, so "just annotating one field"
from a multi-GB table costs a join. Whether the left side shuffles depends
only on key structure: a lookup on a prefix of the left's key is an ordered
join, and anything else re-keys and shuffles the whole left side (next entry).

When it does shuffle, on Dataproc this is where boot-disk sizing bites:
secondary / preemptible workers with small (~40 GB) boot disks and several
executors per node sharing that disk overflow when a shuffle spills > ~2 GB,
killing the stage.

**Plan around it:**
- Look up on a prefix of the left's key where you can.
- **Disk, not just memory, is a scaling limit on Dataproc** — size worker
  boot disks (or use SSDs) for the shuffle, or shrink / shard the join.
- Co-partition tables you join repeatedly: read them all with the intervals
  from `gnomad.utils.file_utils.repartition_for_join` (see
  [Partitioning](partitioning-shuffles-oom.md#partitioning-shuffles-and-oom)).
- Chaining many non-prefix `[key]` lookups chains many shuffles; assemble via
  one keyed join where you can.

**Symptom:** an `annotate` that "should be cheap" spends its wall time in
a shuffle stage; the stage dies disk-full / `FetchFailed` on secondary
workers.

---

## A `ht[key_expr]` lookup on a **prefix** of the left's key doesn't shuffle

**`[Hail]` · Verified on Hail 0.2.134**

Hail decides shuffle-versus-ordered-join from **key structure alone**. If the
index expressions form a prefix of the left table's key, the join is ordered
and the left is never re-keyed. If they don't — *even when every field is a key
field* — Hail re-keys and shuffles the whole left side.

```python
# left key: (locus, alleles, transcript_id, uniprot_id, gene_id)
plp_by_locus_alleles[ht.locus, ht.alleles]           # prefix     -> ordered, no shuffle
an_by_locus_transcript[ht.locus, ht.transcript_id]   # NOT prefix -> shuffles the left
```

The second looks harmless — both fields are in the key — but `alleles` sits
between them, so it is not a prefix.

**Telling them apart before you run anything:** count `TableKeyBy` nodes in the
plan.

```python
str(ht.annotate(x=right[...])._tir).count("(TableKeyBy")
```

A prefix join adds none for the left; a non-prefix join adds three.

Measured on a 292,944,386-row table: two prefix lookups carried the partition
count through unchanged (4109 in, 4109 out) with no ordering shuffle logged and
the whole pass — read, both joins, full rewrite — taking ~7 minutes. The same
shape of lookup on a non-prefix key re-keyed all 292,944,386 rows.

**`hl.or_missing` does not avoid it.** Masking the result leaves the join in the
plan; the shuffle is paid and then discarded:

```python
hl.or_missing(cond, right[ht.a, ht.b].x)   # TableFilter: 0 — still shuffles
```

To keep rows that cannot match out of the shuffle, **filter before the join and
join the result back on the left's own key** (which is a prefix of itself, so
that step is ordered). This matters most when the unmatchable rows share one
key value: a sentinel like `"None"` in a key field passes `hl.is_defined`, so
every null guard misses it, and all of those rows hash to a single partition.

**Symptom:** one task in a stage runs while every other finishes in seconds, and
Spark launches a speculative duplicate that cannot help — both copies read the
same oversized partition.

---

## `other[ht.key]` keeps one row per key when the right table has duplicates

**`[Hail]` · Verified on Hail 0.2.134**

A non-interval `other[key]` lookup is a left join against a right table with
*distinct* keys — Hail lowers it to `TableLeftJoinRightDistinct`. If `other`
has several rows with the same key, the lookup returns one of them and
silently ignores the rest. `other.index(key, all_matches=True)` (marked
experimental) returns an array of every match.

```python
# right has two rows keyed 1 ("a", "b") and one keyed 2 ("c")
left.annotate(x=right[left.k].x)                           # k=1 -> "a" only
left.annotate(xs=right.index(left.k, all_matches=True).x)  # k=1 -> ["a", "b"]
```

**Symptom:** counts a little low with no error. A lookup table with 39
duplicated keys silently dropped 45 rows from a downstream count.

**Rule of thumb:** before relying on a lookup, check that the right side's key
is unique (`other.distinct().count() == other.count()`), or use
`all_matches=True` and handle the array.

---

## `.order_by()` destroys the key — use `add_index` to rekey cheaply

**`[Hail]`**

`ht.order_by(expr)` returns an unkeyed table. Any subsequent
`ht[key_expr]` join will fail (`Table is not keyed`), and even
`.key_by(...)` after the sort requires a *full pass* to rekey. This is
easy to trip over when computing ranks: sort by score, annotate
`hl.scan.count()` to get 0-based ranks, then try to join the rank back
onto the original table.

**Cheap-rekey pattern**: assign a stable integer index *before* the
sort, then key by that integer after the sort. Integer keys are
zero-overhead to rejoin because you already have the key in-hand from
either side:

```python
# BAD -- re-keying after order_by triggers a full-table shuffle
ranked = ht.order_by(ht.score).annotate(rank=hl.scan.count())
ranked = ranked.key_by(*ht.key)   # slow: has to re-shard by the original key
ht = ht.annotate(rank=ranked[ht.key].rank)

# GOOD -- add an integer index before sorting, use it to rejoin
ht = ht.add_index("_row_idx")
ranked = (
    ht.order_by(ht.score)
      .annotate(rank=hl.scan.count())
      .key_by("_row_idx")
)
ht = ht.annotate(rank=ranked[ht._row_idx].rank)
ht = ht.drop("_row_idx")
```

**Missing values sort last.** `order_by` puts every missing value after the
non-missing ones (verified on 0.2.134), so in the pattern above rows with a
missing `score` take the top ranks. Rank only the rows with a defined value —
filter first, and use the filtered row count as the bin denominator. Ranking
with the missing rows left in is what gave the v4.1 and v4.1.1 LOEUF releases
wrong rank and bin fields (fixed in v4.1.2).

**When to reach for this**: computing per-row ranks, quantile bins, or
any other order-dependent annotation you need to fold back onto the
original table. Also useful when you have many independent "sort by
column X, compute stat, join back" passes: run each on a `select`ed
narrow table (`ht.select("_row_idx", _val=X)`) plus the shared index,
checkpoint the tiny result, then join each back in one wide pass.
Avoids checkpointing the full wide table once per rank iteration.

**Symptom:** unexpected `Table is not keyed` after an `order_by`; or a
per-column rank loop that scales super-linearly with the number of
columns because each iteration is checkpointing the full wide table.

---

## `hl.scan.*` is the aggregator over everything *before* the current row

**`[Hail]` · Verified on Hail 0.2.134**

`hl.scan.<agg>` mirrors `hl.agg.<agg>` but produces a per-row running
value over all **preceding** rows in the table's current order,
**excluding** the current row:

```python
ht = hl.utils.range_table(5)
ht = ht.annotate(cum=hl.scan.sum(ht.idx), rank=hl.scan.count())
# cum   -> [0, 0, 1, 3, 6]     (row i sees rows 0..i-1)
# rank  -> [0, 1, 2, 3, 4]     (0-based row index)
# hl.agg.sum(ht.idx) for contrast -> 10 (one value for the whole table)
```

The exclusivity is the part people get wrong: `hl.scan.sum` on row 0 is
the aggregator's *empty* value, not that row's own value. For an
inclusive running total, add the current row back: `hl.scan.sum(x) + x`.

Scans follow the table's current row order, so an `order_by` before the
scan is what makes a rank meaningful — and `order_by` drops the key (see
the `.order_by()` entry above). `hl.scan.count()` is the cheap way to
get a 0-based row index for ranking.

`hl.scan` is thinly documented; the aggregators it accepts are the same
set as `hl.agg`.

---

## A Python function in `.map()` runs once, at plan-build time — not per row

**`[Hail]` · Verified on Hail 0.2.134**

Hail compiles expressions to JVM bytecode; it cannot call back into your
Python per row. A Python callable passed to `.map()` / `.filter()` is
invoked **once**, while the query plan is built, with an *expression* as
its argument — its job is to return an expression:

```python
calls = []
def py(x):
    calls.append(1)     # runs once for the whole array, not per element
    return x + 1
hl.eval(hl.array([1, 2, 3]).map(py))   # [2, 3, 4]; len(calls) == 1
```

That works only because `x + 1` is expressible in Hail. Anything with
real Python semantics inside — a regex library call, a lookup against a
Python object, an `if` on the *value* — either fails at plan-build time
or silently captures a constant.

**What to reach for instead:**
- Express it in Hail (`hl.if_else`, `hl.case`, `hl.rbind`, the string /
  math function library).
- `hl.experimental.define_function(f, *param_types)` — packages an
  expression-returning function as a reusable Hail-native function.
- Genuinely non-expressible per-row Python → run it on the driver over a
  *small* aggregate (size it first — that's the "don't pull a large table
  to the driver" problem), or, for a large table, drop to Spark and run it
  on the executors (next entry).

---

## Arbitrary Python per row at scale: drop to Spark and come back

**`[Hail]` · Round trip verified locally on Hail 0.2.134; gotchas from production runs (cited)**

When per-row work needs real Python — scoring with a pickled XGBoost model,
parsing thousands of structure files with Biopython — and the table is too big
for the driver, hand the rows to Spark, run the Python on the executors, and
convert back. This needs Hail's Spark backend (Dataproc or local):
`to_spark`, `from_spark`, and `hl.spark_context()` don't exist on
Query-on-Batch.

```python
import pickle

from pyspark.sql import Row
from pyspark.sql.types import DoubleType, StringType, StructField, StructType

spark = hl.utils.java.Env.spark_session()
sc = hl.spark_context()

# 1. Out of Hail: flat columns, and the locus as a string (to_spark() would
#    split it into locus.contig / locus.position columns).
ft = ht.key_by()
ft = ft.select(locus_str=hl.str(ft.locus), **{c: ft[c] for c in features})
ft.to_spark().write.mode("overwrite").parquet(in_path)
df = spark.read.parquet(in_path)  # plain Spark DataFrame, no Hail lineage

# 2. Python on the executors: broadcast the model once, then score each
#    partition (batch its rows into a pandas frame for a vectorized predict).
bc_model = sc.broadcast(pickle.dumps(model))


def score_partition(rows):
    m = pickle.loads(bc_model.value)
    for r in rows:
        x = [[r[c] for c in features]]
        yield Row(locus_str=r.locus_str, score=float(m.predict_proba(x)[0, 1]))


schema = StructType(
    [StructField("locus_str", StringType()), StructField("score", DoubleType())]
)
scored_df = spark.createDataFrame(df.rdd.mapPartitions(score_partition), schema)
scored_df.write.mode("overwrite").parquet(out_path)

# 3. Back into Hail.
scored = hl.Table.from_spark(spark.read.parquet(out_path))
scored = scored.annotate(
    locus=hl.parse_locus(scored.locus_str, reference_genome="GRCh38")
)
scored = scored.key_by("locus").drop("locus_str")
```

**Starting from files instead of rows:** `spark.read.format("text").load(glob,
wholetext=True)` gives one row per file, with the path in
`_metadata.file_path`. Parse each file in a pandas UDF or `mapPartitions`, then
take the same parquet + `from_spark` step back. This is how per-protein
sequence, distance-matrix, pLDDT, and pAE tables were built from thousands of
AlphaFold2 structure files.

**Gotchas:**
- **Write parquet on both sides.** Working directly on the Hail-backed
  DataFrame that `to_spark()` returns can hang (seen in production scoring
  runs). On the way back, `from_spark` without an
  intermediate parquet write "sometimes has missing and/or duplicate rows"
  (`gnomad/variant_qc/random_forest.py`).
- **Packages must be on the workers, and `--packages` doesn't put them
  there.** hailctl's init script installs `--packages` (and Hail itself) on
  the master only. Install executor-side packages with a Dataproc init action
  instead (`hailctl dataproc start --init=...`); see the packages entry in
  [Dataproc / operational](dataproc-operational.md#dataproc--operational).
- **The executors don't have Hail.** Any module the executor has to import to
  run your function — the module that defines it, and the modules its globals
  come from — must import cleanly without `hail`. Otherwise every task fails
  with `ModuleNotFoundError: No module named 'hail'` and the stage aborts after
  Spark's retries, even though the same code ran fine locally. Keep
  executor-side helpers in a module that doesn't import Hail.
- **`from_spark` can't take a `BinaryType` column.** Its docstring lists
  `BinaryType`, but Hail's type parser has no binary type, so the conversion
  fails with a misleading `IncompleteParseError ... begins with 'uct{...'`.
  Convert binary columns before `from_spark` — `pyspark.sql.functions.base64`
  turns one into a string, or decode it to an array. The conversion runs at the
  end of the job, so this loses the whole run (one SHAP run lost 6 hours).
- **Local runs: pin the workers' Python.** Spark starts executors with the
  `python3` on `PATH`; if that isn't the driver's interpreter the job fails with
  `PYTHON_VERSION_MISMATCH`. Set `PYSPARK_PYTHON` to the driver's Python before
  starting. (On Dataproc every node shares one Python.)
- **Partition count sets the per-task overhead.** Each task unpickles the model
  once, so many tiny partitions spend their time loading it and a few huge ones
  run out of memory. Batch rows within a partition rather than converting a
  whole partition to pandas at once.

---

## `hl.agg.count_distinct` does not exist

**`[Hail]` · Verified on Hail 0.2.134**

`hasattr(hl.agg, "count_distinct")` is `False`. The name is intuitive enough
that it gets written, and gets caught only at runtime — pylint won't flag it,
because `hl.agg` resolves dynamically.

```python
hl.len(hl.agg.collect_as_set(ht.gene_id))   # what you want instead
```

That builds the whole set in aggregator state, so it is fine for tens of
thousands of distinct values (gene or transcript IDs) and not for millions.

---

## `hl.literal` on an empty collection can't infer its type

**`[Hail]` · Verified on Hail 0.2.134**

`hl.literal(set())` raises `ExpressionException: Hail cannot impute type`. It
bites when the collection is *usually* non-empty and the empty case only shows
up on a small input or an edge-case test — a filter that matched nothing, an
optional gene panel.

```python
hl.literal(genes, hl.tset(hl.tstr))   # always pass the type
```

---

[← back to the index](README.md)
