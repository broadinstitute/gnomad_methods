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
- Dataproc: an `hl.literal(<PBT-index-set>)` embedded in a per-pair
  count expression got replicated into every count task and killed the
  job at task-scheduling time with a serialized-task-too-big error.

**Fix pattern:** keep the large collection in the *data plane* (a
Table you join or semi-join against), not the *code plane*.

```python
# BAD: literal in a row-scoped expression → replicated to every task
mt = mt.annotate_rows(
    is_pbt=hl.literal(pbt_index_set).contains(mt.sample_idx)
)

# GOOD: put it in a Table, join once, checkpoint if the set is big
pbt_ht = hl.Table.parallelize(
    [{"sample_idx": s} for s in pbt_index_set],
    hl.tstruct(sample_idx=hl.tint32),
    key="sample_idx",
)
mt = mt.annotate_rows(is_pbt=hl.is_defined(pbt_ht[mt.sample_idx]))
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

Concrete case: a per-population count built
```python
hl.dict([
    (pop, hl.struct(raw=_count_from_sets(...), adj=_count_from_sets(...)))
    for pop in pops
])
```
inside one `vp.select(gt_counts_by_pop=…)`. `_count_from_sets` is ~dozens
of set intersections; 2 pops compiled fine, 3 pops (`afr,nfe,eas`) blew
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

**`[Hail]`**

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

## `ht[key_expr]` is a Spark join (shuffle), not an indexed point lookup

**`[Hail]`**

Hail's Spark backend compiles `other[ht.key]` /
`ht.annotate(x=other[ht.key].x)` into a **Spark join** — broadcast if the
right side is small, otherwise sort-merge with a shuffle to disk — not an
indexed random-access lookup. So "just annotating one field" from a
multi-GB table triggers a full shuffle. On Dataproc this is where
boot-disk sizing bites: secondary / preemptible workers with small
(~40 GB) boot disks and several executors per node sharing that disk
overflow when a shuffle spills > ~2 GB, killing the stage.

**Plan around it:**
- A right-side table > ~10 MB won't broadcast; it sort-merge joins
  (shuffle). Keep frequently-joined annotation tables small, or
  `checkpoint` + co-partition them.
- **Disk, not just memory, is a scaling limit on Dataproc** — size worker
  boot disks (or use SSDs) for the shuffle, or shrink / shard the join.
- Chaining many `[key]` lookups chains many joins/shuffles; assemble via
  one keyed join where you can.

**Symptom:** an `annotate` that "should be cheap" spends its wall time in
a shuffle stage; the stage dies disk-full / `FetchFailed` on secondary
workers.

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

## There are no Python UDFs — a Python function in `.map()` runs once, at plan-build time

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
- Genuinely non-expressible per-row Python → restructure so the Python
  runs on the driver over a *small* aggregate, or export and process
  outside Hail. That's the "don't pull a large table to the driver"
  problem; size it before you reach for it.

---

[← back to the index](README.md)
