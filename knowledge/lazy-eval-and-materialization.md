# Lazy evaluation and materialization

## Hail is lazy — an action re-runs everything upstream of it

**`[Hail]`**

Hail builds an unevaluated query plan; nothing executes until an
*action* forces it. The common actions: `count()`, `collect()`,
`take()`, `show()`, `aggregate()`, `to_pandas()`, `export()`,
`write()`, `checkpoint()`. (`describe()` is not one — it prints the
schema off the plan.)

The consequence that costs real money: an action on an
**un-checkpointed** table re-executes its entire lineage, and a second
action re-executes it *again*. A "harmless" logging line in the middle
of a pipeline —

```python
ht = expensive_join(...)
logger.info("%d rows", ht.count())   # <-- runs the whole join
ht = ht.filter(...)
ht.write(out)                        # <-- runs the whole join a SECOND time
```

— doubles the stage. Nothing about it looks expensive at the call site.

**AI assistants are especially prone to this**: they insert `count()` /
`show()` calls for progress logging and validity-checking, which is
harmless on a toy table and doubles wall-clock on a real one. When
reviewing generated code, grep for `.count()` and `.show()` and check
that each one sits *after* a checkpoint.

**Rule of thumb:** checkpoint first, then count. If you want a count
without keeping the table, that's still a full pass — decide it's worth
it deliberately.

---

## `ht.count()` after `checkpoint` is free — use it, don't avoid it

**`[Hail]` · Verified on Hail 0.2.134**

`ht.count()` on an unmaterialized table forces a full pass through the
query plan, so it's expensive on large HTs and the usual guidance is
"never count for logging on a big table." That guidance flips
completely once you've checkpointed:

```python
ht = build_expensive_ht(...).checkpoint("gs://scratch/x.ht")
n = ht.count()                   # <-- reads metadata from the checkpoint,
                                 #     doesn't re-run the upstream query
logger.info("checkpoint has %d rows", n)
```

A checkpointed HT is a materialized on-disk table in Hail's native
format (not Parquet) with a per-partition index; `count()` answers from
that metadata rather than scanning, and its wall time stays flat as the
table grows. You get a free logging / validity-check point. (`_force_count()`
is the private counterpart that *does* force the full pass — useful when
you actually want to materialize or time one.)

**Rule of thumb:** place `logger.info("%d rows", ht.count())` calls
*after* checkpoints, not before them. Also useful for the "did my
filter drop what I expected" check — `ht_after.count()` after a
checkpoint is essentially free; `ht_before.count()` before the
checkpoint costs a full extra pass.

**It's also the place to add a `naive_coalesce`.** `ht.n_partitions()` is free
here too, so you can see rows per partition. If a heavy filter left the
partitions sparse, coalesce *after* the checkpoint: the expensive work upstream
of it keeps its full parallelism (see
[Partitioning](partitioning-shuffles-oom.md#partitioning-shuffles-and-oom)).

**Symptom of the anti-pattern:** a pipeline where a "quick" logging
`.count()` doubles the wall-clock of a stage, because it's on an
un-checkpointed intermediate and forces a re-execution of everything
upstream.

---

## cache vs persist vs checkpoint

**`[Hail]` · Verified on Hail 0.2.128, 0.2.134, and 0.2.139 (source)**

`cache()` and `persist()` are **not** in-memory alternatives to a
checkpoint. `Table.persist` never passes its storage level on, and the
backend's `persist` (the Spark backend doesn't override it) calls
`checkpoint()` on a `persist_Table*` path under `hl.tmp_dir()`.

| | what it actually does | survives executor loss | deleted by |
|---|---|---|---|
| `ht.cache()` | `persist("MEMORY_ONLY")`, so the row below | if `tmp_dir` is durable (e.g. GCS) | `unpersist()` only |
| `ht.persist(level)` | ignores `level`; checkpoints under `hl.tmp_dir()` | if `tmp_dir` is durable | `unpersist()` only |
| `ht.checkpoint(path)` | `write()` + read back at `path` | yes | you (it's your path) |

`tmp_dir` defaults to `/tmp` on the default filesystem; with
`hl.init(tmp_dir="gs://…")` every `cache()` is a GCS write.

- **All three cost a full write.** Pick `checkpoint(path)` when you want
  to choose the location or reuse it across runs; `cache()` only saves
  you naming a path.
- **A library function that calls `.cache()` is checkpointing on the
  caller's behalf**, into their `tmp_dir`, and nothing removes it when
  the function returns unless it calls `unpersist()`.

`checkpoint()` takes `_read_if_exists=True` — private, but an
established gnomad_qc idiom (`_read_if_exists=not overwrite`) that lets
a re-run of a partially-completed pipeline skip already-written steps.

---

## `.cache()` before a shuffle is load-bearing, not just a perf optimization

**`[Hail]` · anecdotal · production runs**

Materializing the input to a `group_by` / big join with `.cache()` (or
`.checkpoint()`) isn't only about speed — it empirically avoids
intermittent shuffle failures (fetch-failures, timeouts) on large
datasets, because the shuffle reads a stable materialized source instead
of recomputing the upstream lineage under retry. Deleting a `.cache()`
that looks redundant (the intermediate is used once) reintroduces flaky
shuffle errors that only surface at scale.

**Symptom:** a `group_by` that intermittently dies with `FetchFailed` /
shuffle-fetch timeouts and sometimes succeeds with no code change — often
worse under executor churn.

**Rule of thumb:** keep `.cache()` before shuffle-heavy ops even when the
intermediate is used once; treat it as load-bearing, not decorative. For
very large intermediates prefer `.checkpoint("gs://…")` over `.cache()`:
both write, but `checkpoint` lets you pick a durable location, while
`.cache()` lands in `hl.tmp_dir()` (see the previous entry).

---

## A skip-if-exists wrapper has to defer *building* the expression, not just the write

**`[Hail]` · Verified on Hail 0.2.134**

`checkpoint(_read_if_exists=not overwrite)` skips the write. It does not skip
whatever ran while you were assembling the table handed to it — and in Python
that argument is evaluated first, before the wrapper can decide anything.

Most Hail table-building is lazy, so this is usually harmless. It stops being
harmless when a builder does real work while constructing the expression:

- an internal `.checkpoint()` — a helper that materialises an intermediate
- `hl.import_table(..., impute=True)` — scans the file to infer column types
- `.collect()` / `.aggregate()` — e.g. pulling a small table to the driver to
  rebroadcast it as a literal

Any of these run *before* the skip is considered, and their results are then
thrown away.

```python
_checkpoint(build_ht(src), path, resume=True)         # build_ht() already ran
_checkpoint(lambda: build_ht(src), path, resume=True) # only runs if needed
```

**Symptom:** a resumed run logs that it reused a checkpoint, but the log above
that line shows it wrote a large intermediate table first. Seen as a 402M-row
write on every resume of a pipeline whose first step checkpointed internally.

---

[← back to the index](README.md)
