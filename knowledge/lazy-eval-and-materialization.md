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
`show()` calls for progress logging and sanity-checking, which is
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
table grows. You get a free logging / sanity-check point. (`_force_count()`
is the private counterpart that *does* force the full pass — useful when
you actually want to materialize or time one.)

**Rule of thumb:** place `logger.info("%d rows", ht.count())` calls
*after* checkpoints, not before them. Also useful for the "did my
filter drop what I expected" check — `ht_after.count()` after a
checkpoint is essentially free; `ht_before.count()` before the
checkpoint costs a full extra pass.

**Symptom of the anti-pattern:** a pipeline where a "quick" logging
`.count()` doubles the wall-clock of a stage, because it's on an
un-checkpointed intermediate and forces a re-execution of everything
upstream.

---

## cache vs persist vs checkpoint

**`[Hail]` · Verified on Hail 0.2.134**

| | what it is | survives executor loss | survives the session |
|---|---|---|---|
| `ht.cache()` | alias for `persist("MEMORY_ONLY")` | no | no |
| `ht.persist()` | Spark block storage, default `MEMORY_AND_DISK` | no | no |
| `ht.checkpoint(path)` | `write()` + read back — a real on-disk table | yes | yes |

- **Small intermediate reused a few times in one session** → `cache()`.
  Also load-bearing before a shuffle (next entry).
- **Large intermediate, anything feeding a big shuffle, or anything
  you'd hate to recompute** → `checkpoint("gs://…")`. Durable across
  executor loss, and makes `count()` free.
- **`persist()`** only when you specifically want `MEMORY_AND_DISK`
  without paying a write.

`checkpoint()` takes `_read_if_exists=True` — private, but an
established gnomad_qc idiom (`_read_if_exists=not overwrite`) that lets
a re-run of a partially-completed pipeline skip already-written steps.

---

## `.cache()` before a shuffle is load-bearing, not just a perf optimization

**`[Hail]`**

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
very large intermediates prefer `.checkpoint("gs://…")` (spills to GCS,
survives executor loss) over in-memory `.cache()`.

---

[← back to the index](README.md)
