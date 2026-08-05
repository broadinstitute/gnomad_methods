# Canary tests

Cheap sanity checks that catch whole classes of the bugs above. Run
after any refactor that touches the annotation join graph or the
aggregation semantics.

```python
# 1. Per-sample count arrays: non-negative, and the sum can't exceed
#    the cohort size (catches set-encoding / complement-storage bugs,
#    which otherwise surface only as a downstream EM that never
#    converges)
assert ht.aggregate(hl.agg.all(hl.sum(ht.gt_counts_adj) <= n_samples))

# 2. Aggregation totals match input (catches null-key group_by drops)
n_input = ht.count()
result = ht.aggregate(hl.agg.group_by(key, hl.agg.count()))
n_output = sum(result.values())  # includes None-key bucket
assert n_input == n_output

# 3. No unexplained nulls in critical annotations
n_null = ht.aggregate(hl.agg.count_where(hl.is_missing(ht.filters)))
# n_null should match your understanding of the coverage boundary

# 4. Coverage-boundary check for permissive vs release-scoped HTs
release_n = ht.aggregate(hl.agg.count_where(ht.filters.length() == 0))
permissive_n = ht.aggregate(hl.agg.count_where(
    ht.filters.difference(hl.set(["AC0"])).length() == 0
))
# permissive_n >= release_n; the gap is your AC0-in-release count
```

---

## Quick tests

Canary tests above check a *production* output. This is the cheaper
loop that comes first: a pure transform can be tested on a table you
type out by hand, locally, in seconds — no cluster, no GCS.

```python
ht = hl.Table.parallelize(
    [
        {"locus": ..., "af": 0.5,  "source": {"clinvar_plp"}},
        {"locus": ..., "af": None, "source": None},          # <- the null row
        {"locus": ..., "af": 0.0,  "source": set()},          # <- the empty-collection row
    ],
    hl.tstruct(locus=hl.tlocus("GRCh38"), af=hl.tfloat64, source=hl.tset(hl.tstr)),
)
assert my_transform(ht).aggregate(...) == expected
```

**Always include a null row, an empty-collection row, and a
boundary-value row.** Those are precisely the inputs the
[Missingness](missingness.md#missingness) section is about, and they're the ones a
hand-written happy-path fixture omits.

**This is a good task to hand to an AI assistant** — generating small
fixtures and their expected values is exactly the kind of mechanical
work it does well, and unlike a Dataproc run you can check the answer
immediately. Ask for the null / empty / boundary cases explicitly;
they're the ones that get skipped.

**When you do need real data, take a slice of it.**
`ht._filter_partitions(range(2))` reads only the first two partitions of
a real table — enough to exercise a transform against genuine schema and
genuine nulls without paying for the whole table. Private API; see
[hidden APIs](hidden-apis.md#hidden-and-undocumented-apis-we-rely-on).

**Keep pipeline code testable this way:** utils functions as pure
transforms (tables in → tables out) with file I/O confined to the CLI
entrypoints. A transform that reads its own inputs can only be tested
against real data on a cluster.

---

[← back to the index](README.md)
