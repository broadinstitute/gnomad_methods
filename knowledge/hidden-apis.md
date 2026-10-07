# Hidden and undocumented APIs we rely on

Private (underscore) parameters are unsupported and can change between
Hail versions — but several are load-bearing in gnomad_qc, so you'll
meet them in code review. Re-check these on a Hail upgrade.

| API | what it does |
|---|---|
| `read_table` / `read_matrix_table(_n_partitions=n)` | set partitioning at read, no shuffle |
| `read_*(_intervals=…, _filter_intervals=True)` | read only the given intervals |
| `checkpoint(_read_if_exists=not overwrite)` | skip already-written steps on a re-run |
| `expr.aggregate(..., _localize=False)` | keep the result as an expression instead of pulling it to the driver |
| `ht._force_count()` | force a real full pass (vs the metadata-only `count()`) |
| `ht._filter_partitions(range(n))` | read just the first `n` partitions of a real table — cheap iteration on real data |

Supported-but-easy-to-miss `hl.experimental` helpers the pipelines lean
on: `sparse_split_multi`, `densify`, `filtering_allele_frequency`,
`pc_project`, `import_gtf`, `get_gene_intervals`, `read_expression` /
`write_expression`, `define_function`. Note the input-validation caveat
in the experimental-kernel entry above before using the numeric ones.

---

## Backend feature flags: `hl._set_flags` / `hl._get_flags` / `hl._with_flags`

**`[Hail]` · Verified on Hail 0.2.139**

Hail carries backend feature flags that are not surfaced anywhere in the
public API or the docs. They are set as **strings** and cleared with `None`:

```python
hl._set_flags(use_new_shuffle="1")   # on
hl._set_flags(use_new_shuffle=None)  # off
hl._get_flags("use_new_shuffle")     # {'use_new_shuffle': '1'}
```

**They are undiscoverable by inspection.** `hl._get_flags()` with no
arguments returns `{}` — you only get back the flags you name. The way to
see the list is to set an invalid one on purpose: `_set_flags` raises a
`FatalError` whose message enumerates every valid flag (28 of them in
0.2.139). `_get_flags` on an unknown name is less helpful — it throws a raw
`Py4JJavaError: NoSuchElementException` instead.

**Use the context manager, not a pair of calls.** `hl._with_flags(**flags)`
reads the previous values and restores them in a `finally`:

```python
with hl._with_flags(use_new_shuffle="1"):
    ht.write(path)
```

Hand-rolled `_set_flags(x="1")` … `_set_flags(x=None)` pairs are common in
older code and have two problems: an exception between them leaks the flag
into the rest of the session, and they reset to `None` rather than to
whatever was set before. If you need to support a Hail older than the one
you verified `_with_flags` on, use `try` / `finally` rather than a bare pair.

**The shuffle-related names** are `use_new_shuffle`,
`shuffle_cutoff_to_local_sort` and `shuffle_max_branch_factor`.

**When `use_new_shuffle` earns its keep** — *`[gnomAD data]` · anecdotal,
from the v4 variant-pair-list work*: it selects Hail's own shuffler over
Spark's. On a mostly-preemptible fleet, Spark keeps shuffle blocks on
executor-local disk, so a reclaimed worker takes its blocks with it and the
stage dies with `FetchFailedException` — the failure that no partition
tuning fixes (see
[Partitioning, shuffles, and OOM](partitioning-shuffles-oom.md#partitioning-shuffles-and-oom)).
The recipe that got that step through genome-wide was all three together:
`use_new_shuffle` on the `group_by`, a materialisation before it, and fewer
preemptible workers. Setting the flag alone did not carry it; neither did
resizing alone.

---

[← back to the index](README.md)
