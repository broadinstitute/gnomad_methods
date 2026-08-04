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

[← back to the index](README.md)
