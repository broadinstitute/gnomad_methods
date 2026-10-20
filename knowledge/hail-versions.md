# Hail versions

## Newer Hail writes files that older Hail can't read

**`[Hail]` · Verified against the Hail change log and `FileFormat.version` in the source at the 0.2.134, 0.2.138, and 0.2.139 tags**

Hail reads every older native file format, but not newer ones: a Table or
MatrixTable (and so a VDS) written by a release that introduced a new format
version can't be opened by any earlier release (see the FAQ at the top of the
[change log](https://hail.is/docs/0.2/change_log.html)).

| format version | introduced in |
|---|---|
| 1.4.0 | 0.2.34 |
| 1.5.0 | 0.2.55 |
| 1.6.0 | 0.2.82 |
| 1.7.0 | 0.2.119 |
| 1.8.0 | 0.2.139 |

**1.8.0 is easy to miss.** The 0.2.139 change log has no "File Format" note
for it; the change appears only as a New Features bullet about index files
([hail#15445](https://github.com/hail-is/hail/pull/15445)).

**To see what wrote a table,** read `hail_version` and `file_version` from its
`metadata.json.gz`. `file_version` packs the format as
`major << 16 | minor << 8 | patch` (67328 is 1.7.0, 67584 is 1.8.0).

```python
import gzip
import json

import hailtop.fs as hfs

with hfs.open(f"{path}/metadata.json.gz", "rb") as f:
    meta = json.loads(gzip.decompress(f.read()))
meta["hail_version"], meta["file_version"]
```

**Symptom:** an `incompatible file format` error when an older Hail reads a
table. In gnomAD, a densify pipeline moved from 0.2.139 to 0.2.128 couldn't
read its own 0.2.139-written intermediates and had to regenerate all of them
(gnomad_production#2206).

**Rule of thumb:** pin one Hail version for a whole pipeline, and check
`hail_version` on its inputs before downgrading. One step run on 0.2.139
makes its outputs unreadable to 0.2.138 and earlier.

---

## Known problems by Hail version

**`[Hail]` · Fix versions checked against the Hail release tags**

Bugs and regressions that have cost gnomAD work. The change log lists fixes
but not always which versions were affected, and some of these aren't in it
at all. Add a row only with a link to its source.

| version | problem | fixed in | source |
|---|---|---|---|
| 0.2.139 | Scans with a large aggregator state, including `hl.experimental.densify`, crash with `NotImplementedError: SimpleSStream.settableTupleTypes` once a table has more partitions than `branching_factor` (default 50). Raising `branching_factor` avoids the crash but can OOM the driver on a large densify; gnomAD's run moved to 0.2.128 instead. | fix merged ([hail#15704](https://github.com/hail-is/hail/pull/15704)), not yet released | [hail#15700](https://github.com/hail-is/hail/issues/15700); gnomad_production#2206 |
| 0.2.139 | A Table filter using an array literal's `.contains()` leaks driver memory until the driver OOMs. Set-based membership isn't affected. gnomad_methods excludes 0.2.139 for this. No upstream Hail issue found. | — | gnomad_methods#845 |
| 0.2.138 | On the `batch` (Query-on-Batch) backend, requester-pays configuration is ignored, and `hl.vep` fails on requester-pays buckets. | 0.2.139 ([hail#15468](https://github.com/hail-is/hail/pull/15468), [hail#15470](https://github.com/hail-is/hail/pull/15470)) | [hail#15465](https://github.com/hail-is/hail/issues/15465) |
| 0.2.134–0.2.138 | Filtering the samples of a VDS (or columns of a MatrixTable) after a per-entry transformation materializes the full entries array for every column. A performance regression that came with the VDS `LEN` field ([hail#14675](https://github.com/hail-is/hail/pull/14675)). | 0.2.139 ([hail#15543](https://github.com/hail-is/hail/pull/15543)) | change log |
| ≤ 0.2.137 | `hl.vep` bound a parameter incorrectly, producing incorrect results. The change log doesn't say which versions were affected. | 0.2.138 ([hail#15337](https://github.com/hail-is/hail/pull/15337)) | change log |
| 0.2.128 | `import hail` needs `pkg_resources`, which `setuptools` 81 and later don't provide; add `setuptools<81` to the cluster's `--packages`. 0.2.134 and later don't import it. | — | gnomad_production#2206 |

**Rule of thumb:** check this table before choosing a Hail version for a
pipeline, and before upgrading one mid-project.

---

[← back to the index](README.md)
