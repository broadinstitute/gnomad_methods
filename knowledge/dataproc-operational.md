# Dataproc / operational

## Always `hl.copy_log()` to GCS in a `finally` — the local log dies with the cluster

**`[Hail]`/ops**

Hail writes its driver log to a local path on the Dataproc master
(`hl.init(log=…)`, default under `/tmp`). When the cluster is deleted —
or the job crashes, or your local `hailctl dataproc submit` loses its
network connection — that log is gone, which is exactly when you want
it. gnomad_qc's convention (`get_logging_path` + `hl.copy_log`) copies it
to GCS; put the copy in a `finally` so a *crashed* run's log is
preserved too:

```python
try:
    main(args)
finally:
    hl.copy_log(f"{tmp_dir}/logs/{script}.{label}.log")
```

`hailctl dataproc submit` streaming is **not** a durable record — a
dropped local connection stops the stream but the job keeps running on
the cluster, so the streamed stdout is not what you post-mortem from.
The `finally`-copied log is. This caught a run that hung for ~3h45m and
then died on an executor loss: the streamed client had long since
disconnected, and the GCS log was the only trace.

**Where to find more:** `get_logging_path` / `qc_temp_prefix` in
`gnomad_qc/v4/resources/basics.py`.

---

## Diagnosing a Dataproc job that "just died" (no Python traceback)

**`[Hail]`/ops**

A job that ERRORs with no Python traceback in the streamed output
usually failed JVM-side. Pull the driver log from GCS:

```bash
gcloud dataproc jobs describe <id> --format='value(driverOutputResourceUri)'
gsutil cat '<uri>*'   # quote the * so the local shell doesn't expand it
```

Find the *first* `Caused by`. Two common ones:

- **`ChunkFetchFailureException: … Executor is not registered (execId=N)`**
  — an executor died (OOM, or spot / preemptible-VM preemption) and a
  later task couldn't fetch its shuffle blocks. Not a code bug. Mitigate
  with non-preemptible workers, more executor memory, or by shrinking the
  shuffle.
- **A single task running for hours while the stage sits at `(N-1)/N`
  complete** — data skew: one partition inherited a pathological share of
  work (e.g. one key with a set ≈ n_samples in size). It often ends in
  the executor-loss above. Fix by sharding / repartitioning the skewed
  key, not by adding memory.

**Symptom:** `gcloud … describe` shows `ERROR`, but the streamed stdout
ended cleanly hours earlier and the progress bar was frozen at
`(N + 1) / N`.

**Gotcha while polling job state:** a `describe` loop that does
`2>/dev/null` turns a persistent auth / network failure into an
indistinguishable "empty state," so it can spin for hours (or days)
reporting nothing while the job already finished. Don't swallow stderr;
treat an empty state as *unknown*, not "still running."

---

## Standalone Dataproc scripts need `hl.init(tmp_dir="gs://...")`

**`[Hail]`/ops**

`hailctl dataproc submit` does **not** wire a GCS scratch path into
`hl.init` for you the way gnomad_qc's shared entrypoints do (via
`qc_temp_prefix`). If your script calls `hl.init()` bare, Hail defaults
its scratch directory to local HDFS `/tmp` on the *driver* node — tens
of GB of local disk, not shared with executors. Any checkpoint,
`hl.utils.new_temp_file()`, or large shuffle spill that outgrows that
fails with:

- `org.apache.hadoop.ipc.RemoteException: could only be replicated to 0
  nodes instead of minReplication`
- `java.io.EOFException: Premature end of file`

Both look like spot-VM preemption or an executor loss. They aren't —
they're the driver running out of local scratch. On a small cluster
this can fire at surprisingly modest sizes (a few-GB checkpoint on a
`n1-standard-4` driver).

**Fix:**
```python
hl.init(
    tmp_dir="gs://<your-scratch-bucket>",
    default_reference="GRCh38",
)
```

**Rule of thumb:** every standalone script that runs on Dataproc should
pass `tmp_dir=gs://...` to `hl.init`. gnomad_qc's shared entrypoints do
this via `qc_temp_prefix`; scripts that don't inherit from that
scaffolding need it explicitly.

**Symptom:** a checkpoint or shuffle write on a moderately-sized HT
fails with an HDFS `minReplication` / `Premature end of file` error;
the failed stage in the DAG is the checkpoint or shuffle-write itself
(not the compute upstream); job history shows no executor loss.

---

## Installing packages on a Dataproc cluster: driver-time pip vs `--packages` at start

**`[Hail]`/ops**

`hailctl dataproc submit`'s image ships Hail, numpy, pandas, scipy,
sklearn — but no ML / DL / bio libraries beyond that. A script that
imports `xgboost` / `lightgbm` / `shap` / `torch` / a specific
`pyarrow` version at module top-level dies with `ModuleNotFoundError`
on load, even if it worked locally in your conda env.

**Fix — pip-install at driver import time**, before any transitive
import that needs the package:

```python
import subprocess, sys
try:
    import xgboost  # noqa: F401
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "xgboost"])
```

Works for driver-only workloads (fit-a-model-on-collected-pandas,
driver-side scoring, etc.). For executor-side use — e.g. broadcasting
a model into a `mapPartitions` scorer — the package has to be on the
**workers**, and only cluster creation can put it there:

```bash
hailctl dataproc start CLUSTER \
  --packages xgboost,shap \        # --pkgs is the same flag; installed on all nodes
  --properties spark:spark.driver.maxResultSize=0 \   # arbitrary Spark config
  --metadata VEP_CONFIG_URI=…      # consumed by the init actions
```

`hailctl dataproc submit` has **no** `--packages` equivalent, so a
driver-time pip-install never reaches the workers. Decide which side
needs the package *before* you start the cluster; otherwise it's a
cluster rebuild.

**Symptom:** job dies with `ModuleNotFoundError` on `import <package>`
at load time, before any Hail work starts. `hailctl dataproc submit`'s
streamed output ends at the traceback with no Spark stage having run.
If the driver imports fine and executors fail later, the package went
on the wrong side.

---

## Anything that pulls a large HT / MT to the driver silently stalls — filter Hail-side first

**`[Hail]`**

`hl.Table.to_pandas()`, `.collect()` (row / expression / global),
`hl.agg.collect(...)` inside an `.aggregate(...)` — all of these
materialize their result *into driver memory*. On a large table
(tens of millions of rows, or a per-row array collected across many
rows), the Spark stage backing the collect sits at N-1/N for hours,
sometimes returns, sometimes cascades into executor-loss. There's
no OOM message; the driver just becomes unresponsive.

**Fix — narrow it Hail-side first.** The specific technique depends
on what you actually need on the driver:

```python
# 1. You need a small subset of rows → semi_join (NOT hl.set().contains,
#    which broadcasts the whole key set into every task; see the
#    hl.literal / large-collection entry above).
key_ht = hl.import_table(keys_tsv, ...).key_by("locus", "alleles").select()
df = ht.semi_join(key_ht).to_pandas()

# 2. You need per-partition summaries, not per-row values → aggregate to
#    a small struct before pulling.
stats = ht.aggregate(hl.struct(
    n=hl.agg.count(),
    mean_af=hl.agg.mean(ht.af),
))  # <- one struct, safe to collect

# 3. You need per-gene numbers → group_by + aggregate before to_pandas.
per_gene = ht.group_by(ht.gene_id).aggregate(n=hl.agg.count()).to_pandas()

# 4. You need a Python list of a small keyed field → export to GCS,
#    read back locally (uses Spark to materialize, doesn't blow driver RAM).
ht.select(ht.gene_id).export("gs://…/genes.tsv.bgz")
```

Rule of thumb: **anything you're about to `.to_pandas()` /
`.collect()` / `.aggregate(hl.agg.collect(...))` needs to be small
enough to fit in driver RAM** — a large HT / MT never is. If the
shape isn't obviously ≪ 10⁶ rows, use one of the Hail-side patterns
above.

**Symptom:** a Spark stage sitting at `(N-1)/N` for hours after a
`to_pandas()` / `collect()` line; no traceback, no OOM message; you
eventually kill the job. If the executor-loss cascade fires, the
driver *does* die — but hours later, with a
`ChunkFetchFailureException` that reads like a preemptible-VM issue
rather than the real "driver was collecting a huge HT" cause.

---

## Size the driver deliberately — and keep scoring on the workers

**`[gnomAD data]` · anecdotal; hailctl defaults verified**

`hailctl dataproc start` defaults: master `n1-highmem-8`, workers
`n1-standard-8` with **40 GB boot disks**, and
`--master-memory-fraction 0.9` (that fraction of master memory goes to
the JVM, the remainder to Python).

- **A bigger driver is usually cheap.** You pay for one node, not N.
  `-m n1-highmem-16` on a job that keeps stalling driver-side often
  costs less than the wall-clock lost to retries.
- **Driver-side Python needs that fraction lowered.** If the driver
  runs pandas / sklearn, the default 0.9 has already handed the memory
  to the JVM.
- **The 40 GB worker boot disk is the shuffle-spill ceiling**
  (`--worker-boot-disk-size`, `--num-worker-local-ssds`); see the
  `ht[key_expr]`-is-a-join entry.

**Don't collect in order to score.** The established gnomAD pattern for
model scoring hands the table to Spark ML and takes it back — it never
pulls the feature matrix to the driver:

```python
df = ht.key_by().select(*features).to_spark()   # gnomad_methods variant_qc/random_forest.py
# ... pyspark.ml pipeline: fit / transform, distributed across workers ...
rf_ht = hl.Table.from_spark(rf_df)
```

**Gotcha, documented in that same file:** write the Spark DataFrame to
parquet and re-read it before `hl.Table.from_spark` — without the
intermediate write the resulting HT "sometimes has missing and/or
duplicate rows."

---

## Use `gcloud storage`, not `gsutil`

**ops · verified against Google's docs, July 2026**

Google's current guidance: *"gsutil is not the recommended CLI for Cloud
Storage. Use `gcloud storage` commands in the Google Cloud CLI
instead."* gsutil is legacy and minimally maintained, and **after March
2027 it will no longer ship with the Google Cloud CLI**. `gcloud
storage` also "requires less manual optimization in order to achieve the
fastest upload and download rates" — which matters when you're moving
multi-TB HTs, where gsutil needed `-m` plus tuning to keep up.

Translation is mostly mechanical:

| gsutil | gcloud storage |
|---|---|
| `gsutil -m cp -r src dst` | `gcloud storage cp --recursive src dst` |
| `gsutil ls 'gs://…/*'` | `gcloud storage ls 'gs://…/*'` |
| `gsutil du -sh gs://…` | `gcloud storage du --summarize --readable-sizes gs://…` |
| `gsutil -m rsync -r a b` | `gcloud storage rsync --recursive a b` |

Two habits worth keeping either way: **quote globs** (`'gs://…/*'`) so
the local shell doesn't expand them, and prefer `rsync` over re-copying
when syncing a directory of checkpoints.

---

## Query on Batch (QoB) vs Dataproc

**`[Hail]`/ops**

Query on Batch is the same Hail you already write (`import hail as hl`)
with a different executor: instead of a Spark cluster you provision with
`hailctl dataproc start`, queries run on a shared, Hail-maintained Batch
service. Switching is one argument:

```python
hl.init(
    backend="batch",
    worker_cores=1, worker_memory="standard",
    driver_cores=1, driver_memory="standard",
)
```

| | Dataproc | QoB |
|---|---|---|
| user-provisioned cluster | yes | no |
| latency, short queries | low | low–medium (depends on cluster size) |
| all operations supported | yes | not all BlockMatrix ops |
| per-job logs visibility | no | yes |
| per-job cost visibility | no | yes |

The cost and log visibility is the practical difference: on Dataproc you
find out what a job cost by looking at the cluster bill afterwards; QoB
attributes it per job. Against that, a shared service means you don't
control the machine shape the way `--master-machine-type` lets you (see
the driver-sizing entry above).

**`[internal only]`** — the operating guidance for this group is to
prefer QoB where it works, and to raise cases where it doesn't with the
Hail team rather than quietly falling back to Dataproc.

**Not everything in this section transfers.** The Dataproc entries above
(boot disks, `--packages`, `hl.copy_log`, driver sizing) are
cluster-provisioning concerns that QoB takes over. The *Hail-level*
entries — shuffles, partitioning, laziness, missingness — apply
identically on both backends.

---

## Cloud storage settings quietly dominate the bill

**ops · current as of July 2025**

Four settings that cost money without appearing anywhere in your code:

- **Multi-regional buckets.** Data replicated across data centers costs
  substantially more to store and to move. Use regional buckets unless
  someone has deliberately decided otherwise — and check inherited or
  requester-pays buckets, which are where the surprises live.
- **Compute/data region mismatch.** A job running in one region reading
  a bucket in another pays inter-region transfer on every byte. Hail's
  Batch default region was effectively "any region" **before 0.2.135**,
  which made this easy to hit by accident; 0.2.135 sets a default.
- **Soft delete is ON by default on new buckets** — you keep paying for
  every object for 7 days *after* deleting it. Turn it off on scratch
  buckets, or the "cleanup" you just ran doesn't show up on the bill for
  a week.
- **Temp directories accumulate.** A remote temp dir wants a lifecycle
  policy of ~7 days or less, should be used *only* by Hail, and should
  be one-per-project so read/write costs bill to the right place.

`hailctl batch init` creates a temp bucket with these settings already
right (lifecycle on, soft-delete off, correctly labeled), and
`hailctl config profile` switches billing project, region, and temp dir
together — worth having one profile per project rather than editing
config by hand.

Deleting an old temp directory is `gcloud storage rm -R <path>`. Check
what you're pointing at first; it is not recoverable.

---

[← back to the index](README.md)
