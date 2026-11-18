# Sequencing Lane Demultiplexing API

FastAPI + SQLite service that splits mixed paired-end gzip FASTQ uploads into
per-sample files. Designed to run as **two container instances on one shared
volume** with a **128 MiB memory limit**.

## What it does

* Accepts a lane number, a dual-index barcode table (1–8 samples, two
  8-base barcodes each) and paired `R1`/`R2` gzip FASTQ files.
  Decompressed input across both mates is capped at **512 MiB**.
* Validates both mates: valid gzip, well-formed FASTQ, equal record counts,
  identical read names in the same order, matching sequence/quality lengths.
  Any corruption fails the **whole lane**.
* Assigns a pair to a sample when the combined R1+R2 barcode matches with
  **at most one substitution** and the candidate sample is **unique**.
  Everything else (no hit, >1 substitution, or a tie) goes to quarantine.
* Assigned reads have both 8-base barcodes and their corresponding quality
  characters trimmed; quarantined reads are kept intact.
* A lane becomes visible only when **all** per-sample paired files, the
  quarantine paired files and the count table exist (atomic publish:
  files staged, renamed, then a single DB transaction).
* Same lane + identical input (hash of barcode table and both files) reuses
  the existing task; same lane with different input is rejected with `409`.
* Failed tasks can be retried (`POST /jobs/{id}/retry`); uploads are kept
  until a lane is published.
* Snapshots freeze the published lane list at creation time. Sample/
  quarantine downloads through a snapshot only ever concatenate files from
  that fixed list — lanes published later never mix in.
* **On-demand read-name indexes** let the data-lake team recover the full
  paired-end record(s) behind a suspicious read name from a published lane
  without rescanning the FASTQ. See [Read-name indexes](#read-name-indexes).

## Running (Docker Compose, two instances on a shared volume)

```bash
docker compose up --build
# instance A: http://localhost:8001   instance B: http://localhost:8002
# OpenAPI demo UI:  http://localhost:8001/docs  (also on 8002)
```

Both containers mount the named volume `demux-shared` at `/data`, are capped
at 128 MiB (`mem_limit`), and automatically recover on restart.

### Quick demo

```bash
python scripts/make_demo_fastq.py /tmp/demo          # makes R1/R2.fastq.gz
# In http://localhost:8001/docs -> POST /jobs:
#   lane: demoLane
#   barcodes: [{"sample_id":"sampleA","barcode_r1":"AAAAAAAA","barcode_r2":"CCCCCCCC"},
#              {"sample_id":"sampleB","barcode_r1":"GGGGGGGG","barcode_r2":"TTTTTTTT"}]
#   r1/r2: the generated gzip files
curl http://localhost:8001/lanes/demoLane             # counts
curl -OJ http://localhost:8001/lanes/demoLane/samples/sampleA/r1
curl -OJ http://localhost:8001/lanes/demoLane/quarantine/r1
```

## API surface

| Method | Path | Purpose |
|---|---|---|
| POST | `/jobs` | Submit a lane (multipart: `lane`, `barcodes` JSON, `r1`, `r2`) |
| GET | `/jobs`, `/jobs/{id}` | Task status (`queued`/`processing`/`published`/`failed`) |
| POST | `/jobs/{id}/retry` | Requeue a failed task |
| GET | `/lanes`, `/lanes/{lane}` | Published lanes and per-bucket counts |
| GET | `/lanes/{lane}/samples/{sample_id}/{r1\|r2}` | Download a sample mate |
| GET | `/lanes/{lane}/quarantine/{r1\|r2}` | Download a quarantine mate |
| POST | `/snapshots` | Freeze the current published lane list |
| GET | `/snapshots`, `/snapshots/{id}` | List/inspect snapshots |
| GET | `/snapshots/{id}/samples/{sample_id}/{mate}` | Snapshot-scoped concatenated download |
| GET | `/snapshots/{id}/quarantine/{mate}` | Snapshot-scoped quarantine download |
| POST | `/lanes/{lane}/samples/{sample_id}/index` | Build (or rebuild) a sample's read-name index |
| POST | `/lanes/{lane}/quarantine/index` | Build the quarantine bucket's index |
| GET | `/lanes/{lane}/samples/{sample_id}/index` | Index build status (`building`/`ready`/`failed`) |
| GET | `/lanes/{lane}/quarantine/index` | Quarantine index build status |
| GET | `/lanes/{lane}/indexes` | List every index task for a lane |
| GET | `/lanes/{lane}/samples/{sample_id}/reads?read_name=...` | Retrieve all pairs with that name |
| GET | `/lanes/{lane}/quarantine/reads?read_name=...` | Quarantine-scoped name lookup |
| GET | `/health` | Liveness, instance id, current job, queue depth |

## Read-name indexes

When anomaly comparison flags a read name, the engineer first asks the
service to materialise an **on-disk index** for the relevant published
bucket (a per-sample pair of files, or the quarantine pair), then queries
by name.

### Lifecycle and guarantees

* **On-demand, background builds.** `POST .../index` enqueues an index task
  (`202`) and returns immediately; the background worker builds it without
  blocking uploads, queries or downloads for any other lane. A repeat POST
  while one exists returns the same task (`200`); POSTing a `failed` task
  requeues a rebuild.
* **Sequential verification, whole-bucket failure.** A build streams both
  published gzip mates record by record and checks, in order: (1) every
  file exists and is a complete, well-formed gzip FASTQ; (2) R1 and R2
  read-name tokens match in the same order and record counts are equal;
  (3) the indexed pair count equals the published count table. Any missing
  file or mismatch fails the entire index — the bucket never gets a partial
  index.
* **Record size limits.** Only records whose read-name token is at most
  **256 bytes** and whose sequence/quality lines (either mate) are at most
  **1000 bytes** are indexed. An oversized record aborts the build with an
  explicit error.
* **Queryable only once fully on disk.** Builds run in a hidden
  `indexes/.build.<id>/` staging directory and are atomically `rename`d to
  `indexes/<id>/` only after the data is fsynced and a `COMPLETE` marker is
  written. A lookup while absent/building returns `409` (or `404` if no index
  was ever requested); a failed build returns `409` with the reason.
* **No rescans on repeated queries.** Each index stores the records as plain
  four-line FASTQ in `r1.store`/`r2.store` plus an SQLite table mapping
  `name → ordinal → byte offsets`. A lookup seeks straight to the records;
  the source gzip is never opened again.
* **Duplicate names return every pair, in order.** A name shared by several
  read segments yields all of them, ordered by their original ordinal.
* **Crash and concurrency safety.** A build holds an exclusive `flock` on its
  staging directory; startup recovery fails interrupted tasks and removes
  their debris only when the lock is free, never touching a live build or a
  `ready` index. Index tasks use an atomic `queued → building` claim, so the
  two instances can never build the same `(lane, bucket)` index twice, and a
  failed build leaves no usable half-product behind.

### Example

```bash
# 1. request the index (202 = queued/building)
curl -X POST http://localhost:8001/lanes/demoLane/samples/sampleA/index
# 2. poll until "status": "ready"
curl http://localhost:8001/lanes/demoLane/samples/sampleA/index
# 3. fetch every pair carrying the suspicious name (works on either instance)
curl "http://localhost:8002/lanes/demoLane/samples/sampleA/reads?read_name=read7"
```

The lookup response contains `matches: [{ordinal, r1:{name,sequence,plus,
quality}, r2:{...}}, ...]` in original order, plus `total_matches`. An
unknown name is a successful empty result (`total_matches: 0`).

## Concurrency and crash safety

* **SQLite WAL** with `BEGIN IMMEDIATE` job claiming: exactly one instance
  ever picks up a given job; the other keeps serving requests. Read-name
  index tasks use the same pattern with an atomic `queued → building`
  transition, so two instances never build the same `(lane, bucket)` index.
* Each instance processes one task at a time in a background thread, so a
  slow lane or a slow index build never blocks uploads, queries or downloads
  on the event loop. The worker alternates between the job and index queues
  so neither workload can starve the other.
* A job runs under an exclusive `flock` on its work directory. Startup
  recovery removes scratch/staging debris only when no lock is held, and
  never deletes directories that have a published-lane DB row — so another
  instance's in-flight work and all published data are preserved.
* An index build holds an exclusive `flock` on its hidden
  `.build.<id>/lock`. On restart, interrupted builds are marked failed and
  their staging directories removed only when that lock is free; `ready`
  indexes (and any live, locked build) are never deleted.
* Publish order: write/verify outputs in `work/`, move to
  `published/lanes/.tmp.<job>`, atomic `rename` to the final lane directory,
  then insert metadata in one DB transaction. Inputs are deleted only after
  the commit, so a crash at any point leaves no half-visible lane.

## Memory

All I/O is streamed record-by-record (1 MiB network chunks, line-based gzip
reads, one open gzip writer per bucket). No file is ever buffered whole.
The automated streaming test processes hundreds of MiB of decompressed
data with RSS growth below 64 MiB.

Index builds follow the same profile: one record pair is in flight at a
time, records are appended to plain stores, and the SQLite index commits in
1024-record batches. A 200k-pair bucket builds with RSS growth below 64 MiB
(see `test_index_build_constant_memory_on_large_bucket`), comfortably inside
the 128 MiB container limit alongside the running API.

## Layout

```
app/            FastAPI app, SQLite layer, streaming splitter, worker,
                on-disk read-name indexer (app/indexing.py)
scripts/        demo FASTQ generator, two-instance integration check
tests/          pytest suite (happy path, corruption, retry, recovery,
                barcode rules, size budget, streaming memory, index build,
                lookup, limits, index crash/concurrency recovery)
Dockerfile      python:3.11-slim, non-root user, /data volume, healthcheck
docker-compose.yml  two 128 MiB instances sharing one named volume
```

On the shared volume, indexes live under `indexes/`: a hidden
`.build.<index_id>/` staging directory (with a `lock`, partial stores and
SQLite) becomes `indexes/<index_id>/` containing `r1.store`, `r2.store`,
`index.sqlite`, `META`, the held `lock` and the `COMPLETE` marker only once
durable.

## Tests

```bash
python -m pytest tests/ -q
python scripts/integration_two_instances.py   # real two-process check
```

The integration script additionally exercises index build idempotency,
cross-instance lookup, four-way concurrent POSTs collapsing onto a single
shared index task, and continued lane submission/download while indexes
build. Its SIGKILL section runs only where process signals are permitted;
the equivalent index crash recovery is covered deterministically by the
pytest suite elsewhere.
