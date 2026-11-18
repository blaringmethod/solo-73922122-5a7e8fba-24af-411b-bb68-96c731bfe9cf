"""On-demand, on-disk read-name indexes for published lane buckets.

A bucket is one per-sample pair of files or the quarantine pair. Building an
index streams the *published* gzip FASTQ mates record by record, verifies
FASTQ integrity / R1-R2 name order / published counts, and materialises:

    <INDEX_DIR>/<index_id>/
        r1.store, r2.store   plain FASTQ records (four lines each)
        index.sqlite         ordinal + read name + byte offsets per record
        META                 build metadata
        COMPLETE            durable marker, written last

Only the presence of COMPLETE *and* the atomic rename into the final
directory make an index queryable: builds always happen in a hidden
``.build.<id>`` directory, and any failure removes it entirely, so no
half-built index is ever visible. Repeated lookups never rescan the FASTQ:
the SQLite index maps a name straight to store offsets.

Memory profile is constant: one record pair in flight at a time, batched
SQLite commits, no whole-file buffers.
"""
from __future__ import annotations

import errno
import fcntl
import gzip
import json
import os
import shutil
import sqlite3
import time
import zlib
from dataclasses import dataclass
from pathlib import Path

from . import config, db

QUARANTINE_BUCKET = ""  # stored in index_tasks.bucket for the quarantine pair


class IndexBuildError(Exception):
    """Raised when an index cannot be built; the task fails as a whole."""


class IndexNotReadyError(Exception):
    """Raised on lookup while the index is absent, building or failed."""

    def __init__(self, message: str, status: str):
        super().__init__(message)
        self.status = status


@dataclass
class _MateRecord:
    header: bytes  # full header line without newline, including '@'
    name_token: bytes  # first token after '@'
    seq: bytes
    plus: bytes
    qual: bytes


def _build_dir(index_id: str) -> Path:
    return config.INDEX_DIR / f".build.{index_id}"


def final_dir(index_id: str) -> Path:
    return config.INDEX_DIR / index_id


def _lock_path(index_id: str) -> Path:
    return _build_dir(index_id) / "lock"


def index_db_path(directory: Path) -> Path:
    return directory / "index.sqlite"


def _lock_is_free(path: Path) -> bool:
    """True iff nobody currently holds flock on path."""
    if not path.exists():
        return True
    fh = open(path, "r")
    try:
        try:
            fcntl.flock(fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            if exc.errno in (errno.EACCES, errno.EAGAIN):
                return False
            raise
        fcntl.flock(fh.fileno(), fcntl.LOCK_UN)
        return True
    finally:
        fh.close()


def _acquire_build_lock(index_id: str):
    bd = _build_dir(index_id)
    bd.mkdir(parents=True, exist_ok=True)
    fh = open(_lock_path(index_id), "w")
    try:
        fcntl.flock(fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as exc:
        fh.close()
        if exc.errno in (errno.EACCES, errno.EAGAIN):
            raise RuntimeError(
                f"index {index_id} is being built by another instance") from exc
        raise
    fh.write(f"{os.getpid()}\n")
    fh.flush()
    return fh


def _read_mate_record(gf: gzip.GzipFile, mate: str,
                      ordinal: int) -> _MateRecord | None:
    """Read and limit-check one record from one mate. None at clean EOF."""
    header_line = gf.readline()
    if header_line == b"":
        return None
    seq_line = gf.readline()
    plus_line = gf.readline()
    qual_line = gf.readline()
    if not header_line.startswith(b"@") or seq_line == b"" \
            or not plus_line.startswith(b"+") or qual_line == b"":
        raise IndexBuildError(
            f"{mate}: record {ordinal} is not a well-formed FASTQ record")
    header = header_line.rstrip(b"\r\n")
    seq = seq_line.rstrip(b"\r\n")
    plus = plus_line.rstrip(b"\r\n")
    qual = qual_line.rstrip(b"\r\n")
    if not seq or seq_line.count(b"\n") == 0:
        raise IndexBuildError(f"{mate}: record {ordinal} is truncated")
    parts = header[1:].split()
    if not parts:
        raise IndexBuildError(f"{mate}: record {ordinal} has an empty read name")
    name_token = parts[0]
    if len(name_token) > config.MAX_READ_NAME_BYTES:
        raise IndexBuildError(
            f"{mate}: record {ordinal} read name is {len(name_token)} bytes, "
            f"limit is {config.MAX_READ_NAME_BYTES}")
    if len(seq) > config.MAX_SEQ_QUAL_BYTES:
        raise IndexBuildError(
            f"{mate}: record {ordinal} sequence line is {len(seq)} bytes, "
            f"limit is {config.MAX_SEQ_QUAL_BYTES}")
    if len(qual) > config.MAX_SEQ_QUAL_BYTES:
        raise IndexBuildError(
            f"{mate}: record {ordinal} quality line is {len(qual)} bytes, "
            f"limit is {config.MAX_SEQ_QUAL_BYTES}")
    if len(seq) != len(qual):
        raise IndexBuildError(
            f"{mate}: record {ordinal} sequence/quality length mismatch")
    return _MateRecord(header, name_token, seq, plus, qual)


def _open_index_db(directory: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(index_db_path(directory), timeout=30,
                           isolation_level=None)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA busy_timeout=30000")
    conn.execute(
        "CREATE TABLE reads ("
        "ord INTEGER PRIMARY KEY, "
        "name TEXT NOT NULL, "
        "r1_off INTEGER NOT NULL, "
        "r2_off INTEGER NOT NULL)")
    conn.execute("CREATE INDEX idx_reads_name ON reads(name)")
    return conn


def build_index(task_row, instance_id: str) -> None:
    """Build one bucket index end to end. Any error fails the whole task and
    leaves no usable index directory behind."""
    index_id = task_row["index_id"]
    # Defensive: only the claiming owner may build, so a stale claim can never
    # overwrite another instance's live build.
    if task_row["status"] != "building" or task_row["instance_id"] != instance_id:
        raise RuntimeError(
            f"instance {instance_id} does not own index {index_id}")
    lane = task_row["lane"]
    bucket = task_row["bucket"]
    count_row = db.get_conn().execute(
        "SELECT * FROM lane_counts WHERE lane=? AND sample_id=?",
        (lane, bucket)).fetchone()
    if count_row is None:
        raise IndexBuildError(
            f"published bucket {bucket or 'quarantine'} vanished from lane {lane}")
    r1_path = Path(count_row["r1_file"])
    r2_path = Path(count_row["r2_file"])
    expected_pairs = count_row["read_pairs"]
    if not r1_path.is_file() or not r2_path.is_file():
        raise IndexBuildError("published R1/R2 file is missing")

    bd = _build_dir(index_id)
    final = final_dir(index_id)
    lock = _acquire_build_lock(index_id)
    try:
        # Debris from a previous crashed attempt of this same task. We hold
        # the exclusive flock on bd/lock, so remove everything *except* the
        # lock file itself: unlinking it would let a sibling acquire a fresh
        # lock on a new inode while this build is still running.
        for entry in bd.iterdir():
            if entry.name == "lock":
                continue
            if entry.is_dir():
                shutil.rmtree(entry, ignore_errors=True)
            else:
                entry.unlink(missing_ok=True)
        if final.exists():
            shutil.rmtree(final, ignore_errors=True)

        store1 = open(bd / "r1.store", "wb")
        store2 = open(bd / "r2.store", "wb")
        iconn = _open_index_db(bd)
        pending: list[tuple[int, str, int, int]] = []
        pairs = 0
        try:
            with open(r1_path, "rb") as f1, open(r2_path, "rb") as f2:
                try:
                    with gzip.GzipFile(fileobj=f1) as g1, \
                            gzip.GzipFile(fileobj=f2) as g2:
                        while True:
                            try:
                                m1 = _read_mate_record(g1, "R1", pairs)
                                m2 = _read_mate_record(g2, "R2", pairs)
                            except (OSError, EOFError, zlib.error) as exc:
                                raise IndexBuildError(
                                    f"gzip stream corrupted: {exc}") from exc
                            if m1 is None and m2 is None:
                                break
                            if m1 is None or m2 is None:
                                raise IndexBuildError(
                                    "R1 and R2 contain different numbers of "
                                    "records")
                            if m1.name_token != m2.name_token:
                                raise IndexBuildError(
                                    f"record {pairs}: R1/R2 read names differ "
                                    f"in order "
                                    f"({m1.name_token.decode('ascii', 'replace')}"
                                    f" vs "
                                    f"{m2.name_token.decode('ascii', 'replace')})")
                            off1 = store1.tell()
                            store1.write(m1.header + b"\n" + m1.seq + b"\n"
                                         + m1.plus + b"\n" + m1.qual + b"\n")
                            off2 = store2.tell()
                            store2.write(m2.header + b"\n" + m2.seq + b"\n"
                                         + m2.plus + b"\n" + m2.qual + b"\n")
                            pending.append(
                                (pairs, m1.name_token.decode("latin-1"),
                                 off1, off2))
                            pairs += 1
                            if len(pending) >= config.INDEX_COMMIT_EVERY:
                                _flush_batch(iconn, pending)
                                pending.clear()
                                db.index_progress(index_id, instance_id, pairs)
                                db.index_heartbeat(index_id, instance_id)
                except IndexBuildError:
                    raise
                except (OSError, EOFError, zlib.error) as exc:
                    raise IndexBuildError(
                        f"gzip stream corrupted: {exc}") from exc
            if pending:
                _flush_batch(iconn, pending)
                pending.clear()
        finally:
            store1.close()
            store2.close()
            iconn.close()

        # Sequential cross-check against the published count table.
        if pairs != expected_pairs:
            raise IndexBuildError(
                f"indexed {pairs} read pairs but the published count table "
                f"records {expected_pairs}")

        meta = {
            "index_id": index_id,
            "lane": lane,
            "bucket": bucket,
            "read_pairs": pairs,
            "built_at": db.utcnow(),
            "source_r1": str(r1_path),
            "source_r2": str(r2_path),
            "limits": {
                "max_read_name_bytes": config.MAX_READ_NAME_BYTES,
                "max_seq_qual_bytes": config.MAX_SEQ_QUAL_BYTES,
            },
        }
        (bd / "META").write_text(json.dumps(meta, indent=2))
        # fsync data files before the completion marker / rename.
        for name in ("r1.store", "r2.store", "index.sqlite"):
            fd = os.open(bd / name, os.O_RDONLY)
            try:
                os.fsync(fd)
            finally:
                os.close(fd)
        (bd / "COMPLETE").touch()
        os.rename(bd, final)
        try:
            db.mark_index_ready(index_id, instance_id, pairs)
        except Exception:
            # The on-disk index is complete (COMPLETE present) even though the
            # status row could not be updated. Leave the durable directory in
            # place: startup recovery reconciles the row against COMPLETE
            # instead of deleting a perfectly usable index.
            raise
    except Exception:
        # Remove only the hidden staging directory. A final directory carries
        # a COMPLETE marker and must never be torn down as a "half-product".
        shutil.rmtree(bd, ignore_errors=True)
        raise
    finally:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        finally:
            lock.close()


def _flush_batch(iconn: sqlite3.Connection,
                 rows: list[tuple[int, str, int, int]]) -> None:
    iconn.execute("BEGIN")
    try:
        iconn.executemany(
            "INSERT INTO reads (ord, name, r1_off, r2_off) VALUES (?,?,?,?)",
            rows)
        iconn.execute("COMMIT")
    except Exception:
        iconn.execute("ROLLBACK")
        raise


def _record_at(store: Path, offset: int) -> dict[str, str]:
    with open(store, "rb") as fh:
        fh.seek(offset)
        header = fh.readline().rstrip(b"\r\n")
        seq = fh.readline().rstrip(b"\r\n")
        plus = fh.readline().rstrip(b"\r\n")
        qual = fh.readline().rstrip(b"\r\n")
    return {
        "name": header.decode("latin-1"),
        "sequence": seq.decode("latin-1"),
        "plus": plus.decode("latin-1"),
        "quality": qual.decode("latin-1"),
    }


def lookup_reads(lane: str, bucket: str, name: str) -> dict:
    """Return every pair carrying ``name`` in original ordinal order.

    Never rescans the FASTQ: the on-disk SQLite index answers directly.
    """
    row = db.get_index_task_by_bucket(lane, bucket)
    if row is None:
        raise IndexNotReadyError(
            "no index has been built for this bucket yet", "absent")
    if row["status"] in ("queued", "building"):
        raise IndexNotReadyError(
            "index is still building; queries are allowed only once it is "
            "ready on disk", "building")
    if row["status"] == "failed":
        raise IndexNotReadyError(
            f"index build failed: {row['error'] or 'unknown error'}", "failed")
    directory = final_dir(row["index_id"])
    if not (directory / "COMPLETE").is_file():
        # Durable marker missing despite a ready row: treat as failed so the
        # caller can request a rebuild rather than reading a broken index.
        raise IndexNotReadyError(
            "on-disk index is incomplete; rebuild required", "failed")
    try:
        name_bytes = name.encode("latin-1")
    except UnicodeEncodeError as exc:
        raise IndexBuildError(f"read name cannot be expressed as bytes: {exc}")
    if not name_bytes or len(name_bytes) > config.MAX_READ_NAME_BYTES:
        raise IndexBuildError(
            f"read name must be 1..{config.MAX_READ_NAME_BYTES} bytes")

    conn = sqlite3.connect(index_db_path(directory), timeout=30)
    try:
        hits = conn.execute(
            "SELECT ord, r1_off, r2_off FROM reads WHERE name=? ORDER BY ord",
            (name,)).fetchall()
    finally:
        conn.close()
    matches = [
        {
            "ordinal": ord,
            "r1": _record_at(directory / "r1.store", off1),
            "r2": _record_at(directory / "r2.store", off2),
        }
        for ord, off1, off2 in hits
    ]
    return {
        "lane": lane,
        "bucket": bucket or "__quarantine__",
        "read_name": name,
        "matches": matches,
        "total_matches": len(matches),
    }


# ---------------- startup recovery ----------------

def recover_indexes_on_startup() -> list[str]:
    """Clean interrupted index builds.

    Rules mirror job recovery: flock ownership protects a live build, and a
    final index directory is never removed while its DB row is 'ready'.
    """
    notes: list[str] = []
    conn = db.get_conn()
    now = time.time()

    # 1. Tasks stuck in 'building' whose owner lock is free. Normally the
    #    owner died mid-build -> fail + clean. The exception is a build whose
    #    on-disk directory reached COMPLETE but whose status update did not
    #    commit (crash between rename and mark_ready): that index is durable
    #    and is reconciled back to 'ready' instead of being discarded.
    #    'queued' tasks are not owned yet and simply get claimed normally.
    for task in conn.execute(
            "SELECT * FROM index_tasks WHERE status='building'").fetchall():
        lock = _lock_path(task["index_id"])
        if not _lock_is_free(lock):
            continue  # another instance is actively building this index
        directory = final_dir(task["index_id"])
        if (directory / "COMPLETE").is_file():
            notes.append(
                f"reconciling completed on-disk index {task['index_id']}")
            conn.execute(
                "UPDATE index_tasks SET status='ready', instance_id=NULL,"
                "error=NULL, updated_at=? WHERE index_id=?",
                (db.utcnow(), task["index_id"]))
            shutil.rmtree(_build_dir(task["index_id"]), ignore_errors=True)
            continue
        notes.append(
            f"recovering interrupted index build {task['index_id']}")
        conn.execute(
            "UPDATE index_tasks SET status='failed', instance_id=NULL,"
            "error='interrupted (owner process disappeared)', updated_at=? "
            "WHERE index_id=?", (db.utcnow(), task["index_id"]))
        shutil.rmtree(_build_dir(task["index_id"]), ignore_errors=True)

    # 2. Orphan .build.* directories with no live building task attached.
    building_ids = {
        r["index_id"] for r in conn.execute(
            "SELECT index_id FROM index_tasks WHERE status='building'").fetchall()
    }
    if config.INDEX_DIR.exists():
        for child in config.INDEX_DIR.iterdir():
            if not child.name.startswith(".build."):
                continue
            index_id = child.name[len(".build."):]
            if index_id in building_ids:
                continue  # row exists; step 1 already checked the lock
            lock = child / "lock"
            if lock.exists() and not _lock_is_free(lock):
                continue  # a live instance owns it even without a visible row
            if now - child.stat().st_mtime < 5:
                continue  # just created by a racing sibling
            notes.append(f"removing orphan index build directory {child.name}")
            shutil.rmtree(child, ignore_errors=True)

    # 3. Ready rows whose durable directory vanished -> mark failed so the
    #    bucket can be rebuilt. Never delete a directory that is still valid.
    for task in conn.execute(
            "SELECT * FROM index_tasks WHERE status='ready'").fetchall():
        directory = final_dir(task["index_id"])
        if not (directory / "COMPLETE").is_file():
            notes.append(
                f"index {task['index_id']} lost its on-disk data; marked failed")
            conn.execute(
                "UPDATE index_tasks SET status='failed', instance_id=NULL,"
                "error='on-disk index data missing', updated_at=? "
                "WHERE index_id=?", (db.utcnow(), task["index_id"]))
            shutil.rmtree(directory, ignore_errors=True)

    # 4. Failed tasks must not leave final directories behind.
    for task in conn.execute(
            "SELECT * FROM index_tasks WHERE status='failed'").fetchall():
        leftover = final_dir(task["index_id"])
        if leftover.exists():
            notes.append(
                f"removing leftover directory of failed index {task['index_id']}")
            shutil.rmtree(leftover, ignore_errors=True)
    return notes
