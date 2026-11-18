"""On-demand read-name index: build, lookup, limits, failure and recovery."""
import gzip
import io
import json
import resource
import time
from pathlib import Path

import pytest

BARCODES = [
    {"sample_id": "s1", "barcode_r1": "AAAAAAAA", "barcode_r2": "CCCCCCCC"},
    {"sample_id": "s2", "barcode_r1": "GGGGGGGG", "barcode_r2": "TTTTTTTT"},
]


def _post(client, lane, barcodes, r1_path, r2_path):
    with open(r1_path, "rb") as f1, open(r2_path, "rb") as f2:
        return client.post(
            "/jobs",
            data={"lane": lane, "barcodes": json.dumps(barcodes)},
            files={"r1": ("R1.fastq.gz", f1, "application/gzip"),
                   "r2": ("R2.fastq.gz", f2, "application/gzip")},
        )


def _wait(client, job_id, timeout=20):
    deadline = time.time() + timeout
    while time.time() < deadline:
        row = client.get(f"/jobs/{job_id}").json()
        if row["status"] in ("published", "failed"):
            return row
        time.sleep(0.05)
    raise AssertionError("timeout waiting for job")


def _wait_index(client, path, timeout=20):
    deadline = time.time() + timeout
    last = None
    while time.time() < deadline:
        r = client.get(path)
        assert r.status_code == 200, r.text
        last = r.json()
        if last["status"] in ("ready", "failed"):
            return last
        time.sleep(0.05)
    raise AssertionError(f"index never finished: {last}")


def _publish_lane(client, helpers, tmp_path, lane, r1_records, r2_records,
                  barcodes=BARCODES):
    p1, p2 = helpers.make_pair(r1_records, r2_records, tmp_path)
    resp = _post(client, lane, barcodes, p1, p2)
    assert resp.status_code == 201, resp.text
    job = _wait(client, resp.json()["job"]["job_id"])
    assert job["status"] == "published", job
    return job


def test_build_and_lookup_returns_pairs_in_order(client, helpers, tmp_path,
                                                  monkeypatch):
    # s1 records, including a duplicated read name (two pairs share 'dup').
    r1 = [
        (b"alpha", b"AAAAAAAATT", b"I" * 10),
        (b"dup", b"AAAAAAAATT", b"I" * 10),
        (b"dup", b"AAAAAAAATT", b"I" * 10),  # same name, later ordinal
        (b"beta", b"AAAAAAAATT", b"I" * 10),
    ]
    r2 = [
        (b"alpha", b"CCCCCCCCGG", b"J" * 10),
        (b"dup", b"CCCCCCCCGG", b"J" * 10),
        (b"dup", b"CCCCCCCCGG", b"J" * 10),
        (b"beta", b"CCCCCCCCGG", b"J" * 10),
    ]
    _publish_lane(client, helpers, tmp_path, "laneIdx", r1, r2)

    # Hold the worker inside the build so the task stays 'building' while we
    # exercise the idempotent POST and the 409 lookup rules deterministically.
    import threading
    from app import indexing
    entered = threading.Event()
    release = threading.Event()
    real_build = indexing.build_index

    def stalled_build(task, instance_id):
        if task["lane"] == "laneIdx":
            entered.set()
            release.wait(timeout=10)
        return real_build(task, instance_id)

    monkeypatch.setattr(indexing, "build_index", stalled_build)
    monkeypatch.setattr("app.worker.indexing.build_index", stalled_build)

    # Lookup before any index exists -> 404.
    r = client.get("/lanes/laneIdx/samples/s1/reads?read_name=alpha")
    assert r.status_code == 404

    r = client.post("/lanes/laneIdx/samples/s1/index")
    assert r.status_code == 202
    task = r.json()
    assert task["status"] == "building" and task["sample_id"] == "s1"
    assert entered.wait(timeout=10)

    # A second POST while building is idempotent: same index, HTTP 200.
    r2resp = client.post("/lanes/laneIdx/samples/s1/index")
    assert r2resp.status_code == 200
    assert r2resp.json()["index_id"] == task["index_id"]

    # Queries are refused (409) until the index is fully on disk.
    r = client.get("/lanes/laneIdx/samples/s1/reads?read_name=alpha")
    assert r.status_code == 409 and "building" in r.json()["detail"].lower()

    release.set()
    ready = _wait_index(client, "/lanes/laneIdx/samples/s1/index")
    assert ready["status"] == "ready"
    assert ready["read_pairs"] == 4

    # Now ready, lookups work.
    r = client.get("/lanes/laneIdx/samples/s1/reads?read_name=alpha")
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["total_matches"] == 1
    m = body["matches"][0]
    assert m["ordinal"] == 0
    assert m["r1"]["name"] == "@alpha" and m["r2"]["name"] == "@alpha"
    # Sample mates are trimmed at publish time (8 barcode bases removed).
    assert m["r1"]["sequence"] == "TT"
    assert m["r2"]["quality"] == "J" * 2

    # Duplicated name: both pairs returned in original order, no rescan.
    r = client.get("/lanes/laneIdx/samples/s1/reads?read_name=dup")
    body = r.json()
    assert body["total_matches"] == 2
    assert [m["ordinal"] for m in body["matches"]] == [1, 2]

    # Unknown name is a successful empty result, not an error.
    r = client.get("/lanes/laneIdx/samples/s1/reads?read_name=ghost")
    assert r.status_code == 200 and r.json()["total_matches"] == 0

    # Lane-level index listing.
    listed = client.get("/lanes/laneIdx/indexes").json()
    assert len(listed) == 1 and listed[0]["index_id"] == task["index_id"]


def test_repeated_lookup_does_not_rescan_fastq(client, helpers, tmp_path):
    r1 = [(b"n%d" % i, b"AAAAAAAATT", b"I" * 10) for i in range(50)]
    r2 = [(b"n%d" % i, b"CCCCCCCCGG", b"J" * 10) for i in range(50)]
    _publish_lane(client, helpers, tmp_path, "laneRepeat", r1, r2)
    client.post("/lanes/laneRepeat/samples/s1/index")
    _wait_index(client, "/lanes/laneRepeat/samples/s1/index")

    from app import indexing
    # The source gzip files must never be opened on the lookup path.
    real_gzip_open = gzip.open
    def boom(*a, **k):  # pragma: no cover - must never be called
        raise AssertionError("lookup rescanned a gzip FASTQ")
    gzip.open = boom
    try:
        for i in range(5):
            body = client.get(
                "/lanes/laneRepeat/samples/s1/reads?read_name=n7").json()
            assert body["total_matches"] == 1
            assert body["matches"][0]["ordinal"] == 7
    finally:
        gzip.open = real_gzip_open


def test_quarantine_index_and_lookup(client, helpers, tmp_path):
    r1 = [(b"q1", b"NNNNNNNNTT", b"I" * 10)]
    r2 = [(b"q1", b"ACGTACGTGG", b"J" * 10)]
    _publish_lane(client, helpers, tmp_path, "laneQ", r1, r2)

    assert client.post("/lanes/laneQ/quarantine/index").status_code == 202
    ready = _wait_index(client, "/lanes/laneQ/quarantine/index")
    assert ready["status"] == "ready" and ready["quarantine"] is True
    body = client.get(
        "/lanes/laneQ/quarantine/reads?read_name=q1").json()
    assert body["total_matches"] == 1
    assert body["bucket"] == "__quarantine__"
    # Quarantine records are stored intact (no barcode trimming).
    assert body["matches"][0]["r1"]["sequence"] == "NNNNNNNNTT"


def test_index_on_unknown_lane_or_bucket_is_404(client):
    assert client.post("/lanes/nope/samples/s1/index").status_code == 404
    assert client.post("/lanes/nope/quarantine/index").status_code == 404


def test_oversized_read_name_fails_build(client, helpers, tmp_path, data_dir):
    long_name = b"x" * 257  # read name token > 256 bytes
    r1 = [(long_name, b"AAAAAAAATT", b"I" * 10)]
    r2 = [(long_name, b"CCCCCCCCGG", b"J" * 10)]
    _publish_lane(client, helpers, tmp_path, "laneBigName", r1, r2)

    client.post("/lanes/laneBigName/samples/s1/index")
    task = _wait_index(client, "/lanes/laneBigName/samples/s1/index")
    assert task["status"] == "failed"
    assert "read name" in task["error"]

    # No usable half-product remains on disk (only hidden build dirs would
    # exist transiently; a failed build leaves neither final nor build dirs).
    entries = list(data_dir.joinpath("indexes").iterdir())
    assert entries == []

    # Lookup on a failed index -> 409, and a rebuild request requeues it.
    r = client.get("/lanes/laneBigName/samples/s1/reads?read_name=anything")
    assert r.status_code == 409
    r = client.post("/lanes/laneBigName/samples/s1/index")
    assert r.status_code == 202 and r.json()["status"] == "building"


def test_oversized_sequence_fails_and_boundary_accepted(
        client, helpers, tmp_path):
    # Sample files have 8 barcode bases trimmed at publish time, so the
    # indexed mate is 8 bases shorter than the submission: submit 1009 to
    # produce a 1001-byte published sequence line (over the limit).
    big = 1009
    r1 = [(b"big", b"A" * big, b"I" * big)]
    r2 = [(b"big", b"C" * (big - 1), b"J" * (big - 1))]
    _publish_lane(client, helpers, tmp_path, "laneBigSeq", r1, r2)
    client.post("/lanes/laneBigSeq/samples/s1/index")
    task = _wait_index(client, "/lanes/laneBigSeq/samples/s1/index")
    assert task["status"] == "failed"
    assert "1000" in task["error"] and (
        "sequence" in task["error"] or "quality" in task["error"])

    # Published sequence line exactly at the 1000-byte boundary is accepted.
    d2 = tmp_path / "d2"; d2.mkdir()
    ok = 1008
    r1b = [(b"ok", b"A" * ok, b"I" * ok)]
    r2b = [(b"ok", b"C" * ok, b"J" * ok)]
    _publish_lane(client, helpers, d2, "laneBoundary", r1b, r2b)
    client.post("/lanes/laneBoundary/samples/s1/index")
    ready = _wait_index(client, "/lanes/laneBoundary/samples/s1/index")
    assert ready["status"] == "ready", ready
    body = client.get(
        "/lanes/laneBoundary/samples/s1/reads?read_name=ok").json()
    assert body["matches"][0]["r1"]["sequence"] == "A" * 1000


def test_unknown_bucket_on_published_lane_is_404(client, helpers, tmp_path):
    r1 = [(b"a", b"AAAAAAAATT", b"I" * 10)]
    r2 = [(b"a", b"CCCCCCCCGG", b"J" * 10)]
    _publish_lane(client, helpers, tmp_path, "laneBuckets", r1, r2)
    assert client.post(
        "/lanes/laneBuckets/samples/ghost/index").status_code == 404
    assert client.get(
        "/lanes/laneBuckets/samples/ghost/index").status_code == 404
    assert client.get(
        "/lanes/laneBuckets/samples/ghost/reads?read_name=a").status_code == 404


def test_r1_r2_name_order_mismatch_fails_index(client, helpers, tmp_path,
                                               data_dir):
    # Publish a good lane, then tamper with one published mate so the index
    # build detects an R1/R2 name mismatch (whole index must fail).
    r1 = [(b"a", b"AAAAAAAATT", b"I" * 10), (b"b", b"AAAAAAAATT", b"I" * 10)]
    r2 = [(b"a", b"CCCCCCCCGG", b"J" * 10), (b"b", b"CCCCCCCCGG", b"J" * 10)]
    _publish_lane(client, helpers, tmp_path, "laneTamper", r1, r2)

    from app import db as dbmod
    row = dbmod.get_conn().execute(
        "SELECT r1_file FROM lane_counts WHERE lane='laneTamper' "
        "AND sample_id='s1'").fetchone()
    target = Path(row["r1_file"])
    buf = io.BytesIO()
    with gzip.open(target, "rb") as g:
        buf.write(g.read())
    tampered = buf.getvalue().replace(b"@b\n", b"@z\n", 1)
    with gzip.open(target, "wb") as g:
        g.write(tampered)

    client.post("/lanes/laneTamper/samples/s1/index")
    task = _wait_index(client, "/lanes/laneTamper/samples/s1/index")
    assert task["status"] == "failed" and "names differ" in task["error"]
    # No half-built index is queryable.
    assert client.get(
        "/lanes/laneTamper/samples/s1/reads?read_name=a").status_code == 409


def test_published_count_mismatch_fails_index(client, helpers, tmp_path,
                                              monkeypatch):
    r1 = [(b"a", b"AAAAAAAATT", b"I" * 10), (b"b", b"AAAAAAAATT", b"I" * 10)]
    r2 = [(b"a", b"CCCCCCCCGG", b"J" * 10), (b"b", b"CCCCCCCCGG", b"J" * 10)]
    _publish_lane(client, helpers, tmp_path, "laneCount", r1, r2)

    # Lie about the published count: indexed pairs must then disagree.
    from app import db as dbmod
    dbmod.get_conn().execute(
        "UPDATE lane_counts SET read_pairs=99 WHERE lane='laneCount' "
        "AND sample_id='s1'")
    client.post("/lanes/laneCount/samples/s1/index")
    task = _wait_index(client, "/lanes/laneCount/samples/s1/index")
    assert task["status"] == "failed" and "count table" in task["error"]


def test_missing_published_file_fails_index(client, helpers, tmp_path):
    r1 = [(b"a", b"AAAAAAAATT", b"I" * 10)]
    r2 = [(b"a", b"CCCCCCCCGG", b"J" * 10)]
    _publish_lane(client, helpers, tmp_path, "laneMissing", r1, r2)
    from app import db as dbmod
    row = dbmod.get_conn().execute(
        "SELECT r1_file FROM lane_counts WHERE lane='laneMissing' "
        "AND sample_id='s1'").fetchone()
    Path(row["r1_file"]).unlink()
    client.post("/lanes/laneMissing/samples/s1/index")
    task = _wait_index(client, "/lanes/laneMissing/samples/s1/index")
    assert task["status"] == "failed" and "missing" in task["error"].lower()


def test_lookup_name_validation(client, helpers, tmp_path):
    r1 = [(b"a", b"AAAAAAAATT", b"I" * 10)]
    r2 = [(b"a", b"CCCCCCCCGG", b"J" * 10)]
    _publish_lane(client, helpers, tmp_path, "laneVal", r1, r2)
    client.post("/lanes/laneVal/samples/s1/index")
    _wait_index(client, "/lanes/laneVal/samples/s1/index")
    # Missing parameter -> 422 from FastAPI; overlong value -> 400.
    assert client.get(
        "/lanes/laneVal/samples/s1/reads").status_code == 422
    assert client.get(
        "/lanes/laneVal/samples/s1/reads?read_name=" + "x" * 257
    ).status_code == 400


def test_recovery_clears_interrupted_index_build(data_dir):
    import fcntl
    from app import db as dbmod, indexing

    dbmod.get_conn().execute(
        "INSERT INTO lane_publishes (lane, job_id, published_at, seq) "
        "VALUES ('recLane','j1',?,1)", (dbmod.utcnow(),))
    dbmod.get_conn().execute(
        "INSERT INTO lane_counts (lane, sample_id, r1_file, r2_file, "
        "read_pairs) VALUES ('recLane','s1','/tmp/a','/tmp/b',0)")
    index_id = "idxdead"
    dbmod.create_index_task(
        {"index_id": index_id, "lane": "recLane", "bucket": "s1",
         "created_at": dbmod.utcnow()})
    dbmod.get_conn().execute(
        "UPDATE index_tasks SET status='building', instance_id='ghost' "
        "WHERE index_id=?", (index_id,))

    bd = data_dir / "indexes" / f".build.{index_id}"
    bd.mkdir(parents=True)
    (bd / "index.sqlite").write_bytes(b"half-built")
    (bd / "lock").write_text("ghost\n")

    notes = indexing.recover_indexes_on_startup()
    assert any(index_id in n for n in notes)
    assert not bd.exists()
    row = dbmod.get_index_task(index_id)
    assert row["status"] == "failed" and "interrupted" in row["error"]
    # Idempotent second run.
    assert indexing.recover_indexes_on_startup() == []


def test_recovery_never_removes_live_index_build(data_dir):
    import fcntl
    from app import db as dbmod, indexing

    index_id = "idxlive"
    dbmod.get_conn().execute(
        "INSERT INTO lane_publishes (lane, job_id, published_at, seq) "
        "VALUES ('liveLane','j1',?,1)", (dbmod.utcnow(),))
    dbmod.get_conn().execute(
        "INSERT INTO lane_counts (lane, sample_id, r1_file, r2_file, "
        "read_pairs) VALUES ('liveLane','s1','/tmp/a','/tmp/b',0)")
    dbmod.create_index_task(
        {"index_id": index_id, "lane": "liveLane", "bucket": "s1",
         "created_at": dbmod.utcnow()})
    dbmod.get_conn().execute(
        "UPDATE index_tasks SET status='building', instance_id='live' "
        "WHERE index_id=?", (index_id,))
    bd = data_dir / "indexes" / f".build.{index_id}"
    bd.mkdir(parents=True)
    lock_fh = open(bd / "lock", "w")
    fcntl.flock(lock_fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    try:
        marker = bd / "inflight"
        marker.write_text("x")
        notes = indexing.recover_indexes_on_startup()
        assert marker.exists()  # flock owner protects the live build
        assert not any(index_id in n for n in notes)
        assert dbmod.get_index_task(index_id)["status"] == "building"
    finally:
        fcntl.flock(lock_fh.fileno(), fcntl.LOCK_UN)
        lock_fh.close()


def test_recovery_reconciles_completed_index_with_stale_status(data_dir):
    """A durable COMPLETE directory whose row is still 'building' (crash
    between rename and status commit) must be reconciled to ready, not
    deleted or left unqueryable."""
    from app import db as dbmod, indexing

    dbmod.get_conn().execute(
        "INSERT INTO lane_publishes (lane, job_id, published_at, seq) "
        "VALUES ('cLane','j1',?,1)", (dbmod.utcnow(),))
    dbmod.get_conn().execute(
        "INSERT INTO lane_counts (lane, sample_id, r1_file, r2_file, "
        "read_pairs) VALUES ('cLane','s1','/tmp/a','/tmp/b',0)")
    index_id = "idxdone"
    dbmod.create_index_task(
        {"index_id": index_id, "lane": "cLane", "bucket": "s1",
         "created_at": dbmod.utcnow()})
    dbmod.get_conn().execute(
        "UPDATE index_tasks SET status='building', instance_id='ghost' "
        "WHERE index_id=?", (index_id,))
    final = data_dir / "indexes" / index_id
    final.mkdir(parents=True)
    (final / "COMPLETE").touch()
    (final / "r1.store").write_text("@a\nA\n+\nI\n")
    (final / "r2.store").write_text("@a\nC\n+\nJ\n")
    iconn = indexing._open_index_db(final)
    indexing._flush_batch(iconn, [(0, "a", 0, 0)])
    iconn.close()

    notes = indexing.recover_indexes_on_startup()
    assert any("reconciling completed" in n for n in notes)
    row = dbmod.get_index_task(index_id)
    assert row["status"] == "ready"
    assert (final / "COMPLETE").exists()  # usable index preserved


def test_recovery_repairs_ready_row_without_on_disk_index(data_dir):
    from app import db as dbmod, indexing

    dbmod.get_conn().execute(
        "INSERT INTO lane_publishes (lane, job_id, published_at, seq) "
        "VALUES ('rLane','j1',?,1)", (dbmod.utcnow(),))
    dbmod.get_conn().execute(
        "INSERT INTO lane_counts (lane, sample_id, r1_file, r2_file, "
        "read_pairs) VALUES ('rLane','s1','/tmp/a','/tmp/b',0)")
    index_id = "idxorphan"
    dbmod.create_index_task(
        {"index_id": index_id, "lane": "rLane", "bucket": "s1",
         "created_at": dbmod.utcnow()})
    dbmod.get_conn().execute(
        "UPDATE index_tasks SET status='ready' WHERE index_id=?",
        (index_id,))
    notes = indexing.recover_indexes_on_startup()
    assert any(index_id in n for n in notes)
    assert dbmod.get_index_task(index_id)["status"] == "failed"


def test_index_build_constant_memory_on_large_bucket(
        client, helpers, tmp_path, data_dir):
    pairs = 200_000  # ~16 MiB decompressed per mate, many distinct names
    r1 = [(b"read%07d" % i, b"AAAAAAAATT", b"I" * 10) for i in range(pairs)]
    r2 = [(b"read%07d" % i, b"CCCCCCCCGG", b"J" * 10) for i in range(pairs)]
    _publish_lane(client, helpers, tmp_path, "laneMem", r1, r2)

    rss_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    assert client.post("/lanes/laneMem/samples/s1/index").status_code == 202
    ready = _wait_index(client, "/lanes/laneMem/samples/s1/index",
                        timeout=120)
    assert ready["status"] == "ready", ready
    rss_after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    growth_mib = (rss_after - rss_before) / 1024
    print(f"\nindexed {pairs} pairs, RSS growth {growth_mib:.1f} MiB")
    assert growth_mib < 64, f"RSS grew {growth_mib:.1f} MiB during indexing"

    # Spot lookups at both ends of the ordinal space.
    for i in (0, pairs - 1):
        body = client.get(
            f"/lanes/laneMem/samples/s1/reads?read_name=read{i:07d}").json()
        assert body["total_matches"] == 1
        assert body["matches"][0]["ordinal"] == i
