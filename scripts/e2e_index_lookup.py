#!/usr/bin/env python3
"""Standalone two-process E2E: /docs, index build, cross-instance lookup."""
import json
import os
import shutil
import socket
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from scripts.make_demo_fastq import main as make_demo  # noqa: E402


def free_port():
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    p = s.getsockname()[1]
    s.close()
    return p


def req(method, url, data=None, ct="application/json"):
    r = urllib.request.Request(
        url, data=data, method=method,
        headers={"Content-Type": ct} if data else {})
    try:
        with urllib.request.urlopen(r, timeout=20) as x:
            return x.status, x.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read()


def main():
    pa, pb = free_port(), free_port()
    d = Path("/tmp/demux-e2e")
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    make_demo(str(d / "fx"), pairs=200)
    r1, r2 = d / "fx/R1.fastq.gz", d / "fx/R2.fastq.gz"
    env = {**os.environ, "DEMUX_DATA_DIR": str(d), "DEMUX_POLL_INTERVAL": "0.2"}
    procs = [
        subprocess.Popen(
            [sys.executable, "-m", "uvicorn", "app.main:app",
             "--host", "127.0.0.1", "--port", str(p)],
            cwd=str(ROOT), env=env,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for p in (pa, pb)]
    try:
        for p in (pa, pb):
            deadline = time.time() + 20
            while time.time() < deadline:
                try:
                    if req("GET", f"http://127.0.0.1:{p}/health")[0] == 200:
                        break
                except OSError:
                    pass
                time.sleep(0.3)
            else:
                raise RuntimeError(f"instance on {p} never became healthy")

        for p in (pa, pb):
            assert req("GET", f"http://127.0.0.1:{p}/docs")[0] == 200
            spec = json.loads(req(
                "GET", f"http://127.0.0.1:{p}/openapi.json")[1])
            idx_paths = {x for x in spec["paths"]
                         if "index" in x or x.endswith("/reads")}
            expected = {
                "/lanes/{lane}/samples/{sample_id}/index",
                "/lanes/{lane}/quarantine/index",
                "/lanes/{lane}/indexes",
                "/lanes/{lane}/samples/{sample_id}/reads",
                "/lanes/{lane}/quarantine/reads",
            }
            assert expected <= idx_paths, idx_paths
        print("docs/openapi OK; all index paths documented")

        boundary, parts = "----b", []

        def field(n, v):
            parts.append(
                f"--{boundary}\r\nContent-Disposition: form-data; "
                f'name="{n}"\r\n\r\n{v}\r\n'.encode())

        def ff(n, path):
            parts.append(
                (f"--{boundary}\r\nContent-Disposition: form-data; "
                 f'name="{n}"; filename="{path.name}"\r\n'
                 "Content-Type: application/gzip\r\n\r\n").encode())
            parts.append(path.read_bytes())
            parts.append(b"\r\n")

        bc = [{"sample_id": "sampleA", "barcode_r1": "AAAAAAAA",
               "barcode_r2": "CCCCCCCC"}]
        field("lane", "e2eLane")
        field("barcodes", json.dumps(bc))
        ff("r1", r1)
        ff("r2", r2)
        parts.append(f"--{boundary}--\r\n".encode())
        code, resp = req(
            "POST", f"http://127.0.0.1:{pa}/jobs", b"".join(parts),
            f"multipart/form-data; boundary={boundary}")
        assert code == 201, (code, resp)
        jid = json.loads(resp)["job"]["job_id"]
        deadline = time.time() + 60
        while time.time() < deadline:
            row = json.loads(req(
                "GET", f"http://127.0.0.1:{pa}/jobs/{jid}")[1])
            if row["status"] in ("published", "failed"):
                break
            time.sleep(0.3)
        assert row["status"] == "published", row
        print("lane published")

        code, resp = req(
            "POST", f"http://127.0.0.1:{pa}/lanes/e2eLane/samples/sampleA/index")
        assert code == 202, (code, resp)
        deadline = time.time() + 60
        while time.time() < deadline:
            st = json.loads(req(
                "GET",
                f"http://127.0.0.1:{pb}/lanes/e2eLane/samples/sampleA/index")[1])
            if st["status"] in ("ready", "failed"):
                break
            time.sleep(0.3)
        assert st["status"] == "ready", st
        out = json.loads(req(
            "GET",
            f"http://127.0.0.1:{pb}/lanes/e2eLane/samples/sampleA/reads"
            "?read_name=read0")[1])
        assert out["total_matches"] == 1
        assert out["matches"][0]["r1"]["name"] == "@read0 mate1"
        print("cross-instance lookup OK:",
              out["matches"][0]["r1"]["name"],
              "->", out["matches"][0]["r1"]["sequence"][:12])

        idirs = [x for x in (d / "indexes").iterdir()
                 if x.is_dir() and not x.name.startswith(".build")]
        assert len(idirs) == 1
        need = {"r1.store", "r2.store", "index.sqlite", "META", "COMPLETE"}
        assert need <= {q.name for q in idirs[0].iterdir()}
        print("on-disk layout OK:",
              sorted(q.name for q in idirs[0].iterdir()))
        print("E2E PASSED")
    finally:
        for pr in procs:
            try:
                pr.terminate()
                pr.wait(timeout=15)
            except PermissionError:
                pass
            except Exception:
                try:
                    pr.kill()
                except PermissionError:
                    pass


if __name__ == "__main__":
    main()
