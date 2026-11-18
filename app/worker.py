"""Background worker: one task at a time per instance.

Runs in a daemon thread inside the uvicorn process. While a slow lane or a
slow index build is running the event loop keeps serving uploads, queries
and downloads; the poller simply does not claim another task until the
current one finishes.

Demultiplexing jobs and on-demand index builds live in separate queues. The
worker alternates which queue it probes first, so a long index build cannot
starve lane processing (or vice versa), and both instances share both
queues through SQLite.
"""
from __future__ import annotations

import logging
import threading
import uuid

from . import config, db, indexing, processing

log = logging.getLogger("demux.worker")


def new_instance_id() -> str:
    return f"{uuid.getnode():x}-{uuid.uuid4().hex[:8]}"


class JobWorker(threading.Thread):
    def __init__(self, instance_id: str):
        super().__init__(name=f"demux-worker-{instance_id[:8]}", daemon=True)
        self.instance_id = instance_id
        self._stop = threading.Event()
        self.current_job: str | None = None
        self._jobs_first = True

    def stop(self) -> None:
        self._stop.set()

    def run(self) -> None:
        log.info("worker %s started", self.instance_id)
        while not self._stop.is_set():
            claimed = self._claim_one()
            if claimed is None:
                self._stop.wait(config.POLL_INTERVAL)
                continue
            kind, task = claimed
            self.current_job = (
                f"job:{task['job_id']}" if kind == "job"
                else f"index:{task['index_id']}")
            if kind == "job":
                self._run_job(task)
            else:
                self._run_index(task)
            self.current_job = None

    def _claim_one(self):
        """Alternate queue priority, then fall back to the other queue."""
        order = ("job", "index") if self._jobs_first else ("index", "job")
        self._jobs_first = not self._jobs_first
        for kind in order:
            try:
                if kind == "job":
                    task = db.claim_next_job(self.instance_id)
                else:
                    task = db.claim_next_index_task(self.instance_id)
            except Exception:
                log.exception("failed to claim a %s task", kind)
                task = None
            if task is not None:
                return kind, task
        return None

    def _run_job(self, job) -> None:
        log.info("instance %s claimed job %s lane %s",
                 self.instance_id, job["job_id"], job["lane"])
        try:
            processing.process_job(job, self.instance_id)
            log.info("job %s published", job["job_id"])
        except processing.FastqFormatError as exc:
            log.warning("job %s failed on format: %s", job["job_id"], exc)
            self._fail_job(job, f"FASTQ format error: {exc}")
        except processing.SizeExceededError as exc:
            log.warning("job %s failed on size: %s", job["job_id"], exc)
            self._fail_job(job, str(exc))
        except Exception:
            log.exception("job %s crashed", job["job_id"])
            self._fail_job(job, "internal processing error")

    def _run_index(self, task) -> None:
        log.info("instance %s claimed index task %s lane %s bucket %r",
                 self.instance_id, task["index_id"], task["lane"],
                 task["bucket"] or "quarantine")
        try:
            indexing.build_index(task, self.instance_id)
            log.info("index %s ready", task["index_id"])
        except indexing.IndexBuildError as exc:
            log.warning("index %s failed: %s", task["index_id"], exc)
            self._fail_index(task, str(exc))
        except Exception:
            log.exception("index %s crashed", task["index_id"])
            self._fail_index(task, "internal indexing error")

    def _fail_job(self, job, message: str) -> None:
        try:
            db.mark_failed(job["job_id"], self.instance_id, message)
        except Exception:
            log.exception("could not mark job %s failed", job["job_id"])

    def _fail_index(self, task, message: str) -> None:
        try:
            db.mark_index_failed(task["index_id"], self.instance_id, message)
        except Exception:
            log.exception(
                "could not mark index %s failed", task["index_id"])
