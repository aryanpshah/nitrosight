"""Background threads, so lookups and image recognition never freeze the window."""

from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import QThread, pyqtSignal

from osr import OsrService, Recognition
from pipeline import ScanResult, scan
from resolve import FDALimits, PubChemClient


class ScanWorker(QThread):
    """Assess one compound."""

    done = pyqtSignal(object)      # ScanResult
    failed = pyqtSignal(str)

    def __init__(self, query: str, kind: str = "auto",
                 client: PubChemClient | None = None,
                 limits: FDALimits | None = None,
                 ruleset=None):
        super().__init__()
        self._query = query
        self._kind = kind
        self._client = client
        self._limits = limits
        self._ruleset = ruleset

    def run(self) -> None:
        try:
            self.done.emit(scan(self._query, self._client, self._limits, self._kind,
                                self._ruleset))
        except Exception as e:  # noqa: BLE001
            self.failed.emit(f"{type(e).__name__}: {e}")


class RecogniseWorker(QThread):
    """Read a structure from a picture. Slow the first time, because the recognition
    model has to load; quick afterwards while the service stays resident."""

    status = pyqtSignal(str)
    done = pyqtSignal(object)      # Recognition
    failed = pyqtSignal(str)

    def __init__(self, image: Path, service: OsrService):
        super().__init__()
        self._image = image
        self._service = service

    def run(self) -> None:
        try:
            if not self._service.ready:
                self.status.emit(
                    "Loading the structure recognition model, this takes about a "
                    "minute the first time...")
            result = self._service.recognise(self._image)
            self.done.emit(result)
        except Exception as e:  # noqa: BLE001
            self.failed.emit(str(e))


class BatchWorker(QThread):
    """Assess a list of compounds, reporting progress as it goes."""

    progress = pyqtSignal(int, int, str)       # done, total, current label
    row_done = pyqtSignal(int, object)         # index, ScanResult
    finished_all = pyqtSignal(int)
    failed = pyqtSignal(str)

    def __init__(self, queries: list[str], kind: str = "auto", ruleset=None):
        super().__init__()
        self._queries = queries
        self._kind = kind
        self._ruleset = ruleset
        self._cancelled = False

    def cancel(self) -> None:
        self._cancelled = True

    def run(self) -> None:
        try:
            client = PubChemClient(verbose=False)
            limits = FDALimits()
            total = len(self._queries)
            completed = 0

            for i, query in enumerate(self._queries):
                if self._cancelled:
                    break
                self.progress.emit(completed, total, query[:60])
                try:
                    result = scan(query, client, limits, self._kind, self._ruleset)
                except Exception as e:  # noqa: BLE001
                    result = ScanResult(query=query, error=f"{type(e).__name__}: {e}")
                self.row_done.emit(i, result)
                completed += 1

            client.save()
            self.finished_all.emit(completed)
        except Exception as e:  # noqa: BLE001
            self.failed.emit(f"{type(e).__name__}: {e}")
