"""
Client for the optical structure recognition service.

Recognition needs TensorFlow, which has no Python 3.14 build, so it runs in its own
environment (.venv-osr) as a separate process. Loading the model takes close to a minute,
so that process is kept resident and answers requests through a small file-based
protocol rather than being restarted per image.

Nothing here scores anything. Recognition turns a picture into a SMILES string and stops;
the CPCA scorer takes over from there. That separation is deliberate: the fallible,
model-driven step never touches the number.
"""

from __future__ import annotations

import json
import os
import subprocess
import time
import uuid
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).parent
VENV_OSR = ROOT / ".venv-osr"


def venv_python(venv: Path = VENV_OSR) -> Path:
    """The interpreter inside a virtual environment, on either layout.

    Windows puts it in Scripts/python.exe; POSIX puts it in bin/python. Returning
    whichever exists keeps image recognition working on Linux and macOS as well as
    Windows. When neither is present the Windows path is returned so that
    is_available() reports False and the setup hint is shown.
    """
    for candidate in (venv / "Scripts" / "python.exe", venv / "bin" / "python"):
        if candidate.exists():
            return candidate
    return venv / "Scripts" / "python.exe"


WORKER_PY = venv_python()
WORKER = ROOT / "tools" / "osr_worker.py"
WORK_DIR = ROOT / "data" / "osr_service"

STARTUP_TIMEOUT = 300.0      # first load pulls in TensorFlow and the model weights
REQUEST_TIMEOUT = 180.0


@dataclass
class Recognition:
    smiles: str | None = None
    error: str | None = None

    @property
    def ok(self) -> bool:
        return bool(self.smiles)


def is_available() -> bool:
    """Is the recognition environment installed?"""
    return WORKER_PY.exists() and WORKER.exists()


if os.name == "nt":
    _SETUP_COMMANDS = ("    py -3.13 -m venv .venv-osr\n"
                       "    .venv-osr\\Scripts\\python -m pip install decimer")
else:
    _SETUP_COMMANDS = ("    python3.13 -m venv .venv-osr\n"
                       "    .venv-osr/bin/python -m pip install decimer")

SETUP_HINT = (
    "Structure recognition is not installed.\n\n"
    "It needs its own environment because TensorFlow has no Python 3.14 build:\n"
    + _SETUP_COMMANDS)


class OsrService:
    """A resident recogniser. Starting it is slow; each image after that is not."""

    def __init__(self, work_dir: Path = WORK_DIR):
        self.work_dir = work_dir
        self._proc: subprocess.Popen | None = None

    @property
    def running(self) -> bool:
        return self._proc is not None and self._proc.poll() is None

    @property
    def ready(self) -> bool:
        return self.running and (self.work_dir / "ready.flag").exists()

    def start(self) -> None:
        if self.running:
            return
        if not is_available():
            raise RuntimeError(SETUP_HINT)

        self.work_dir.mkdir(parents=True, exist_ok=True)
        for leftover in list(self.work_dir.glob("req_*.json")) + \
                list(self.work_dir.glob("res_*.json")):
            leftover.unlink(missing_ok=True)
        (self.work_dir / "stop.flag").unlink(missing_ok=True)
        (self.work_dir / "ready.flag").unlink(missing_ok=True)

        self._proc = subprocess.Popen(
            [str(WORKER_PY), str(WORKER), "--serve", "--dir", str(self.work_dir)],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))

    def wait_until_ready(self, timeout: float = STARTUP_TIMEOUT) -> bool:
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            if self.ready:
                return True
            if self._proc is not None and self._proc.poll() is not None:
                return False
            time.sleep(0.25)
        return False

    def stop(self) -> None:
        if self._proc is None:
            return
        try:
            (self.work_dir / "stop.flag").write_text("stop", encoding="utf-8")
            self._proc.wait(timeout=10)
        except Exception:  # noqa: BLE001
            self._proc.kill()
        finally:
            self._proc = None

    def recognise(self, image: Path, timeout: float = REQUEST_TIMEOUT) -> Recognition:
        if not self.running:
            self.start()
        if not self.wait_until_ready():
            return Recognition(error="Recognition service failed to start.")

        request_id = uuid.uuid4().hex[:12]
        req = self.work_dir / f"req_{request_id}.json"
        res = self.work_dir / f"res_{request_id}.json"
        req.write_text(json.dumps({"image": str(image)}), encoding="utf-8")

        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            if res.exists():
                try:
                    payload = json.loads(res.read_text(encoding="utf-8"))
                except json.JSONDecodeError:
                    time.sleep(0.1)     # still being written
                    continue
                res.unlink(missing_ok=True)
                return Recognition(smiles=payload.get("smiles"),
                                   error=payload.get("error"))
            if self._proc is not None and self._proc.poll() is not None:
                return Recognition(error="Recognition service stopped unexpectedly.")
            time.sleep(0.15)

        req.unlink(missing_ok=True)
        return Recognition(error=f"Recognition timed out after {timeout:g}s.")
