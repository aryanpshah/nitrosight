"""
Cached, rate-limited PubChem client.

PubChem asks for no more than 5 requests/second and 400/minute. We stay well under
both. Every response is cached to disk, so re-running a resolution pass costs nothing
for queries already seen and the job is resumable after an interruption.
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
USER_AGENT = "impurity-scanner/0.1 (NDSRI screening research tool)"
PROPERTIES = "SMILES,ConnectivitySMILES,MolecularFormula,IUPACName"

CACHE_PATH = Path(__file__).parent.parent / "data" / "pubchem_cache.json"

MIN_INTERVAL = 0.25   # seconds between requests -> 4/sec
MAX_RETRIES = 3


@dataclass
class Hit:
    query: str
    cid: int | None
    smiles: str | None
    formula: str | None
    iupac_name: str | None
    status: str          # "ok" | "not_found" | "error"
    message: str = ""

    @property
    def found(self) -> bool:
        return self.status == "ok" and bool(self.smiles)


class PubChemClient:
    def __init__(self, cache_path: Path = CACHE_PATH, verbose: bool = True):
        self.cache_path = cache_path
        self.verbose = verbose
        self._cache: dict[str, dict] = {}
        self._last_request = 0.0
        self._live_requests = 0
        self._load()

    # cache

    def _load(self) -> None:
        if self.cache_path.exists():
            self._cache = json.loads(self.cache_path.read_text(encoding="utf-8"))
            if self.verbose:
                print(f"  cache: {len(self._cache)} entries loaded")

    def save(self) -> None:
        self.cache_path.parent.mkdir(parents=True, exist_ok=True)
        self.cache_path.write_text(
            json.dumps(self._cache, indent=1, sort_keys=True), encoding="utf-8")

    @property
    def live_requests(self) -> int:
        return self._live_requests

    # fetching

    def _throttle(self) -> None:
        elapsed = time.monotonic() - self._last_request
        if elapsed < MIN_INTERVAL:
            time.sleep(MIN_INTERVAL - elapsed)
        self._last_request = time.monotonic()

    def _fetch(self, query: str) -> dict:
        url = f"{BASE}/name/{urllib.parse.quote(query, safe='')}/property/{PROPERTIES}/JSON"
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})

        for attempt in range(1, MAX_RETRIES + 1):
            self._throttle()
            self._live_requests += 1
            try:
                with urllib.request.urlopen(req, timeout=30) as r:
                    props = json.loads(r.read())["PropertyTable"]["Properties"][0]
                return {
                    "status": "ok",
                    "cid": props.get("CID"),
                    "smiles": props.get("SMILES") or props.get("ConnectivitySMILES"),
                    "formula": props.get("MolecularFormula"),
                    "iupac_name": props.get("IUPACName"),
                }
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    return {"status": "not_found", "message": "no CID matches this name"}
                if e.code in (429, 503) and attempt < MAX_RETRIES:
                    time.sleep(2.0 * attempt)   # server is throttling; back off
                    continue
                return {"status": "error", "message": f"HTTP {e.code}"}
            except Exception as e:  # noqa: BLE001
                if attempt < MAX_RETRIES:
                    time.sleep(1.5 * attempt)
                    continue
                return {"status": "error", "message": f"{type(e).__name__}: {e}"}

        return {"status": "error", "message": "retries exhausted"}

    def lookup(self, query: str) -> Hit:
        """Resolve a name or CAS number to a structure. Cached across runs."""
        key = query.strip().lower()
        if not key:
            return Hit(query, None, None, None, None, "not_found", "empty query")

        if key not in self._cache:
            self._cache[key] = self._fetch(query.strip())

        d = self._cache[key]
        return Hit(
            query=query,
            cid=d.get("cid"),
            smiles=d.get("smiles"),
            formula=d.get("formula"),
            iupac_name=d.get("iupac_name"),
            status=d.get("status", "error"),
            message=d.get("message", ""),
        )
