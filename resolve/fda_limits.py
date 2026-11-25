"""
FDA's published acceptable intake limits.

Where FDA has already published a limit for a specific nitrosamine, that value governs
and a CPCA prediction must not override it. CPCA exists for compounds that *lack* an
established limit; the guidance is explicit that compound-specific data takes precedence.

Data comes from `tools/fetch_fda_limits.py`, which scrapes FDA's limits page. That page
is updated periodically, so the snapshot carries a fetch date and goes stale.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from datetime import date
from pathlib import Path

def _data_path() -> Path:
    """Where the published-limits snapshot lives.

    Running from a clone it sits in data/ at the repository root. In an
    installed copy that directory is outside the package, so a copy is shipped
    as package data inside resolve/ and used as the fallback.
    """
    repo = Path(__file__).parent.parent / "data" / "fda_published_limits.json"
    if repo.exists():
        return repo
    return Path(__file__).parent / "data" / "fda_published_limits.json"


DATA = _data_path()

# Which scraped table a row came from, and what that means.
TABLE_BASIS = {
    0: "published",   # FDA recommended AI limit, usually with a potency category
    1: "surrogate",   # derived by read-across from a surrogate compound
    2: "interim",     # time-limited interim limit
}

STALE_AFTER_DAYS = 90


@dataclass
class PublishedLimit:
    name: str
    ai_ng_day: float
    category: int | None
    basis: str
    source: str = ""
    surrogate: str = ""
    duration: str = ""

    @property
    def is_interim(self) -> bool:
        return self.basis == "interim"


def normalise(name: str) -> str:
    """Collapse naming variants so 'N-Nitroso Nortriptyline' matches
    'N-nitroso-nortriptyline'."""
    low = name.lower()
    low = re.sub(r"\(.*?\)", " ", low)          # drop parenthetical abbreviations
    low = low.replace("nitrosamine", "nitroso")
    low = re.sub(r"[^a-z0-9]+", "", low)
    # Common spelling variants seen across FDA and vendor naming.
    low = low.replace("nitrosoamino", "nitrosamino")
    low = low.replace("desmethyl", "desmethyl").replace("ndesmethyl", "desmethyl")
    return low


class FDALimits:
    def __init__(self, path: Path = DATA):
        payload = json.loads(path.read_text(encoding="utf-8"))
        self.fetched = payload.get("_fetched", "unknown")
        self.source_url = payload.get("_source_url", "")
        self._by_name: dict[str, PublishedLimit] = {}

        for e in payload["entries"]:
            cols = e.get("columns", {})
            cat_raw = (cols.get("potency category") or "").strip()
            entry = PublishedLimit(
                name=e["name"],
                ai_ng_day=e["ai_ng_day"],
                category=int(cat_raw) if cat_raw.isdigit() else None,
                basis=TABLE_BASIS.get(e["table"], "published"),
                source=cols.get("source*") or cols.get("source", ""),
                surrogate=cols.get("surrogate", ""),
                duration=cols.get("estimated duration***", ""),
            )
            key = normalise(entry.name)
            if not key:
                continue
            # Prefer a definitive published limit over an interim one.
            existing = self._by_name.get(key)
            if existing is None or (existing.is_interim and not entry.is_interim):
                self._by_name[key] = entry

    def __len__(self) -> int:
        return len(self._by_name)

    @property
    def days_old(self) -> int:
        try:
            return (date.today() - date.fromisoformat(self.fetched)).days
        except ValueError:
            return 10_000

    @property
    def is_stale(self) -> bool:
        return self.days_old > STALE_AFTER_DAYS

    def lookup(self, name: str) -> PublishedLimit | None:
        return self._by_name.get(normalise(name))

    def staleness_warning(self) -> str | None:
        if not self.is_stale:
            return None
        return (f"FDA published-limits snapshot is {self.days_old} days old "
                f"(fetched {self.fetched}). Re-run tools/fetch_fda_limits.py.")
