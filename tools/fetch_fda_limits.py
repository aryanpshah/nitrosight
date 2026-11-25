"""
Scrape FDA's published nitrosamine AI limits.

FDA's September 2024 revision of "Control of Nitrosamine Impurities in Human Drugs"
removed the static limits table and now points to a web page that is updated
periodically. Those published, compound-specific limits take precedence over a CPCA
prediction, so the tool needs them.

    python tools/fetch_fda_limits.py

Writes data/fda_published_limits.json with a fetch date, because this is a snapshot of
a living page and will go stale.
"""

from __future__ import annotations

import json
import re
import urllib.request
from datetime import date
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).parent.parent
OUT = ROOT / "data" / "fda_published_limits.json"
URL = ("https://www.fda.gov/regulatory-information/search-fda-guidance-documents/"
       "cder-nitrosamine-impurity-acceptable-intake-limits")
UA = {"User-Agent": "impurity-scanner/0.1 (NDSRI screening research tool)"}


class TableParser(HTMLParser):
    """Collect every <table> as a list of rows, each row a list of cell strings."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.tables: list[list[list[str]]] = []
        self._table: list[list[str]] | None = None
        self._row: list[str] | None = None
        self._cell: list[str] | None = None

    def handle_starttag(self, tag, attrs):
        if tag == "table":
            self._table = []
        elif tag == "tr" and self._table is not None:
            self._row = []
        elif tag in ("td", "th") and self._row is not None:
            self._cell = []

    def handle_endtag(self, tag):
        if tag == "table" and self._table is not None:
            self.tables.append(self._table)
            self._table = None
        elif tag == "tr" and self._row is not None:
            if any(c.strip() for c in self._row):
                self._table.append(self._row)
            self._row = None
        elif tag in ("td", "th") and self._cell is not None:
            text = re.sub(r"\s+", " ", "".join(self._cell)).strip()
            self._row.append(text)
            self._cell = None

    def handle_data(self, data):
        if self._cell is not None:
            self._cell.append(data)


NUM = re.compile(r"^\s*([\d,]+(?:\.\d+)?)\s*$")


def as_limit(text: str) -> float | None:
    m = NUM.match(text.replace("ng/day", "").strip())
    return float(m.group(1).replace(",", "")) if m else None


def main() -> int:
    req = urllib.request.Request(URL, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        html = r.read().decode("utf-8", "replace")
    print(f"  fetched {len(html):,} bytes")

    parser = TableParser()
    parser.feed(html)
    print(f"  {len(parser.tables)} tables found")

    entries: list[dict] = []
    for t_i, table in enumerate(parser.tables):
        if len(table) < 2:
            continue
        header = [c.lower() for c in table[0]]
        print(f"    table {t_i}: {len(table)-1} rows | header: {header[:5]}")

        # Find the column holding the AI limit, and a plausible name column.
        ai_col = next((i for i, h in enumerate(header)
                       if "ai" in h and ("limit" in h or "ng" in h)), None)
        if ai_col is None:
            ai_col = next((i for i, h in enumerate(header) if "ng/day" in h), None)
        if ai_col is None:
            continue

        for row in table[1:]:
            if len(row) <= ai_col:
                continue
            limit = as_limit(row[ai_col])
            if limit is None:
                continue
            cells = [c for i, c in enumerate(row) if i != ai_col and c]
            if not cells:
                continue
            entries.append({
                "table": t_i,
                "columns": dict(zip(header, row)),
                "name": max(cells, key=len),
                "ai_ng_day": limit,
            })

    # FDA's revision table records when each compound was added to the page. That is
    # the only dated signal available, and it is what makes a genuine temporal
    # hold-out possible: a compound FDA published after our rules were frozen cannot
    # have influenced them. See tools/holdout.py.
    revisions = []
    for table in parser.tables:
        if len(table) < 2:
            continue
        header = [c.lower() for c in table[0]]
        name_col = next((i for i, h in enumerate(header) if "name" in h), None)
        date_col = next((i for i, h in enumerate(header)
                         if "date" in h and "post" in h), None)
        change_col = next((i for i, h in enumerate(header) if "change" in h), None)
        if name_col is None or date_col is None:
            continue
        for row in table[1:]:
            if len(row) <= max(name_col, date_col):
                continue
            name, posted = row[name_col].strip(), row[date_col].strip()
            if not name or not posted:
                continue
            revisions.append({
                "name": name,
                "date_posted": posted,
                "change": row[change_col].strip() if (
                    change_col is not None and len(row) > change_col) else "",
            })
    print(f"  {len(revisions)} dated revision entries captured")

    payload = {
        "_source_url": URL,
        "_fetched": date.today().isoformat(),
        "_warning": "FDA updates this page periodically. Re-run this script to refresh; "
                    "treat the values as a dated snapshot, not a permanent record.",
        "_count": len(entries),
        "entries": entries,
        "revisions": revisions,
    }
    text = json.dumps(payload, indent=1)
    OUT.write_text(text, encoding="utf-8")
    print(f"\n  {len(entries)} limit entries -> {OUT.relative_to(ROOT)}")

    # resolve/fda_limits.py falls back to this copy when the package is installed
    # rather than run from a clone, so the two must not drift apart.
    # tools/validate_data.py fails if they do.
    packaged = ROOT / "resolve" / "data" / "fda_published_limits.json"
    packaged.parent.mkdir(parents=True, exist_ok=True)
    packaged.write_text(text, encoding="utf-8")
    print(f"  mirrored for packaging  -> {packaged.relative_to(ROOT)}")
    return 0 if entries else 1


if __name__ == "__main__":
    raise SystemExit(main())
