"""Peek at the scraped FDA limits: shape, sample rows, and key compounds."""

import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).parent.parent
d = json.loads((ROOT / "data" / "fda_published_limits.json").read_text(encoding="utf-8"))
entries = d["entries"]

print(f"  fetched {d['_fetched']}  |  {len(entries)} entries")
print(f"  by table: {dict(Counter(e['table'] for e in entries))}\n")

for t in sorted({e["table"] for e in entries}):
    rows = [e for e in entries if e["table"] == t]
    print(f"--- table {t} ({len(rows)} rows) ---")
    print(f"    columns: {list(rows[0]['columns'].keys())}")
    for e in rows[:3]:
        print(f"      {e['name'][:52]:<52} {e['ai_ng_day']:>7g} ng/day")
    print()

print("--- key small-molecule nitrosamines ---")
for needle in ("dimethylamine", "diethylamine", "methylphenylamine",
               "aminobutyric", "isopropylethylamine", "diisopropylamine",
               "dibutylamine"):
    hits = [e for e in entries if needle in e["name"].lower()]
    for e in hits[:2]:
        print(f"    {e['name'][:56]:<56} {e['ai_ng_day']:>7g} ng/day")

print("\n--- compounds my tool put in Category 1 ---")
for needle in ("fluoxetine", "duloxetine", "atomoxetine", "nortriptyline",
               "desipramine", "diphenhydramine", "benzathine"):
    hits = [e for e in entries if needle in e["name"].lower()]
    for e in hits[:2]:
        cat = e["columns"].get("potency category", "")
        print(f"    {e['name'][:52]:<52} {e['ai_ng_day']:>7g} ng/day  cat={cat}")
