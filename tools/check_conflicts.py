"""Provenance of the rows where we disagree with MA.

If the structure came from a verified CAS lookup, a disagreement points at the manual
answer. If we built the structure from the parent drug, our molecule may simply not be
the compound the vendor is selling - and the disagreement proves nothing.
"""

import json
from pathlib import Path

ROOT = Path(__file__).parent.parent
rows = {r["impurity_name"]: r
        for r in json.loads((ROOT / "data" / "resolved_structures.json")
                            .read_text(encoding="utf-8"))}

NAMES = [
    "Bosutinib Nitroso Impurity",
    "Eletriptan Nitroso Impurity",
    "Frovatriptan Nitroso Impurity 1",
    "Levomefolate Nitroso Impurity 1",
]

for name in NAMES:
    r = next((v for k, v in rows.items() if k.startswith(name)), None)
    if r is None:
        print(f"  {name}: not found")
        continue
    print(f"  {r['impurity_name'][:60]}")
    print(f"     confidence : {r.get('confidence')}   method: {r.get('method')}")
    print(f"     cas        : {r.get('cas')}   parent: {r.get('parent_drug')}")
    print(f"     candidates : {r.get('n_candidates', '-')}  "
          f"routes: {r.get('routes', '-')}")
    print(f"     category   : {r.get('category')}  ({r.get('ai_ng_day')} ng/day)")
    if r.get("published_ai_ng_day"):
        print(f"     FDA        : cat {r.get('published_category')} "
              f"({r['published_ai_ng_day']:g} ng/day)  [{r.get('published_basis')}]")
    print()

# Any row where a partially-filled manual answer coexists with a verified CAS structure?
print("  rows MA touched that have a HIGH-confidence (CAS-verified) structure:")
n = 0
for r in rows.values():
    if r.get("confidence") == "high" and r.get("scored"):
        n += 1
print(f"    {n} high-confidence scored rows in total")
