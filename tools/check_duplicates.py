"""Are the repeated names in the results real duplicate catalogue rows, or a bug?"""

import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).parent.parent
rows = json.loads((ROOT / "data" / "resolved_structures.json").read_text(encoding="utf-8"))

print(f"  rows in results file        : {len(rows)}")

by_name = Counter(r["impurity_name"] for r in rows)
dup_names = {n: c for n, c in by_name.items() if c > 1}
print(f"  distinct impurity names     : {len(by_name)}")
print(f"  names appearing >1 time     : {len(dup_names)}")
print(f"  extra rows from repeats     : {sum(dup_names.values()) - len(dup_names)}")
print()

by_cat = Counter(r.get("catalogue_no", "") for r in rows if r.get("catalogue_no"))
dup_cat = {n: c for n, c in by_cat.items() if c > 1}
print(f"  distinct catalogue numbers  : {len(by_cat)}")
print(f"  catalogue numbers repeated  : {len(dup_cat)}")
print()

print("  most repeated names:")
for name, n in Counter(dup_names).most_common(6):
    print(f"    x{n}  {name[:60]}")
print()

# Are repeated rows byte-identical in the source columns?
sample = next(iter(dup_names), None)
if sample:
    same = [r for r in rows if r["impurity_name"] == sample]
    print(f"  '{sample[:45]}' appears {len(same)} times:")
    for r in same:
        print(f"    parent={r['parent_drug'][:18]:<18} cas={r['cas'][:14]:<14} "
              f"cat={r.get('catalogue_no','')}")
