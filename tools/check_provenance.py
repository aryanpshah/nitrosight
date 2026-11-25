"""Exactly where did each scored structure come from?"""

import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).parent.parent
rows = json.loads((ROOT / "data" / "resolved_structures.json").read_text(encoding="utf-8"))

seen: dict[str, dict] = {}
for r in rows:
    seen.setdefault((r.get("catalogue_no") or "").strip()
                    or r["impurity_name"].lower(), r)
catalogue = list(seen.values())
scored = [r for r in catalogue if r.get("scored")]

print(f"  distinct compounds : {len(catalogue)}")
print(f"  scored             : {len(scored)}")
print(f"  not scored         : {len(catalogue) - len(scored)}")
print()

print("  how the SCORED structures were obtained:")
for method, n in Counter(r.get("method") for r in scored).most_common():
    label = {
        "cas_lookup": "CAS number -> PubChem (a real lookup, verified to contain N-N=O)",
        "parent_nitrosation": "parent drug name -> PubChem -> we ADDED the nitroso ourselves",
        "vendor_supplied": "supplied by the vendor",
    }.get(method, method)
    print(f"    {n:>4}  {label}")
print()

built = [r for r in scored if r.get("method") == "parent_nitrosation"]
one_site = sum(1 for r in built if (r.get("n_candidates") or 1) == 1)
print(f"  of the {len(built)} built from a parent drug:")
print(f"    {one_site:>4}  had only one possible nitrogen (unambiguous)")
print(f"    {len(built)-one_site:>4}  had several, so the molecule is a guess")
print()

print("  photos used in any of this: 0")
print(f"  photos extracted and sitting unused: "
      f"{len(list((ROOT / 'data' / 'structure_images').glob('*.png')))}")
