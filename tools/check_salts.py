"""Did salt counter-ions leak into the scoring?

'Carboxylic acid group anywhere on molecule' is +3. If PubChem returned a maleate or
citrate salt, that acid sits on the counter-ion, not on the nitrosamine - and inflates
the score.
"""

import json
from pathlib import Path

ROOT = Path(__file__).parent.parent
rows = json.loads((ROOT / "data" / "resolved_structures.json").read_text(encoding="utf-8"))

multi = [r for r in rows if r.get("smiles") and "." in r["smiles"]]
acid_and_multi = [
    r for r in multi
    if any("Carboxylic acid" in f["name"] for f in (r.get("factors") or []))
]

print(f"  scored rows with a multi-fragment SMILES : {len(multi)}")
print(f"  ... of which scored a carboxylic acid    : {len(acid_and_multi)}")
print()

for r in multi[:12]:
    frags = r["smiles"].split(".")
    acid = any("Carboxylic acid" in f["name"] for f in (r.get("factors") or []))
    print(f"  {r['impurity_name'][:42]:<42} cat={r.get('category')} "
          f"frags={len(frags)} acid={'YES' if acid else 'no'}")
    print(f"      {r['smiles'][:100]}")

print()
print("  specific suspects:")
for needle in ("Carbinoxamine", "Dabigatran"):
    for r in rows:
        if needle.lower() in r["impurity_name"].lower() and r.get("smiles"):
            print(f"    {r['impurity_name'][:50]}")
            print(f"      {r['smiles'][:120]}")
            for f in (r.get("factors") or []):
                print(f"        {f['score']:+d} {f['name']}")
            break
