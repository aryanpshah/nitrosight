"""What counts as doubtful, and how many are there of each kind?

Dr Ali asked for the uncertain structures tabled separately so a chemist can settle
them. Only cases with a specific, answerable question belong there.
"""

import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).parent.parent
rows = json.loads((ROOT / "data" / "resolved_structures.json").read_text(encoding="utf-8"))
index = json.loads((ROOT / "data" / "structure_images_index.json").read_text(encoding="utf-8"))
has_image = {e["impurity_name"].strip().lower() for e in index}

seen: dict[str, dict] = {}
for r in rows:
    seen.setdefault((r.get("catalogue_no") or "").strip()
                    or r["impurity_name"].lower(), r)
catalogue = list(seen.values())


def kind(r: dict) -> str | None:
    if r.get("image_disagrees"):
        return "image disagrees with our verified structure"
    if r.get("confidence") == "conflict":
        return "CAS resolves to a structure with no N-nitroso group"
    if r.get("scored") and r.get("confidence") == "low":
        return "several possible nitrosation sites, we guessed"
    if r.get("method") == "image_recognition":
        return "structure read from a picture, unverified"
    if not r.get("scored"):
        return "no structure at all"
    return None


counts = Counter(kind(r) for r in catalogue)
print(f"  distinct compounds: {len(catalogue)}\n")
for k, n in counts.most_common():
    label = k or "confident, nothing to ask"
    print(f"    {n:>4}  {label}")

answerable = [r for r in catalogue
              if kind(r) and kind(r) != "no structure at all"]
with_pic = sum(1 for r in answerable
               if r["impurity_name"].strip().lower() in has_image)
print()
print(f"  with a specific question a chemist could answer : {len(answerable)}")
print(f"     ... of which we have the structure picture   : {with_pic}")
