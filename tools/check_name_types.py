"""
Some catalogue entries are vendor shorthand ("Afatinib Nitroso Impurity 2"), which says
nothing about which nitrogen carries the nitroso group. Others are full systematic IUPAC
names, which say exactly where it is.

If a meaningful share are systematic, a name-to-structure parser (OPSIN) could resolve
them without needing anything from the vendor.
"""

import json
import re
from pathlib import Path

ROOT = Path(__file__).parent.parent
rows = json.loads((ROOT / "data" / "resolved_structures.json").read_text(encoding="utf-8"))

seen: dict[str, dict] = {}
for r in rows:
    seen.setdefault((r.get("catalogue_no") or "").strip()
                    or r["impurity_name"].lower(), r)
catalogue = list(seen.values())

# Systematic names carry locants and substituent brackets: "N-(2-((2-(dimethyl..."
LOCANT = re.compile(r"\d+-\(|\)-\d|-\d+[HH]-|\d+,\d+-|\(\d+[RSEZ]\)")
IUPAC_WORDS = ("yl)", "amino", "oxy", "carbox", "phenyl", "methyl", "ethyl",
               "pyrimidin", "indol", "piperazin", "yl-", "oate", "amide")


def is_systematic(name: str) -> bool:
    if len(name) < 40:
        return False
    if not LOCANT.search(name):
        return False
    return sum(1 for w in IUPAC_WORDS if w in name.lower()) >= 2


systematic = [r for r in catalogue if is_systematic(r["impurity_name"])]
shorthand = [r for r in catalogue if not is_systematic(r["impurity_name"])]

print(f"  catalogue compounds        : {len(catalogue)}")
print(f"  systematic IUPAC-style name: {len(systematic)}  "
      f"({len(systematic)/len(catalogue):.0%})")
print(f"  vendor shorthand           : {len(shorthand)}")
print()


def bucket(rows_):
    weak = [r for r in rows_
            if not r.get("scored") or r.get("confidence") in ("low", "unresolved",
                                                              "conflict")]
    return len(weak)


print(f"  currently weak or unresolved, by name type:")
print(f"    systematic names : {bucket(systematic):>4} of {len(systematic)}")
print(f"    shorthand names  : {bucket(shorthand):>4} of {len(shorthand)}")
print()

print("  systematic names that are currently weak (OPSIN could likely fix these):")
n = 0
for r in systematic:
    if r.get("scored") and r.get("confidence") not in ("low", "unresolved", "conflict"):
        continue
    print(f"    [{r.get('confidence')}] {r['impurity_name'][:88]}")
    n += 1
    if n >= 8:
        break

print()
print("  shorthand names, for contrast (nothing to parse):")
for r in shorthand[:6]:
    if (r.get("n_candidates") or 0) > 1:
        print(f"    {r['impurity_name'][:52]:<52} "
              f"{r['n_candidates']} possible sites")
