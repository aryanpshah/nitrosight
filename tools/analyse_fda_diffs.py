"""
Drill into where we and FDA's reference implementation disagree.

A disagreement is one of three things, and they need separating before any of it means
anything:

  * we picked a different N-nitroso site to score (then every feature differs, and the
    real difference is one of site selection, not of rules)
  * we implement a rule differently
  * one of us is simply wrong

Alpha-hydrogen counts are the tell for the first case: if those differ, we are almost
certainly looking at different nitrogens.

    python tools/analyse_fda_diffs.py            # overview
    python tools/analyse_fda_diffs.py <feature>  # examples for one feature
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).parent.parent
REPORT = ROOT / "data" / "fda_benchmark_comparison.json"


def main() -> int:
    records = json.loads(REPORT.read_text(encoding="utf-8"))

    same_site = [r for r in records
                 if r.get("fda_alpha_h") and r.get("our_alpha_h")
                 and r["fda_alpha_h"] == r["our_alpha_h"]]
    diff_site = [r for r in records
                 if r.get("fda_alpha_h") and r.get("our_alpha_h")
                 and r["fda_alpha_h"] != r["our_alpha_h"]]
    declined = [r for r in records if r["verdict"] == "we declined"]

    print("=" * 76)
    print("  WHAT KIND OF DISAGREEMENT IS IT?")
    print("=" * 76)
    print(f"    {len(records):>4}  structures differ in some way")
    print(f"    {len(diff_site):>4}  ... alpha-hydrogen counts differ "
          f"(we scored a DIFFERENT nitroso site)")
    print(f"    {len(same_site):>4}  ... same site, so a genuine rule difference")
    print(f"    {len(declined):>4}  ... we declined to score, FDA did not")
    print()

    cat_diff_same_site = [r for r in same_site if r["verdict"] == "category differs"]
    print(f"    of the same-site cases, {len(cat_diff_same_site)} also differ on "
          f"the final category")
    print()

    if len(sys.argv) > 1:
        needle = sys.argv[1].lower()
        hits = [r for r in same_site
                if any(needle in d["feature"].lower() for d in r["feature_diffs"])]
        print("=" * 76)
        print(f"  SAME-SITE EXAMPLES: {needle}")
        print("=" * 76)
        for r in hits[:12]:
            diff = next(d for d in r["feature_diffs"]
                        if needle in d["feature"].lower())
            print(f"\n  {r['name'][:62]}")
            print(f"     {r['smiles'][:96]}")
            print(f"     alpha-H {r['our_alpha_h']} (same)   "
                  f"FDA={diff['fda']}  ours={diff['ours']}")
            print(f"     category  FDA={r['fda_category']}  ours={r['our_category']}"
                  f"    score FDA={r.get('fda_score')} ours={r.get('our_score')}")
        print(f"\n  {len(hits)} same-site cases for this feature")
        return 0

    # Which rules disagree, counting only same-site comparisons.
    print("=" * 76)
    print("  RULE DIFFERENCES (same nitroso site only)")
    print("=" * 76)
    counts = Counter()
    direction = Counter()
    for r in same_site:
        for d in r["feature_diffs"]:
            counts[d["feature"]] += 1
            direction[(d["feature"], "FDA scores it, we do not" if d["fda"] > d["ours"]
                       else "we score it, FDA does not")] += 1
    for feature, n in counts.most_common():
        print(f"\n    {n:>4}  {feature[:68]}")
        for (f, way), m in direction.most_common():
            if f == feature:
                print(f"          {m:>4}  {way}")

    print()
    print("  Rerun with a feature keyword for examples, e.g.:")
    print("    python tools/analyse_fda_diffs.py carboxylic")
    print("    python tools/analyse_fda_diffs.py chains")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
