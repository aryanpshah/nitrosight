"""
Validate the CPCA scorer against FDA's own published potency categories.

FDA's limits page lists a recommended potency category for ~280 named nitrosamines.
Wherever one of those matches a catalogue compound we scored, we have a direct
ground-truth comparison: does our implementation of the CPCA reach the same category
FDA did?

This is a far stronger check than the two worked examples in the guidance.

    python tools/validate_vs_fda.py
    python tools/validate_vs_fda.py --show-agreements
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from resolve import FDALimits  # noqa: E402

RESULTS = ROOT / "data" / "resolved_structures.json"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--show-agreements", action="store_true")
    args = ap.parse_args()

    fda = FDALimits()
    rows = json.loads(RESULTS.read_text(encoding="utf-8"))

    seen: dict[str, dict] = {}
    for r in rows:
        if r.get("scored"):
            seen.setdefault(
                (r.get("catalogue_no") or "").strip() or r["impurity_name"].lower(), r)
    scored = list(seen.values())

    print(f"  FDA snapshot : {fda.fetched}  ({len(fda)} named limits)")
    if warning := fda.staleness_warning():
        print(f"  ! {warning}")
    print(f"  scored rows  : {len(scored)}\n")

    matched, agree, disagree = [], [], []
    for r in scored:
        pub = fda.lookup(r["impurity_name"])
        if pub is None or pub.category is None:
            continue
        matched.append((r, pub))
        (agree if r["category"] == pub.category else disagree).append((r, pub))

    if not matched:
        print("  No catalogue compound matched an FDA entry that carries a category.")
        return 1

    rate = len(agree) / len(matched)
    print("=" * 74)
    print(f"  MATCHED {len(matched)} compounds against FDA's published categories")
    print(f"  AGREEMENT: {len(agree)}/{len(matched)}  ({rate:.0%})")
    print("=" * 74)
    print()

    if disagree:
        print(f"--- {len(disagree)} DISAGREEMENTS ---\n")
        direction = Counter(
            "stricter than FDA" if r["category"] < p.category else "looser than FDA"
            for r, p in disagree)
        for label, n in direction.most_common():
            print(f"    {n:>3}  our prediction is {label}")
        print()
        for r, p in sorted(disagree, key=lambda x: x[0]["category"]):
            arrow = "STRICTER" if r["category"] < p.category else "LOOSER  "
            print(f"  [{arrow}] {r['impurity_name'][:48]}")
            print(f"      ours: cat {r['category']} ({r['ai_ng_day']:g} ng/day, "
                  f"score {r.get('potency_score')})   "
                  f"FDA: cat {p.category} ({p.ai_ng_day:g} ng/day)")
            factors = r.get("factors") or []
            if factors:
                print("      factors: " +
                      "; ".join(f"{f['score']:+d} {f['name'][:38]}" for f in factors))
            elif r.get("triage_exit"):
                print(f"      triage: {r['triage_exit']}")
            print(f"      confidence: {r.get('confidence')} ({r.get('method')})")
            print()

    if args.show_agreements:
        print(f"--- {len(agree)} AGREEMENTS ---")
        for r, p in sorted(agree, key=lambda x: x[0]["category"]):
            print(f"  cat {r['category']}  {r['impurity_name'][:52]:<52} "
                  f"{r['ai_ng_day']:>6g} ng/day  [{r.get('confidence')}]")

    print("=" * 74)
    print("  AGREEMENT BY OUR OWN CONFIDENCE LABEL")
    print("=" * 74)
    for label in ("high", "medium", "low"):
        sub = [(r, p) for r, p in matched if r.get("confidence") == label]
        if not sub:
            continue
        ok = sum(1 for r, p in sub if r["category"] == p.category)
        print(f"    {label:<8} {ok:>3}/{len(sub):<3}  ({ok/len(sub):.0%})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
