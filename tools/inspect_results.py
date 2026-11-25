"""Summarise data/resolved_structures.json: what resolved, what didn't, and why."""

import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).parent.parent
RESULTS = ROOT / "data" / "resolved_structures.json"


def main() -> None:
    show = sys.argv[1] if len(sys.argv) > 1 else "failures"
    rows = json.loads(RESULTS.read_text(encoding="utf-8"))

    if show == "failures":
        bad = [r for r in rows if not r.get("scored")]
        print(f"  {len(bad)} of {len(rows)} rows not scored\n")
        reasons = Counter(
            (r.get("issue") or "no issue recorded").split(";")[0][:70] for r in bad)
        for reason, n in reasons.most_common():
            print(f"  {n:>5}  {reason}")
        print("\n  examples:")
        for r in bad[:15]:
            print(f"    {r.get('confidence','?'):<11} {r['impurity_name'][:40]:<40} "
                  f"parent={r['parent_drug'][:20]}")

    elif show == "strict":
        strict = sorted((r for r in rows if r.get("scored") and r["category"] <= 2),
                        key=lambda r: (r["category"], r["impurity_name"]))
        print(f"  {len(strict)} impurities in Category 1 or 2\n")
        for r in strict:
            print(f"  Cat {r['category']}  {r['ai_ng_day']:>6} ng/day  "
                  f"{r['impurity_name'][:44]:<44} [{r.get('confidence')}]")
            factors = " ".join(f"{f['score']:+d}" for f in r.get("factors", []))
            print(f"          score {r.get('potency_score')}  ({factors})  "
                  f"{r.get('catalogue_no','')}")

    elif show == "uncertain":
        flagged = [r for r in rows
                   if any(f.get("uncertain") for f in r.get("factors", []))]
        print(f"  {len(flagged)} scored rows contain a [needs review] factor\n")
        for r in flagged[:40]:
            names = [f["name"] for f in r["factors"] if f.get("uncertain")]
            print(f"  Cat {r['category']}  {r['impurity_name'][:44]:<44} {'; '.join(names)[:40]}")

    else:
        print("usage: inspect_results.py [failures|strict|uncertain]")


if __name__ == "__main__":
    main()
