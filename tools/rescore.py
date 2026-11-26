"""
Re-score the resolved catalogue without going back to the network.

fetch_structures.py both RESOLVES structures (PubChem lookups, NDSRI enumeration) and
SCORES them. That is the right thing when the catalogue or the resolution logic
changes, but it is the wrong thing when only the scorer has changed: PubChem answers
vary between runs, so a re-run moves compounds for reasons that have nothing to do with
the change being tested. One run in November 2026 silently lost 14 structures that way.

This tool re-scores the structures already in data/resolved_structures.json and leaves
resolution alone, so the difference is attributable to the scorer and nothing else.

    python tools/rescore.py --dry-run     # show what would change
    python tools/rescore.py               # apply and write

Candidates: every row with several nitrosation candidates stores ALL of them,
with each one's formation route (verified on the committed data: all 69 such
rows, counts matching n_candidates). That is what lets this tool record the
name-implied candidate in the named_* fields without going back to the network.
It still does NOT re-pick the screening structure: the stored worst-case winner
is re-scored as it stands, so if a scorer change would make a DIFFERENT
candidate the worst case, re-run fetch_structures.py.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

sys.path.insert(0, str(ROOT / "tools"))

from fetch_structures import (  # noqa: E402
    NAMED_FIELDS, apply_published_limits, named_reading, score_best,
)

RESULTS = ROOT / "data" / "resolved_structures.json"

# Written by score_best; cleared before re-scoring so a row that no longer scores does
# not keep a stale category.
SCORE_FIELDS = ("scored", "scored_smiles", "category", "ai_ng_day", "potency_score",
                "factors", "triage_exit", "warnings")
LIMIT_FIELDS = ("final_ai_ng_day", "final_basis", "published_name",
                "published_ai_ng_day", "published_category", "published_basis",
                "override_note")


def ascii_(text: str) -> str:
    return text.encode("ascii", "replace").decode()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true",
                    help="report the changes without writing the file")
    args = ap.parse_args()

    rows = json.loads(RESULTS.read_text(encoding="utf-8"))
    before = {id(r): (r.get("scored"), r.get("category")) for r in rows}
    snapshot = [dict(r) for r in rows]

    changed, lost, gained = [], 0, 0
    for row, old in zip(rows, snapshot):
        smiles = row.get("smiles")
        if not smiles:
            continue

        for field in SCORE_FIELDS + LIMIT_FIELDS:
            row.pop(field, None)
        row.update(score_best([smiles]))

        was, now = old.get("scored"), row.get("scored")
        if was and not now:
            lost += 1
        elif now and not was:
            gained += 1
        elif now and was and old.get("category") != row.get("category"):
            changed.append((row["impurity_name"], old["category"], row["category"],
                            old.get("ai_ng_day"), row.get("ai_ng_day")))

    apply_published_limits(rows)

    # The identity reading, from the STORED candidates. Done after the limits so
    # the new keys land at the end of each record and no existing key moves.
    # Cleared first, like SCORE_FIELDS, so a re-run never leaves a stale value.
    n_named = n_named_differs = 0
    for row in rows:
        for field in NAMED_FIELDS:
            row.pop(field, None)
        if row.get("candidates"):
            row.update(named_reading(row["impurity_name"], row["candidates"],
                                     row.get("smiles")))
            if "named_smiles" in row:
                n_named += 1
                n_named_differs += bool(row.get("named_differs_from_reported"))

    print(f"\n{'=' * 72}\n  RESCORE\n{'=' * 72}")
    print(f"    {len(rows)} rows, {sum(1 for r in rows if r.get('smiles'))} with a structure")
    print(f"    category changed : {len(changed)}")
    print(f"    newly scored     : {gained}")
    print(f"    no longer scored : {lost}")
    print(f"    name-implied reading recorded : {n_named}"
          f"  ({n_named_differs} point at a different structure from the screening one)")

    if changed:
        print(f"\n{'-' * 72}")
        for name, a, b, ai_a, ai_b in sorted(changed, key=lambda c: (c[1], c[2])):
            arrow = f"Cat {a} -> {b}"
            limits = (f"   {ai_a:g} -> {ai_b:g} ng/day"
                      if ai_a is not None and ai_b is not None else "")
            print(f"    {ascii_(name)[:50]:<52} {arrow}{limits}")

    dist = Counter(r["category"] for r in rows if r.get("scored"))
    print(f"\n{'-' * 72}\n  CATEGORY DISTRIBUTION")
    for cat in sorted(dist):
        print(f"    Category {cat}:  {dist[cat]}")
    tight = sum(n for c, n in dist.items() if c in (1, 2))
    print(f"\n    {tight} impurities at 100 ng/day or tighter")

    if args.dry_run:
        print("\n  --dry-run: nothing written")
        return 0

    RESULTS.write_text(json.dumps(rows, indent=1), encoding="utf-8")
    print(f"\n  written to {RESULTS.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
