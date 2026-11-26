"""
Export the scored catalogue as a deduplicated CSV.

The source workbook contains ~100 exact duplicate rows, so results are collapsed on
catalogue number (falling back to impurity name) before anything is counted or written.

    python tools/export_report.py              # everything, to data/scored_catalogue.csv
    python tools/export_report.py --strict     # Category 1-2 only
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).parent.parent
RESULTS = ROOT / "data" / "resolved_structures.json"

COLUMNS = [
    "catalogue_no", "impurity_name", "parent_drug", "cas",
    "final_ai_ng_day", "final_basis",
    "fda_published_category", "fda_published_ai_ng_day",
    "cpca_category", "cpca_ai_ng_day", "potency_score", "factors",
    "confidence", "method", "needs_review", "smiles", "issue", "override_note",
    # Appended AFTER every existing column, so no existing column changes
    # position and a reader that goes by column index keeps working.
    #
    # selected_by: how the reported structure was chosen among the candidates
    # (only_candidate / worst_case). The rest is the identity reading, recorded
    # alongside the screening category above. named_category is NOT a limit to
    # apply; it says which molecule the vendor name points at, for matching
    # against regulators' tables.
    "selected_by", "named_route", "named_smiles", "named_category", "named_ai_ng_day",
]


def effective_ai(row: dict) -> float | None:
    """The limit that actually governs: FDA's published value if there is one,
    otherwise the CPCA prediction."""
    return row.get("final_ai_ng_day") or (row.get("ai_ng_day") if row.get("scored") else None)


def key(row: dict) -> str:
    return (row.get("catalogue_no") or "").strip() or row["impurity_name"].strip().lower()


def dedupe(rows: list[dict]) -> list[dict]:
    seen: dict[str, dict] = {}
    for row in rows:
        seen.setdefault(key(row), row)
    return list(seen.values())


def flatten(row: dict) -> dict:
    factors = row.get("factors") or []
    return {
        "catalogue_no": row.get("catalogue_no", ""),
        "impurity_name": row.get("impurity_name", ""),
        "parent_drug": row.get("parent_drug", ""),
        "cas": row.get("cas", ""),
        "final_ai_ng_day": effective_ai(row) or "",
        "final_basis": row.get("final_basis", ""),
        "fda_published_category": row.get("published_category") or "",
        "fda_published_ai_ng_day": row.get("published_ai_ng_day") or "",
        "cpca_category": row.get("category", ""),
        "cpca_ai_ng_day": row.get("ai_ng_day", ""),
        "override_note": row.get("override_note", ""),
        "potency_score": row.get("potency_score", ""),
        "factors": "; ".join(f"{f['score']:+d} {f['name']}" for f in factors)
                   or (row.get("triage_exit") or ""),
        "confidence": row.get("confidence", ""),
        "method": row.get("method", ""),
        # How the reported structure was chosen among the candidates:
        # only_candidate / vendor_name_route / worst_case. Without this a
        # consumer cannot tell a determined structure from a worst-case pick.
        "selected_by": row.get("selected_by", ""),
        "named_route": row.get("named_route", "") or "",
        "named_smiles": row.get("named_smiles", "") or "",
        "named_category": row.get("named_category", "") or "",
        "named_ai_ng_day": row.get("named_ai_ng_day", "") or "",
        "needs_review": "yes" if any(f.get("uncertain") for f in factors) else "",
        "smiles": row.get("smiles", "") or "",
        "issue": row.get("issue", ""),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--strict", action="store_true",
                    help="Category 1-2 only")
    args = ap.parse_args()

    rows = dedupe(json.loads(RESULTS.read_text(encoding="utf-8")))
    scored = [r for r in rows if r.get("scored")]

    # "Strict" means the limit that actually governs is 100 ng/day or tighter, after
    # any FDA published value has overridden the prediction.
    strict_rows = [r for r in rows
                   if (ai := effective_ai(r)) is not None and ai <= 100]
    selected = strict_rows if args.strict else rows
    selected.sort(key=lambda r: (effective_ai(r) or 1e9, r.get("impurity_name", "")))

    out = ROOT / "data" / ("category_1_2_hits.csv" if args.strict
                           else "scored_catalogue.csv")
    with out.open("w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        for row in selected:
            writer.writerow(flatten(row))

    # summary
    conf = Counter(r.get("confidence", "?") for r in rows)
    with_limit = [r for r in rows if effective_ai(r) is not None]
    basis = Counter(r.get("final_basis", "?") for r in with_limit)
    limits = Counter(effective_ai(r) for r in with_limit)

    print(f"  distinct impurities         : {len(rows)}")
    print(f"  with a governing AI limit   : {len(with_limit)}  "
          f"({len(with_limit)/len(rows):.0%})")
    print()
    print("  basis for the governing limit:")
    for k, v in basis.most_common():
        print(f"    {k:<22} {v:>4}")
    print()
    print("  governing limits:")
    for ai in sorted(limits):
        bar = "#" * max(1, round(34 * limits[ai] / max(limits.values())))
        print(f"    {ai:>7g} ng/day  {limits[ai]:>4}  {bar}")
    print()
    print(f"  at 100 ng/day or tighter    : {len(strict_rows)}")
    print("    by confidence:", ", ".join(
        f"{k} {v}" for k, v in Counter(r.get("confidence") for r in strict_rows).most_common()))
    flagged = sum(1 for r in strict_rows
                  if any(f.get("uncertain") for f in r.get("factors", [])))
    print(f"    containing a [needs review] factor: {flagged}")
    print()
    print(f"  unresolved / conflict       : {conf['unresolved'] + conf['conflict']}")
    print(f"\n  wrote {len(selected)} rows to {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
