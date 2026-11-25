"""
Compare the scorer against the assessor's hand-scored rows in the source workbook.

The workbook records the three Figure 1 triage answers per impurity:

    I : Does the N-nitrosamine have any hydrogens on its alpha-carbon?
    J : Does it have more than one alpha-hydrogen on one or both sides?
    K : Does it have a tertiary alpha-carbon?

Figure 1 exits early, so the assessor leaves later questions as "0" once an earlier one
answers No. Only questions actually answered Y/N are compared.

    python tools/validate_vs_manual.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import openpyxl
from rdkit import RDLogger

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from cpca import score_smiles  # noqa: E402

RDLogger.DisableLog("rdApp.*")

WB = ROOT / "data" / "WEEK1_Nitroso_Impurities_MA.xlsm"
RESULTS = ROOT / "data" / "resolved_structures.json"
MANUAL = ROOT / "data" / "manual_scores.json"

QUESTIONS = [
    ("I", "gate_any_alpha_h", "any hydrogens on an alpha-carbon?"),
    ("J", "gate_multiple_alpha_h", ">1 alpha-hydrogen on one or both sides?"),
    ("K", "gate_tertiary_alpha_c", "tertiary alpha-carbon?"),
]


def yn(value) -> bool | None:
    """The assessor's answers: Y/y/N/n. '0' or blank means 'not asked' under the early exit."""
    if value is None:
        return None
    s = str(value).strip().lower()
    if s in ("y", "yes"):
        return True
    if s in ("n", "no"):
        return False
    return None


def read_workbook() -> list[dict]:
    """The assessor's answers, straight from the source workbook."""
    ws = openpyxl.load_workbook(WB, data_only=True, keep_vba=True)["Working Sheet"]
    out = []
    for row in ws.iter_rows(min_row=5, max_row=ws.max_row):
        c = {cell.column_letter: cell.value for cell in row}
        name = (c.get("F") or "").strip() if c.get("F") else ""
        if not name:
            continue
        answers = {col: yn(c.get(col)) for col, _, _ in QUESTIONS}
        if all(v is None for v in answers.values()):
            continue
        out.append({"name": name, "answers": answers})
    return out


def main() -> int:
    # The workbook is a commercial catalogue and is not redistributed, so the
    # answers are cached alongside the other derived data. The cache holds only
    # the assessor's own Y/N responses and the compound name, which is what this
    # comparison needs; the catalogue itself stays out of the repository.
    if WB.exists():
        manual_rows = read_workbook()
        MANUAL.write_text(json.dumps(manual_rows, indent=1), encoding="utf-8")
        source = f"read from {WB.name} and cached to {MANUAL.name}"
    elif MANUAL.exists():
        manual_rows = json.loads(MANUAL.read_text(encoding="utf-8"))
        source = f"read from {MANUAL.name}; the source workbook is not distributed"
    else:
        print(f"Neither {WB.name} nor {MANUAL.name} is present, so there is "
              f"nothing to compare against.")
        return 1

    structures = {
        r["impurity_name"]: r.get("smiles")
        for r in json.loads(RESULTS.read_text(encoding="utf-8"))
        if r.get("smiles")
    }
    for r in manual_rows:
        r["smiles"] = structures.get(r["name"])

    have_structure = [r for r in manual_rows if r["smiles"]]

    print("=" * 76)
    print("  COMPARISON WITH THE ASSESSOR'S HAND-SCORED ROWS")
    print("=" * 76)
    print(f"  source: {source}")
    print(f"  rows the assessor answered at least one question on : "
          f"{len(manual_rows)}")
    print(f"  ... for which we resolved a structure               : "
          f"{len(have_structure)}")
    print()

    tallies = {col: {"agree": 0, "differ": 0} for col, _, _ in QUESTIONS}
    conflicts = []

    for r in have_structure:
        result = score_smiles(r["smiles"])
        site = result.governing_site
        if site is None:
            continue
        for col, attr, question in QUESTIONS:
            manual = r["answers"][col]
            if manual is None:
                continue
            ours = getattr(site, attr)
            if ours is None:
                continue
            if manual == ours:
                tallies[col]["agree"] += 1
            else:
                tallies[col]["differ"] += 1
                conflicts.append({
                    "name": r["name"], "question": question,
                    "manual": manual, "ours": ours,
                    "counts": site.alpha_h_counts, "smiles": r["smiles"],
                })

    total_agree = sum(t["agree"] for t in tallies.values())
    total = total_agree + sum(t["differ"] for t in tallies.values())
    if total == 0:
        print("  No comparable answers - the assessor's rows and our resolved structures "
              "do not overlap.")
        return 1

    print("  PER QUESTION")
    for col, _, question in QUESTIONS:
        t = tallies[col]
        n = t["agree"] + t["differ"]
        if n:
            print(f"    {col}  {t['agree']:>3}/{n:<3} ({t['agree']/n:>4.0%})  {question}")
    print()
    print(f"  OVERALL: {total_agree}/{total} answers agree ({total_agree/total:.0%})")
    print()

    if conflicts:
        print("=" * 76)
        print(f"  {len(conflicts)} DISAGREEMENTS")
        print("=" * 76)
        for c in conflicts:
            print(f"  {c['name'][:56]}")
            print(f"     Q: {c['question']}")
            print(f"     assessor said {'Y' if c['manual'] else 'N'}, "
                  f"we compute {'Y' if c['ours'] else 'N'}   "
                  f"(alpha-H counts {c['counts']})")
            print(f"     {c['smiles'][:80]}")
            print()
        print("  NOTE: a disagreement may be the assessor's error, ours, or a structure that")
        print("  resolved to the wrong molecule. Each needs a look at the structure.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
