"""
Refresh the data in Nitrosamine_Potency_Tables.xlsx WITHOUT touching its headers.

export_tables.py rebuilds that workbook from scratch, which destroys hand edits. The
headers in the shipped copy have been edited by hand (columns renamed, the Priority
column removed), so regenerating is not an option. This tool instead opens the existing
workbook, finds each sheet's header row, maps the headers that are actually there onto
the values the catalogue provides, and rewrites only the rows beneath them.

    python tools/update_tables.py --dry-run     # report what would change
    python tools/update_tables.py               # write it

Rules it follows:
  * header rows, sheet titles, sheet order and the explanatory sheet are never touched
  * a column whose header it does not recognise is LEFT ALONE, not guessed at, and is
    named in the report
  * data rows take their formatting from the first existing data row of that sheet
  * a .bak copy is written beside the workbook before anything changes
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from collections import Counter, defaultdict
from copy import copy
from pathlib import Path

import openpyxl
from openpyxl.styles import Font, PatternFill

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from export_tables import (AI, CATEGORY_FILL, CONFIDENCE_LABEL,  # noqa: E402
                           METHOD_LABEL, MEANING, dedupe, derivation)

WORKBOOK = ROOT / "data" / "Nitrosamine_Potency_Tables.xlsx"
RESULTS = ROOT / "data" / "resolved_structures.json"


# One producer per column heading. Headings are matched case- and space-insensitively,
# and several spellings map to the same producer so a renamed column still works.


def _ai(r):
    return r.get("final_ai_ng_day") or r.get("ai_ng_day") or ""


COMPOUND_COLUMNS = {
    "parent drug": lambda r: r.get("parent_drug", ""),
    "impurity name": lambda r: r.get("impurity_name", ""),
    "catalogue no.": lambda r: r.get("catalogue_no", ""),
    "potency category": lambda r: r.get("category") or "",
    "ai limit (ng/day)": _ai,
    "basis for the limit": lambda r: (
        "FDA published limit" if str(r.get("final_basis", "")).startswith("fda")
        else "CPCA prediction"),
    "potency score": lambda r: (r.get("potency_score")
                                if r.get("potency_score") is not None else ""),
    "how the score was reached": derivation,
    "structure confidence": lambda r: CONFIDENCE_LABEL.get(
        r.get("confidence"), r.get("confidence") or ""),
    "how we got the structure": lambda r: METHOD_LABEL.get(
        r.get("method"), r.get("method") or ""),
    "notes": lambda r: r.get("override_note") or r.get("issue") or "",
}


def norm(text) -> str:
    return " ".join(str(text or "").split()).strip().lower()


def find_header_row(ws, wanted: str) -> int | None:
    for row in range(1, min(ws.max_row, 8) + 1):
        if norm(ws.cell(row, 1).value) == wanted:
            return row
    return None


def clear_below(ws, header_row: int) -> None:
    if ws.max_row > header_row:
        ws.delete_rows(header_row + 1, ws.max_row - header_row)


def template_styles(ws, header_row: int, ncols: int) -> list:
    """Formatting of the first existing data row, reused for every new row."""
    return [copy(ws.cell(header_row + 1, c)._style) for c in range(1, ncols + 1)]


def write_compound_sheet(ws, header_row: int, rows: list[dict],
                         report: list[str]) -> int:
    ncols = ws.max_column
    headers = [norm(ws.cell(header_row, c).value) for c in range(1, ncols + 1)]
    styles = template_styles(ws, header_row, ncols)

    unknown = [h for h in headers if h and h not in COMPOUND_COLUMNS]
    for h in unknown:
        report.append(f"    ! {ws.title}: column {h!r} not recognised, left empty")

    clear_below(ws, header_row)

    for i, r in enumerate(rows):
        out_row = header_row + 1 + i
        for col, head in enumerate(headers, start=1):
            producer = COMPOUND_COLUMNS.get(head)
            cell = ws.cell(out_row, col, producer(r) if producer else None)
            cell._style = copy(styles[col - 1])
            if head == "potency category" and r.get("category") in CATEGORY_FILL:
                cell.fill = PatternFill("solid",
                                        fgColor=CATEGORY_FILL[r["category"]])
                cell.font = Font(bold=True)
    return len(rows)


def write_summary(ws, header_row: int, scored: list[dict], rows: list[dict],
                  report: list[str]) -> int:
    counts = Counter(r["category"] for r in scored)
    drugs = defaultdict(set)
    for r in scored:
        drugs[r["category"]].add(r.get("parent_drug") or r["impurity_name"])

    headers = [norm(ws.cell(header_row, c).value) for c in range(1, ws.max_column + 1)]
    styles = template_styles(ws, header_row, ws.max_column)

    # The five category rows are overwritten in place rather than cleared, because a
    # Coverage block sits below them and must survive. Its labels are matched by text
    # so a renamed label is reported rather than silently overwritten.
    coverage = {
        "distinct impurities in the catalogue": len(rows),
        "scored": len(scored),
        "not scored (no usable structure)": sum(1 for r in rows if not r.get("scored")),
        "at 100 ng/day or tighter": sum(1 for r in scored if r["category"] <= 2),
        "at 100 ng/day or lower": sum(1 for r in scored if r["category"] <= 2),
    }
    updated_labels = 0
    for row in range(header_row + len(AI) + 1, ws.max_row + 1):
        label = norm(ws.cell(row, 1).value)
        if label in coverage:
            for col in (2, 3, 4):
                if ws.cell(row, col).value is not None or col == 3:
                    ws.cell(row, col, coverage[label] if col == 3 else None)
            updated_labels += 1
        elif label and label != "coverage":
            report.append(f"    ! Summary: coverage label {label!r} not recognised, "
                          "left as it was")
    if updated_labels:
        report.append(f"    - Summary: {updated_labels} coverage figures refreshed")

    producers = {
        "category": lambda cat: f"Category {cat}",
        "ai limit (ng/day)": lambda cat: AI[cat],
        "impurities": lambda cat: counts.get(cat, 0),
        "parent drugs": lambda cat: len(drugs.get(cat, ())),
        "what the category means": lambda cat: MEANING[cat],
    }
    for h in headers:
        if h and h not in producers:
            report.append(f"    ! Summary: column {h!r} not recognised, left empty")

    for i, cat in enumerate(sorted(AI)):
        out_row = header_row + 1 + i
        for col, head in enumerate(headers, start=1):
            fn = producers.get(head)
            cell = ws.cell(out_row, col, fn(cat) if fn else None)
            cell._style = copy(styles[col - 1])
            if head == "category":
                cell.fill = PatternFill("solid", fgColor=CATEGORY_FILL[cat])
                cell.font = Font(bold=True)
    return len(AI)


def write_by_drug(ws, header_row: int, scored: list[dict],
                  report: list[str]) -> int:
    by_drug: dict[str, list[dict]] = defaultdict(list)
    for r in scored:
        drug = r.get("parent_drug") or "-"
        if drug and drug != "-":
            by_drug[drug].append(r)

    def pack(drug: str, impurities: list[dict]) -> dict:
        worst = min(impurities, key=lambda r: r["category"])
        ais = [_ai(r) for r in impurities]
        strict = sum(1 for a in ais if isinstance(a, (int, float)) and a <= 100)
        return {
            "parent drug": drug,
            "worst category": worst["category"],
            "lowest ai limit (ng/day)": _ai(worst),
            "tightest ai (ng/day)": _ai(worst),
            "impurities at 100 ng/day or lower": strict or "",
            "impurities at 100 ng/day or tighter": strict or "",
            "impurities total": len(impurities),
            "impurity driving it": worst["impurity_name"],
            "categories present": ", ".join(
                str(c) for c in sorted({r["category"] for r in impurities})),
        }

    packed = [pack(d, imps) for d, imps in by_drug.items()]
    packed.sort(key=lambda p: (p["worst category"],
                               -(p["impurities at 100 ng/day or lower"] or 0),
                               p["parent drug"]))

    headers = [norm(ws.cell(header_row, c).value) for c in range(1, ws.max_column + 1)]
    styles = template_styles(ws, header_row, ws.max_column)
    for h in headers:
        if h and h not in packed[0]:
            report.append(f"    ! By parent drug: column {h!r} not recognised, "
                          "left empty")

    clear_below(ws, header_row)

    for i, p in enumerate(packed):
        out_row = header_row + 1 + i
        for col, head in enumerate(headers, start=1):
            cell = ws.cell(out_row, col, p.get(head))
            cell._style = copy(styles[col - 1])
            if head == "worst category" and p["worst category"] in CATEGORY_FILL:
                cell.fill = PatternFill(
                    "solid", fgColor=CATEGORY_FILL[p["worst category"]])
                cell.font = Font(bold=True)
    return len(packed)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if not WORKBOOK.exists():
        raise SystemExit(f"{WORKBOOK} not found. Run tools/export_tables.py once first.")

    rows = dedupe(json.loads(RESULTS.read_text(encoding="utf-8")))
    # Most potent first, with the unscored rows last, which is the existing order.
    rows.sort(key=lambda r: (0 if r.get("scored") else 1,
                             r.get("category") or 99,
                             r.get("impurity_name", "")))
    scored = [r for r in rows if r.get("scored")]

    wb = openpyxl.load_workbook(WORKBOOK)
    report: list[str] = []
    written: list[tuple[str, int, int]] = []

    for ws in wb.worksheets:
        title = ws.title

        # Count the contiguous data block under the header, so the Summary sheet's
        # separate Coverage block is not counted as data rows.
        head = find_header_row(ws, "category") or find_header_row(ws, "parent drug")
        before = 0
        if head:
            row = head + 1
            while row <= ws.max_row and ws.cell(row, 1).value is not None:
                before += 1
                row += 1

        if title == "Summary":
            header_row = find_header_row(ws, "category")
            n = (write_summary(ws, header_row, scored, rows, report)
                 if header_row else None)
        elif title == "By parent drug":
            header_row = find_header_row(ws, "parent drug")
            n = write_by_drug(ws, header_row, scored, report) if header_row else None
        elif title.startswith("Category ") or title == "All compounds":
            header_row = find_header_row(ws, "parent drug")
            if header_row is None:
                n = None
            elif title == "All compounds":
                # Every distinct impurity, including the ones with no usable
                # structure, which is how this sheet has always been built.
                n = write_compound_sheet(ws, header_row, rows, report)
            else:
                cat = int(title.split()[1])
                subset = [r for r in scored if r["category"] == cat]
                n = write_compound_sheet(ws, header_row, subset, report)
        else:
            report.append(f"    - {title}: left untouched")
            continue

        if n is None:
            report.append(f"    ! {title}: header row not found, left untouched")
        else:
            written.append((title, before, n))

    print("=" * 72)
    print("  UPDATE Nitrosamine_Potency_Tables.xlsx  (data rows only)")
    print("=" * 72)
    for title, before, after in written:
        delta = after - before
        sign = f"  ({delta:+d})" if delta else ""
        print(f"    {title:<28} {before:>4} -> {after:>4} rows{sign}")
    if report:
        print()
        for line in report:
            print(line)

    # Two different counts, because an FDA published limit can override the predicted
    # category in either direction. The Coverage block uses the category-based count,
    # which is how this workbook has always reported it.
    by_category = sum(1 for r in scored if r["category"] <= 2)
    by_limit = sum(1 for r in scored
                   if isinstance(_ai(r), (int, float)) and _ai(r) <= 100)
    print(f"\n    {len(scored)} scored impurities")
    print(f"      {by_category:>4} in Potency Category 1 or 2   "
          f"(the figure in the Coverage block)")
    print(f"      {by_limit:>4} whose governing limit is 100 ng/day or tighter, "
          f"after any FDA published limit is applied")

    if args.dry_run:
        print("\n  --dry-run: nothing written")
        return 0

    backup = WORKBOOK.with_suffix(".xlsx.bak")
    shutil.copy2(WORKBOOK, backup)
    wb.save(WORKBOOK)
    print(f"\n  backup  -> {backup.relative_to(ROOT)}")
    print(f"  written -> {WORKBOOK.relative_to(ROOT)}")
    print("  headers, sheet titles and the explanatory sheet were not modified.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
