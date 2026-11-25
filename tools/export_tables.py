"""
Build the reviewer-facing Excel workbook: parent drugs grouped by potency category,
with the full derivation of every score.

The point of the "How the score is derived" sheet is traceability. Every number in this
workbook can be walked back to a numbered table in the FDA guidance. No language model
takes part in the calculation.

    python tools/export_tables.py
"""

from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import openpyxl
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from resolve import FDALimits  # noqa: E402

RESULTS = ROOT / "data" / "resolved_structures.json"
OUT = ROOT / "data" / "Nitrosamine_Potency_Tables.xlsx"

AI = {1: 26.5, 2: 100.0, 3: 400.0, 4: 1500.0, 5: 1500.0}
MEANING = {
    1: "Highest potency. Limit equals the class-specific nitrosamine limit (based on NDEA).",
    2: "Potency comparable to NDMA and NNK.",
    3: "Lower potency than Category 2, e.g. a weakly deactivating feature is present.",
    4: "May be activated by alpha-hydroxylation but predicted to be of low potency.",
    5: "Not predicted to be activated by alpha-hydroxylation (no or hindered alpha-hydrogens).",
}
CATEGORY_FILL = {
    1: "F4CCCC", 2: "FCE5CD", 3: "FFF2CC", 4: "D9EAD3", 5: "D9EAD3",
}

HEADER_FILL = PatternFill("solid", fgColor="1F3864")
HEADER_FONT = Font(bold=True, color="FFFFFF", size=10)
TITLE_FONT = Font(bold=True, size=13)
WRAP = Alignment(wrap_text=True, vertical="top")

COLUMNS = [
    ("Parent drug", 24),
    ("Impurity name", 46),
    ("Catalogue no.", 14),
    ("Potency category", 10),
    ("AI limit (ng/day)", 12),
    ("Basis for the limit", 22),
    ("Potency score", 9),
    ("How the score was reached", 56),
    ("Structure confidence", 14),
    ("How we got the structure", 22),
    ("Notes", 40),
]


def dedupe(rows: list[dict]) -> list[dict]:
    seen: dict[str, dict] = {}
    for r in rows:
        seen.setdefault((r.get("catalogue_no") or "").strip()
                        or r["impurity_name"].lower(), r)
    return list(seen.values())


def derivation(r: dict) -> str:
    factors = r.get("factors") or []
    if factors:
        parts = [f"{f['score']:+d} {f['name']} [{f['reference']}]" for f in factors]
        return " ; ".join(parts) + f"  =  {r.get('potency_score')}"
    if r.get("triage_exit"):
        return r["triage_exit"]
    return r.get("issue", "")


METHOD_LABEL = {
    "cas_lookup": "CAS number, verified",
    "parent_nitrosation": "built from parent drug",
    "vendor_supplied": "supplied by vendor",
    "parent_lookup_failed": "could not resolve",
    "none": "no identifier available",
}

CONFIDENCE_LABEL = {
    "high": "High (verified)",
    "medium": "Medium (inferred, one site)",
    "low": "Low (inferred, several possible sites)",
    "conflict": "Conflict (check the CAS)",
    "unresolved": "Unresolved",
}


def style_header(ws, row: int = 1) -> None:
    for col, (title, width) in enumerate(COLUMNS, start=1):
        cell = ws.cell(row, col, title)
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = Alignment(wrap_text=True, vertical="center")
        ws.column_dimensions[get_column_letter(col)].width = width
    ws.freeze_panes = ws.cell(row + 1, 1)
    ws.auto_filter.ref = f"A{row}:{get_column_letter(len(COLUMNS))}{row}"


def write_row(ws, row: int, r: dict) -> None:
    cat = r.get("category")
    values = [
        r.get("parent_drug", ""),
        r.get("impurity_name", ""),
        r.get("catalogue_no", ""),
        cat or "",
        r.get("final_ai_ng_day") or r.get("ai_ng_day") or "",
        "FDA published limit" if str(r.get("final_basis", "")).startswith("fda")
        else "CPCA prediction",
        r.get("potency_score") if r.get("potency_score") is not None else "",
        derivation(r),
        CONFIDENCE_LABEL.get(r.get("confidence"), r.get("confidence") or ""),
        METHOD_LABEL.get(r.get("method"), r.get("method") or ""),
        r.get("override_note") or r.get("issue") or "",
    ]
    for col, value in enumerate(values, start=1):
        cell = ws.cell(row, col, value)
        cell.alignment = WRAP
        if col == 4 and cat in CATEGORY_FILL:
            cell.fill = PatternFill("solid", fgColor=CATEGORY_FILL[cat])
            cell.font = Font(bold=True)


def sheet_readme(wb, rows: list[dict], scored: list[dict], fda: FDALimits) -> None:
    ws = wb.create_sheet("How the score is derived")
    ws.column_dimensions["A"].width = 118

    lines = [
        ("How the score is derived", "title"),
        ("", ""),
        ("No language model takes part in the calculation. Every number below can be "
         "walked back to a numbered table in the FDA guidance.", ""),
        ("", ""),
        ("THE SOURCE", "head"),
        ("FDA, Recommended Acceptable Intake Limits for Nitrosamine Drug "
         "Substance-Related Impurities (NDSRIs), Guidance for Industry, August 2023.", ""),
        ("A copy is in the guidelines/ folder of the project. The scoring rules are "
         "Appendix A (pages 18 to 20) and the decision flowchart is Figure 1 (page 10).", ""),
        ("", ""),
        ("THE CALCULATION", "head"),
        ("1. The structure is read with RDKit, an open-source chemistry toolkit. The "
         "N-nitroso group is located, and the hydrogens on each neighbouring "
         "(alpha) carbon are counted.", ""),
        ("2. Three triage questions from Figure 1 are answered. Any one of them can "
         "send the compound straight to Category 5:", ""),
        ("       a. Are there any hydrogens on the alpha-carbons?", ""),
        ("       b. Does any alpha-carbon carry more than one hydrogen?", ""),
        ("       c. Is there a tertiary alpha-carbon?", ""),
        ("3. Otherwise a Potency Score is summed:", ""),
        ("       Potency Score = alpha-Hydrogen Score (Table A)", ""),
        ("                     + Deactivating Feature Score (Table B, positive)", ""),
        ("                     + Activating Feature Score (Table C, negative)", ""),
        ("4. The score maps to a category and a limit (Figure 1):", ""),
        ("       score 1 or less -> Category 1 -> 26.5 ng/day", ""),
        ("       score 2         -> Category 2 -> 100 ng/day", ""),
        ("       score 3         -> Category 3 -> 400 ng/day", ""),
        ("       score 4 or more -> Category 4 -> 1500 ng/day", ""),
        ("5. Where FDA has already published a limit for that specific nitrosamine, "
         "the published value is used instead of the prediction. CPCA is only for "
         "compounds that lack an established limit.", ""),
        ("", ""),
        ("The 'How the score was reached' column on every sheet lists each "
         "contributing factor with the table it came from, so any row can be "
         "checked by hand against the guidance.", ""),
        ("", ""),
        ("HOW FAR IT HAS BEEN CHECKED", "head"),
        ("The two worked examples in the FDA guidance (N-nitroso-l-nebivolol and "
         "N-nitroso-ketamine) are reproduced exactly, factor for factor.", ""),
        ("FDA publishes its own potency category for about 280 named nitrosamines. "
         "119 of those appear in this catalogue. Agreement:", ""),
        ("       where the structure is verified          29 of 29    (100%)", ""),
        ("       where it was inferred from the parent    67 of 73    (92%)", ""),
        ("       where several sites were possible         8 of 17    (47%)", ""),
        ("", ""),
        ("That pattern is the important one. Where the structure is known, the "
         "implementation matches FDA every time. The error is in knowing which "
         "molecule we are looking at, not in the scoring.", ""),
        ("", ""),
        ("ON THE COHORT OF CONCERN", "head"),
        ("Every compound in these tables is already in it. ICH M7(R2) names "
         "N-nitroso compounds, alongside aflatoxin-like and alkyl-azoxy structures, as "
         "the cohort of concern: structural classes potent enough that the ordinary "
         "threshold of toxicological concern does not protect against them.", ""),
        ("That is why these limits are in nanograms per day rather than the 1.5 "
         "micrograms per day ICH M7 allows for an unstudied impurity. The strictest "
         "category here, 26.5 ng/day, is about 57 times tighter.", ""),
        ("So the question is not whether a compound is cohort of concern. It is how "
         "potent it is within that cohort, which is what the category says.", ""),
        ("", ""),
        ("WHAT TO TREAT WITH CAUTION", "head"),
        ("Rows marked Low confidence were built from the parent drug where more than "
         "one nitrogen could carry the nitroso group. The impurity name does not say "
         "which, so the worst case is reported. These need the real structure.", ""),
        ("Rows marked Conflict have a CAS number that resolves to a structure with no "
         "N-nitroso group at all, so either the CAS or the name is wrong.", ""),
        ("This is a screening aid for a qualified assessor. It is not a regulatory "
         "determination.", ""),
        ("", ""),
        (f"FDA published-limits snapshot taken {fda.fetched}. "
         f"Catalogue compounds: {len(rows)}. Scored: {len(scored)}.", ""),
    ]

    for i, (text, kind) in enumerate(lines, start=1):
        cell = ws.cell(i, 1, text)
        cell.alignment = Alignment(wrap_text=True, vertical="top")
        if kind == "title":
            cell.font = TITLE_FONT
        elif kind == "head":
            cell.font = Font(bold=True, size=11, color="1F3864")


def sheet_summary(wb, scored: list[dict], rows: list[dict]) -> None:
    ws = wb.create_sheet("Summary")
    ws.column_dimensions["A"].width = 16
    ws.column_dimensions["B"].width = 18
    ws.column_dimensions["C"].width = 10
    ws.column_dimensions["D"].width = 14
    ws.column_dimensions["E"].width = 72

    ws.cell(1, 1, "Potency categories across the catalogue").font = TITLE_FONT

    headers = ["Category", "AI limit (ng/day)", "Impurities", "Parent drugs",
               "What the category means"]
    for col, title in enumerate(headers, start=1):
        cell = ws.cell(3, col, title)
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = Alignment(wrap_text=True, vertical="center")

    counts = Counter(r["category"] for r in scored)
    drugs = defaultdict(set)
    for r in scored:
        if r.get("parent_drug") and r["parent_drug"] not in {"-", "N/A"}:
            drugs[r["category"]].add(r["parent_drug"])

    for i, cat in enumerate(sorted(counts), start=4):
        ws.cell(i, 1, f"Category {cat}").font = Font(bold=True)
        ws.cell(i, 1).fill = PatternFill("solid", fgColor=CATEGORY_FILL[cat])
        ws.cell(i, 2, AI[cat])
        ws.cell(i, 3, counts[cat])
        ws.cell(i, 4, len(drugs[cat]))
        ws.cell(i, 5, MEANING[cat]).alignment = WRAP

    row = 4 + len(counts) + 2
    ws.cell(row, 1, "Coverage").font = Font(bold=True, size=11, color="1F3864")
    unresolved = sum(1 for r in rows if not r.get("scored"))
    for offset, (label, value) in enumerate([
        ("Distinct impurities in the catalogue", len(rows)),
        ("Scored", len(scored)),
        ("Not scored (no usable structure)", unresolved),
        ("At 100 ng/day or tighter", sum(1 for r in scored if r["category"] <= 2)),
    ], start=1):
        ws.cell(row + offset, 1, label)
        ws.cell(row + offset, 3, value).font = Font(bold=True)
        ws.merge_cells(start_row=row + offset, start_column=1,
                       end_row=row + offset, end_column=2)


PRIORITY = {
    1: ("Study first", "FF9999"),
    2: ("Study first", "FFCC99"),
    3: ("Worth a look", "FFEB99"),
    4: ("Low priority", "DFF0D8"),
    5: ("Low priority", "DFF0D8"),
}


def sheet_by_drug(wb, scored: list[dict]) -> None:
    """One row per parent drug, ordered for deciding where method development is worth
    the effort: a drug carrying a high-potency impurity earns attention."""
    ws = wb.create_sheet("By parent drug")
    ws.cell(1, 1, "Parent drugs ranked by their most potent impurity. A drug near the "
                  "top carries an impurity needing control at 100 ng/day or tighter, "
                  "so it is worth study and method development.").font = TITLE_FONT
    ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=8)

    headers = ["Priority", "Parent drug", "Worst category", "Tightest AI (ng/day)",
               "Impurities at 100 ng/day or tighter", "Impurities total",
               "Impurity driving it", "Categories present"]
    widths = [14, 28, 12, 16, 15, 12, 46, 16]
    for col, (title, width) in enumerate(zip(headers, widths), start=1):
        cell = ws.cell(3, col, title)
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = Alignment(wrap_text=True, vertical="center")
        ws.column_dimensions[get_column_letter(col)].width = width
    ws.freeze_panes = "A4"
    ws.auto_filter.ref = f"A3:{get_column_letter(len(headers))}3"

    by_drug = defaultdict(list)
    for r in scored:
        drug = r.get("parent_drug") or ""
        if drug and drug not in {"-", "N/A"}:
            by_drug[drug].append(r)

    ordered = sorted(
        by_drug.items(),
        key=lambda kv: (min(r["category"] for r in kv[1]),
                        -sum(1 for r in kv[1] if r["category"] <= 2),
                        kv[0]))

    for i, (drug, impurities) in enumerate(ordered, start=4):
        worst = min(impurities, key=lambda r: r["category"])
        cat = worst["category"]
        cats = sorted({r["category"] for r in impurities})
        strict = sum(1 for r in impurities if r["category"] <= 2)

        label, colour = PRIORITY[cat]
        p = ws.cell(i, 1, label)
        p.fill = PatternFill("solid", fgColor=colour)
        p.font = Font(bold=True)

        ws.cell(i, 2, drug)
        c = ws.cell(i, 3, cat)
        c.fill = PatternFill("solid", fgColor=CATEGORY_FILL[cat])
        c.font = Font(bold=True)
        ws.cell(i, 4, worst.get("final_ai_ng_day") or worst["ai_ng_day"])
        ws.cell(i, 5, strict or "")
        ws.cell(i, 6, len(impurities))
        ws.cell(i, 7, worst["impurity_name"]).alignment = WRAP
        ws.cell(i, 8, ", ".join(str(x) for x in cats))


def main() -> int:
    rows = dedupe(json.loads(RESULTS.read_text(encoding="utf-8")))
    scored = [r for r in rows if r.get("scored")]
    fda = FDALimits()

    wb = openpyxl.Workbook()
    wb.remove(wb.active)

    sheet_readme(wb, rows, scored, fda)
    sheet_summary(wb, scored, rows)
    sheet_by_drug(wb, scored)

    # One sheet per category, most potent first, highest confidence first within each.
    order = {"high": 0, "medium": 1, "low": 2, "conflict": 3, "unresolved": 4}
    for cat in sorted({r["category"] for r in scored}):
        ws = wb.create_sheet(f"Category {cat} ({AI[cat]:g} ng-day)")
        ws.cell(1, 1, f"Potency Category {cat}   |   {AI[cat]:g} ng/day   |   "
                      f"{MEANING[cat]}").font = TITLE_FONT
        ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=len(COLUMNS))
        style_header(ws, row=3)
        members = sorted((r for r in scored if r["category"] == cat),
                         key=lambda r: (order.get(r.get("confidence"), 9),
                                        r.get("parent_drug") or "zzz",
                                        r["impurity_name"]))
        for i, r in enumerate(members, start=4):
            write_row(ws, i, r)

    # Everything, including rows we could not score.
    ws = wb.create_sheet("All compounds")
    style_header(ws, row=1)
    for i, r in enumerate(sorted(rows, key=lambda r: (r.get("category") or 99,
                                                      r.get("parent_drug") or "zzz")),
                          start=2):
        write_row(ws, i, r)

    wb.save(OUT)

    print(f"  sheets: {', '.join(wb.sheetnames)}")
    print(f"\n  wrote {OUT.relative_to(ROOT)}")
    print(f"  {len(rows)} compounds, {len(scored)} scored")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
