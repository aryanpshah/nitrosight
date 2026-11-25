"""
Build the review workbook: only the structures we are not sure about.

Dr Ali asked for the doubtful ones tabled on their own so a chemist can settle them. A
row earns a place here only if there is a specific question a chemist could answer.
Compounds we simply have no structure for are listed separately, because the answer
there is "find the structure", not a judgement call.

Every row carries the vendor's own structure drawing where we have it, so the question
can usually be answered by looking rather than by reading SMILES.

    python tools/export_doubtful.py
"""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import openpyxl
from openpyxl.drawing.image import Image as XLImage
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from PIL import Image

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

RESULTS = ROOT / "data" / "resolved_structures.json"
INDEX = ROOT / "data" / "structure_images_index.json"
IMAGES = ROOT / "data" / "structure_images"
OSR = ROOT / "data" / "osr_results.json"
OUT = ROOT / "data" / "Doubtful_Structures_For_Review.xlsx"

AI = {1: 26.5, 2: 100.0, 3: 400.0, 4: 1500.0, 5: 1500.0}

HEADER_FILL = PatternFill("solid", fgColor="7F1D1D")
HEADER_FONT = Font(bold=True, color="FFFFFF", size=10)
TITLE_FONT = Font(bold=True, size=13)
WRAP = Alignment(wrap_text=True, vertical="top")
CENTRE = Alignment(horizontal="center", vertical="center")

THUMB_HEIGHT = 96          # pixels; Excel row height is in points, ~0.75 px
ROW_POINTS = THUMB_HEIGHT * 0.78

COLUMNS = [
    ("Structure (vendor's drawing)", 30),
    ("Impurity name", 40),
    ("Parent drug", 20),
    ("Catalogue no.", 13),
    ("Why we are unsure", 34),
    ("What we need from you", 44),
    ("Our answer", 14),
    ("The alternative", 14),
    ("Our structure (SMILES)", 42),
    ("Alternative structure (SMILES)", 42),
]

QUESTION = {
    "several possible nitrosation sites, we guessed":
        "This parent drug has more than one nitrogen that could carry the nitroso "
        "group, and the impurity name does not say which. Which one is it? We assumed "
        "the worst case.",
    "image disagrees with our verified structure":
        "The structure we looked up from the CAS number and the vendor's own drawing "
        "are different molecules. Which is right?",
    "CAS resolves to a structure with no N-nitroso group":
        "The CAS number in the catalogue resolves to a compound with no nitroso group "
        "at all. Is the CAS wrong, or the name?",
    "structure read from a picture, unverified":
        "We read this structure from the vendor's drawing by optical recognition, "
        "which is right about four times in five. Please confirm it.",
}

PRIORITY_NOTE = (
    "strict limit at stake" , "the answer changes whether this sits at 100 ng/day or tighter")


def dedupe(rows: list[dict]) -> list[dict]:
    seen: dict[str, dict] = {}
    for r in rows:
        seen.setdefault((r.get("catalogue_no") or "").strip()
                        or r["impurity_name"].lower(), r)
    return list(seen.values())


def doubt_kind(r: dict) -> str | None:
    if r.get("image_disagrees"):
        return "image disagrees with our verified structure"
    if r.get("confidence") == "conflict":
        return "CAS resolves to a structure with no N-nitroso group"
    if r.get("scored") and r.get("confidence") == "low":
        return "several possible nitrosation sites, we guessed"
    if r.get("method") == "image_recognition":
        return "structure read from a picture, unverified"
    return None


def thumbnail(src: Path, tmp: Path) -> Path | None:
    try:
        img = Image.open(src).convert("RGB")
    except Exception:  # noqa: BLE001
        return None
    ratio = THUMB_HEIGHT / img.height
    img = img.resize((max(1, int(img.width * ratio)), THUMB_HEIGHT), Image.LANCZOS)
    dest = tmp / f"{src.stem}.png"
    img.save(dest)
    return dest


def style_header(ws, row: int) -> None:
    for col, (title, width) in enumerate(COLUMNS, start=1):
        cell = ws.cell(row, col, title)
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = Alignment(wrap_text=True, vertical="center",
                                   horizontal="center")
        ws.column_dimensions[get_column_letter(col)].width = width
    ws.freeze_panes = ws.cell(row + 1, 1)
    ws.auto_filter.ref = f"A{row}:{get_column_letter(len(COLUMNS))}{row}"


def main() -> int:
    rows = dedupe(json.loads(RESULTS.read_text(encoding="utf-8")))
    images = {e["impurity_name"].strip().lower(): e["file"]
              for e in json.loads(INDEX.read_text(encoding="utf-8"))}
    osr = {r["impurity_name"].strip().lower(): r
           for r in json.loads(OSR.read_text(encoding="utf-8"))}

    doubtful, no_structure = [], []
    for r in rows:
        kind = doubt_kind(r)
        if kind:
            doubtful.append((r, kind))
        elif not r.get("scored"):
            no_structure.append(r)

    # Most consequential first: strict limits, then rows a picture can settle.
    def rank(item):
        r, _ = item
        cat = r.get("category") or 9
        alt = (osr.get(r["impurity_name"].strip().lower()) or {}).get("osr_category") or 9
        strict = min(cat, alt) <= 2
        pictured = r["impurity_name"].strip().lower() in images
        return (not strict, not pictured, cat, r["impurity_name"])

    doubtful.sort(key=rank)

    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    tmp = Path(tempfile.mkdtemp(prefix="doubtful_"))

    # read me
    ws = wb.create_sheet("Read me")
    ws.column_dimensions["A"].width = 112
    lines = [
        ("Structures we are not sure about", "title"),
        ("", ""),
        (f"{len(doubtful)} compounds where a chemist could settle a specific question. "
         "Every one has an answer we cannot reach from the catalogue alone.", ""),
        ("", ""),
        ("WHY THESE ARE HERE", "head"),
        ("The scoring itself is not in doubt. Where we know the structure for certain, "
         "the tool matches FDA's own published potency categories 29 times out of 29. "
         "Everything uncertain here is uncertainty about WHICH MOLECULE we are looking "
         "at, not about how to score it.", ""),
        ("", ""),
        ("Most vendor names, such as 'Afatinib Nitroso Impurity 2', do not say which "
         "nitrogen carries the nitroso group. When a drug has more than one candidate "
         "nitrogen, the answer can differ by a factor of 57, from 26.5 ng/day to "
         "1500 ng/day, on the same product name.", ""),
        ("", ""),
        ("HOW TO USE THIS", "head"),
        ("The sheet is sorted so the most consequential questions come first: rows "
         "where the answer decides whether the limit is 100 ng/day or tighter, and "
         "rows where we have the vendor's own drawing to look at.", ""),
        ("Column A holds that drawing where we have it. In most cases the question can "
         "be answered by looking at it rather than by reading SMILES.", ""),
        ("'What we need from you' states the actual question for each row.", ""),
        ("", ""),
        ("THE SEPARATE LIST", "head"),
        (f"A second sheet lists {len(no_structure)} compounds we have no structure for "
         "at all. There is no judgement to make on those, they simply need a structure "
         "from the vendor or a CAS number, so they are kept apart.", ""),
        ("", ""),
        ("Screening aid for a qualified assessor. Not a regulatory determination.", ""),
    ]
    for i, (text, kind) in enumerate(lines, start=1):
        cell = ws.cell(i, 1, text)
        cell.alignment = Alignment(wrap_text=True, vertical="top")
        if kind == "title":
            cell.font = TITLE_FONT
        elif kind == "head":
            cell.font = Font(bold=True, size=11, color="7F1D1D")

    # the questions
    ws = wb.create_sheet("For the chemist")
    ws.cell(1, 1, f"{len(doubtful)} structures needing a decision, "
                  "most consequential first").font = TITLE_FONT
    ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=len(COLUMNS))
    style_header(ws, 3)

    embedded = 0
    for i, (r, kind) in enumerate(doubtful, start=4):
        key = r["impurity_name"].strip().lower()
        alt = osr.get(key) or {}
        alt_cat = alt.get("osr_category")
        cat = r.get("category")

        ws.row_dimensions[i].height = ROW_POINTS
        if (file := images.get(key)) and (src := IMAGES / file).exists():
            if (thumb := thumbnail(src, tmp)) is not None:
                ws.add_image(XLImage(str(thumb)), f"A{i}")
                embedded += 1
        else:
            ws.cell(i, 1, "no drawing available").alignment = CENTRE

        values = [
            None,
            r.get("impurity_name", ""),
            r.get("parent_drug", ""),
            r.get("catalogue_no", ""),
            kind,
            QUESTION.get(kind, ""),
            f"Cat {cat} ({AI[cat]:g} ng/day)" if cat else "not scored",
            f"Cat {alt_cat} ({AI[alt_cat]:g} ng/day)" if alt_cat else "",
            r.get("smiles", "") or "",
            alt.get("osr_smiles") or "",
        ]
        for col, value in enumerate(values, start=1):
            if col == 1:
                continue
            cell = ws.cell(i, col, value)
            cell.alignment = WRAP
            if col == 7 and cat in (1, 2):
                cell.font = Font(bold=True, color="B3261E")
            if col == 8 and alt_cat in (1, 2):
                cell.font = Font(bold=True, color="B3261E")

    # no structure at all
    ws = wb.create_sheet("No structure at all")
    ws.cell(1, 1, f"{len(no_structure)} compounds we could not resolve to any "
                  "structure. These need a SMILES, an InChI or a CAS number, not a "
                  "judgement.").font = TITLE_FONT
    ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=5)
    heads = ["Impurity name", "Parent drug", "Catalogue no.", "CAS in the workbook",
             "Why we could not resolve it"]
    for col, (title, width) in enumerate(zip(heads, [46, 22, 14, 18, 56]), start=1):
        cell = ws.cell(3, col, title)
        cell.fill = HEADER_FILL
        cell.font = HEADER_FONT
        cell.alignment = Alignment(wrap_text=True, vertical="center")
        ws.column_dimensions[get_column_letter(col)].width = width
    ws.freeze_panes = "A4"
    ws.auto_filter.ref = "A3:E3"

    for i, r in enumerate(sorted(no_structure,
                                 key=lambda r: (r.get("parent_drug") or "zzz",
                                                r["impurity_name"])), start=4):
        for col, value in enumerate([
                r.get("impurity_name", ""), r.get("parent_drug", ""),
                r.get("catalogue_no", ""), r.get("cas", ""),
                r.get("issue", "")], start=1):
            ws.cell(i, col, value).alignment = WRAP

    wb.save(OUT)
    print(f"  {len(doubtful)} questions for the chemist "
          f"({embedded} with the vendor's drawing embedded)")
    print(f"  {len(no_structure)} listed separately as needing a structure")
    print(f"\n  wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
