"""
Pull the structure images out of the workbook and label them by compound.

Column B of the workbook shows #VALUE! in every row, which looks like the structures are
missing. They are not: 318 PNGs are embedded in the file. Excel stores pictures
separately from cell values, and the cell formula that was meant to display them is what
is broken, not the images.

Each image is anchored to a specific cell, so it can be tied back to the impurity on that
row. That matters because a picture shows exactly which nitrogen carries the nitroso
group, which is precisely what an ambiguous vendor name does not tell us.

    python tools/extract_workbook_images.py

Writes data/structure_images/<catalogue_no>__<slug>.png and an index JSON.
"""

from __future__ import annotations

import json
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import openpyxl

ROOT = Path(__file__).parent.parent
WB = ROOT / "data" / "WEEK1_Nitroso_Impurities_MA.xlsm"
OUT = ROOT / "data" / "structure_images"
INDEX = ROOT / "data" / "structure_images_index.json"

NS = {
    "xdr": "http://schemas.openxmlformats.org/drawingml/2006/spreadsheetDrawing",
    "a": "http://schemas.openxmlformats.org/drawingml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "rel": "http://schemas.openxmlformats.org/package/2006/relationships",
}
EMBED = f"{{{NS['r']}}}embed"


def slug(text: str, limit: int = 48) -> str:
    text = re.sub(r"[^A-Za-z0-9]+", "-", text).strip("-")
    return text[:limit] or "unnamed"


def rels_for(zf: zipfile.ZipFile, part: str) -> dict[str, str]:
    """{rId: target path} for a given part."""
    p = Path(part)
    rels_path = f"{p.parent.as_posix()}/_rels/{p.name}.rels"
    if rels_path not in zf.namelist():
        return {}
    root = ET.fromstring(zf.read(rels_path))
    out = {}
    for rel in root.findall("rel:Relationship", NS):
        target = rel.get("Target", "")
        # Targets are relative to the part's folder.
        resolved = (p.parent / target).as_posix().replace("/./", "/")
        while "/../" in resolved:
            resolved = re.sub(r"[^/]+/\.\./", "", resolved, count=1)
        out[rel.get("Id", "")] = resolved
    return out


def sheet_drawings(zf: zipfile.ZipFile) -> dict[str, str]:
    """{worksheet part: drawing part}"""
    out = {}
    for name in zf.namelist():
        if not re.fullmatch(r"xl/worksheets/sheet\d+\.xml", name):
            continue
        for rid, target in rels_for(zf, name).items():
            if "/drawings/drawing" in target:
                out[name] = target
    return out


def sheet_order(zf: zipfile.ZipFile) -> list[str]:
    """Worksheet parts in workbook order, so they line up with openpyxl's sheetnames."""
    book_rels = rels_for(zf, "xl/workbook.xml")
    root = ET.fromstring(zf.read("xl/workbook.xml"))
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    parts = []
    for sheet in root.find("m:sheets", ns):
        rid = sheet.get(f"{{{NS['r']}}}id")
        parts.append(book_rels.get(rid, ""))
    return parts


def anchors(zf: zipfile.ZipFile, drawing_part: str) -> list[tuple[int, int, str]]:
    """(row0, col0, media path) for every picture in this drawing."""
    rels = rels_for(zf, drawing_part)
    root = ET.fromstring(zf.read(drawing_part))
    out = []
    for anchor in root:
        frm = anchor.find("xdr:from", NS)
        pic = anchor.find("xdr:pic", NS)
        if frm is None or pic is None:
            continue
        blip = pic.find(".//a:blip", NS)
        if blip is None:
            continue
        target = rels.get(blip.get(EMBED, ""), "")
        if not target:
            continue
        out.append((int(frm.find("xdr:row", NS).text),
                    int(frm.find("xdr:col", NS).text),
                    target))
    return out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    wb = openpyxl.load_workbook(WB, data_only=True, keep_vba=True)

    with zipfile.ZipFile(WB) as zf:
        parts = sheet_order(zf)
        drawings = sheet_drawings(zf)
        total_media = sum(1 for n in zf.namelist() if n.startswith("xl/media/"))
        print(f"  {total_media} images embedded in the workbook")
        print(f"  {len(drawings)} sheets carry a drawing\n")

        index, saved, unmatched = [], 0, 0

        for part, sheet_name in zip(parts, wb.sheetnames):
            drawing = drawings.get(part)
            if not drawing:
                continue
            ws = wb[sheet_name]
            found = anchors(zf, drawing)
            if not found:
                continue
            print(f"  sheet {sheet_name!r}: {len(found)} anchored images")

            # Sheets label these columns differently, and the source workbook has
            # typos in its headers ("Nitosoamine impuroty name"), so match loosely.
            headers = {}
            for row in ws.iter_rows(min_row=1, max_row=6):
                for cell in row:
                    if not isinstance(cell.value, str):
                        continue
                    low = cell.value.strip().lower()
                    is_parent = ("parent" in low or low == "source"
                                 or low.startswith("source"))
                    if is_parent:
                        headers.setdefault("parent", cell.column)
                    elif "catalog" in low:
                        headers.setdefault("cat", cell.column)
                    elif "name" in low and "drug" not in low:
                        headers.setdefault("name", cell.column)
                    elif "potency" in low and "categ" in low:
                        headers.setdefault("category", cell.column)
                    elif "ai limit" in low or ("recommended" in low and "ai" in low):
                        headers.setdefault("ai", cell.column)
            if not headers:
                print(f"    (no recognisable headers on {sheet_name!r}, skipping)")
                continue

            def cell_text(row: int, key: str) -> str:
                if key not in headers:
                    return ""
                return str(ws.cell(row, headers[key]).value or "").strip()

            for row0, col0, media in sorted(found):
                excel_row = row0 + 1
                name = cell_text(excel_row, "name")
                cat = cell_text(excel_row, "cat")
                parent = cell_text(excel_row, "parent")

                if not name and not cat:
                    unmatched += 1
                    continue

                filename = f"{slug(cat) or slug(sheet_name)[:6]}__{slug(name)}.png"
                (OUT / filename).write_bytes(zf.read(media))
                index.append({
                    "file": filename, "sheet": sheet_name, "excel_row": excel_row,
                    "impurity_name": name, "catalogue_no": cat, "parent_drug": parent,
                    "sheet_potency_category": cell_text(excel_row, "category"),
                    "sheet_ai_limit": cell_text(excel_row, "ai"),
                    "source_media": media,
                })
                saved += 1

    INDEX.write_text(json.dumps(index, indent=1), encoding="utf-8")

    print(f"\n{'=' * 64}\n  RESULT\n{'=' * 64}")
    print(f"    images matched to a compound : {saved}")
    print(f"    anchored but no row data     : {unmatched}")
    print(f"    saved to                     : {OUT.relative_to(ROOT)}")
    print(f"    index                        : {INDEX.relative_to(ROOT)}")

    if index:
        print("\n  examples:")
        for e in index[:6]:
            print(f"    row {e['excel_row']:>4}  {e['catalogue_no']:<12} "
                  f"{e['impurity_name'][:46]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
