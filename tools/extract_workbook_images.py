"""
Pull the structure images out of the workbook and label them by compound.

Column B of the workbook shows #VALUE! in every row, which looks like the structures are
missing. They are not: 318 PNGs are embedded in the file. Excel stores pictures
separately from cell values, and the cell formula that was meant to display them is what
is broken, not the images.

Each image is anchored to cells, so it can be tied back to the impurity on its row.
That matters because a picture shows exactly which nitrogen carries the nitroso group,
which is precisely what an ambiguous vendor name does not tell us.

WHICH ROW A PICTURE BELONGS TO. A picture is assigned to the row that holds its
VERTICAL CENTRE, computed from the real row heights. An earlier version used the row
holding the picture's TOP-LEFT corner, and that was wrong in a systematic way: pictures
in this workbook routinely start a few points above their own row, so the corner falls
in the row above. Twenty pictures were credited to the previous entry that way, most
with only 1-5% of their height in it (the largest share was 29%; every moved picture
has 71-99% of its height in the row the centre rule picks). Fourteen of them collided
with the previous entry's own picture under the same file name and silently replaced
it - N-Nitroso Midodrine Impurity ended up showing nitroso-migalastat, which is the
picture of the next row - and the next entries lost their pictures altogether.

FILE NAMES ARE NEVER REUSED FOR DIFFERENT CONTENT. When two different pictures would
get the same file name (two rows sharing a name, say), the later one gets a row suffix
instead of overwriting. The content hash of every picture is recorded in the index so
that run_osr.py can tell when the picture behind a file name has changed.

    python tools/extract_workbook_images.py

Writes data/structure_images/<catalogue_no>__<slug>.png and an index JSON.
"""

from __future__ import annotations

import hashlib
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



EMU_PER_POINT = 12700


def centre_row(row_tops: list[float], y_top: float, y_bottom: float) -> int:
    """The 1-based row whose vertical extent contains the midpoint of [y_top, y_bottom].

    `row_tops[r]` is the y position, in points, of the top of excel row r, for
    r = 1..len(row_tops)-1, with the last entry closing the final row. A pure
    function so the placement rule can be tested without the workbook.
    """
    mid = (y_top + y_bottom) / 2
    for r in range(1, len(row_tops) - 1):
        if row_tops[r] <= mid < row_tops[r + 1]:
            return r
    return len(row_tops) - 2


def placed_pictures(zf: zipfile.ZipFile, drawing_part: str, ws) -> list[tuple[int, int, int, str]]:
    """(centre_row, top_left_row, col0, media path) for every picture in this drawing.

    Rows are 1-based excel rows. See the module docstring for why the centre row,
    not the top-left row, is the one a picture belongs to.
    """
    default_h = ws.sheet_format.defaultRowHeight or 15.0
    last = max(ws.max_row, 1) + 200      # pictures can hang below the last data row
    row_tops = [0.0, 0.0]
    for r in range(1, last + 1):
        row_tops.append(row_tops[-1] + (ws.row_dimensions[r].height or default_h))

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
        top_row = int(frm.find("xdr:row", NS).text) + 1
        col0 = int(frm.find("xdr:col", NS).text)
        if top_row >= len(row_tops) - 1:
            out.append((top_row, top_row, col0, target))
            continue
        y_top = row_tops[top_row] + int(frm.find("xdr:rowOff", NS).text) / EMU_PER_POINT
        to = anchor.find("xdr:to", NS)
        if to is not None:
            to_row = min(int(to.find("xdr:row", NS).text) + 1, len(row_tops) - 2)
            y_bottom = row_tops[to_row] + int(to.find("xdr:rowOff", NS).text) / EMU_PER_POINT
        else:                                  # oneCellAnchor: size given by ext
            ext = anchor.find("xdr:ext", NS)
            height = int(ext.get("cy")) / EMU_PER_POINT if ext is not None else 0.0
            y_bottom = y_top + height
        out.append((centre_row(row_tops, y_top, y_bottom), top_row, col0, target))
    return out


def _require_workbook(path):
    """The source workbook is a commercial catalogue and is not redistributed.
    Say so plainly rather than failing with a traceback."""
    if path.exists():
        return
    raise SystemExit(
        f"{path.name} is not in this repository.\n"
        "It is a commercial reference-standard catalogue and belongs to its\n"
        "publisher, so it is not redistributed. Everything derived from it is\n"
        "in data/, so the published results can still be checked. See the\n"
        "README section 'What is not in this repository, and why'."
    )

def main() -> int:
    _require_workbook(WB)
    OUT.mkdir(parents=True, exist_ok=True)
    wb = openpyxl.load_workbook(WB, data_only=True, keep_vba=True)

    with zipfile.ZipFile(WB) as zf:
        parts = sheet_order(zf)
        drawings = sheet_drawings(zf)
        total_media = sum(1 for n in zf.namelist() if n.startswith("xl/media/"))
        print(f"  {total_media} images embedded in the workbook")
        print(f"  {len(drawings)} sheets carry a drawing\n")

        index, saved, unmatched = [], 0, 0
        written: dict[str, str] = {}          # lower-cased file name -> md5 this run
        moved = fallbacks = renamed = 0

        for part, sheet_name in zip(parts, wb.sheetnames):
            drawing = drawings.get(part)
            if not drawing:
                continue
            ws = wb[sheet_name]
            found = placed_pictures(zf, drawing, ws)
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

            for centre, top_left, col0, media in sorted(found):
                excel_row, assigned_by = centre, "centre"
                if not cell_text(excel_row, "name") and not cell_text(excel_row, "cat") \
                        and (cell_text(top_left, "name") or cell_text(top_left, "cat")):
                    # Centred on a row with no data but the corner sits on one:
                    # fall back rather than lose the picture, and count it.
                    excel_row, assigned_by = top_left, "top_left_fallback"
                    fallbacks += 1
                if excel_row != top_left and assigned_by == "centre":
                    moved += 1
                name = cell_text(excel_row, "name")
                cat = cell_text(excel_row, "cat")
                parent = cell_text(excel_row, "parent")

                if not name and not cat:
                    unmatched += 1
                    continue

                data = zf.read(media)
                digest = hashlib.md5(data).hexdigest()
                filename = f"{slug(cat) or slug(sheet_name)[:6]}__{slug(name)}.png"
                # Keyed case-insensitively: on Windows and macOS "labetalol" and
                # "Labetalol" are the SAME file, so a case-only difference is no
                # protection against overwriting.
                if written.get(filename.lower(), digest) != digest:
                    # Same name, DIFFERENT picture: never overwrite. Suffix the row.
                    stem = filename[:-4]
                    filename, n = f"{stem}__r{excel_row}.png", 2
                    while written.get(filename.lower(), digest) != digest:
                        filename, n = f"{stem}__r{excel_row}_{n}.png", n + 1
                    renamed += 1
                written[filename.lower()] = digest
                (OUT / filename).write_bytes(data)
                index.append({
                    "file": filename, "sheet": sheet_name, "excel_row": excel_row,
                    "impurity_name": name, "catalogue_no": cat, "parent_drug": parent,
                    "sheet_potency_category": cell_text(excel_row, "category"),
                    "sheet_ai_limit": cell_text(excel_row, "ai"),
                    "source_media": media,
                    "md5": digest,
                    "assigned_by": assigned_by,
                    "top_left_row": top_left,
                })
                saved += 1

    INDEX.write_text(json.dumps(index, indent=1), encoding="utf-8")

    print(f"\n{'=' * 64}\n  RESULT\n{'=' * 64}")
    print(f"    images matched to a compound : {saved}")
    print(f"    anchored but no row data     : {unmatched}")
    print(f"    placed on a different row than its top-left corner : {moved}")
    print(f"    fell back to the top-left row (centre row empty)  : {fallbacks}")
    print(f"    given a row suffix to avoid overwriting            : {renamed}")
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
