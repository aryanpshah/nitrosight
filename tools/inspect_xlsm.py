"""Dump the structure and contents of the supplied nitrosamine workbook."""

import sys
from pathlib import Path

import openpyxl

PATH = Path(sys.argv[1]) if len(sys.argv) > 1 else (
    Path(__file__).parent.parent / "data" / "WEEK1_Nitroso_Impurities_MA.xlsm"
)
MAX_ROWS = int(sys.argv[2]) if len(sys.argv) > 2 else 60
MAX_COLS = 14

# Pass a path to inspect any workbook. Without one this defaults to the source
# catalogue, which is not redistributed, so explain that rather than crashing.
if not PATH.exists():
    raise SystemExit(
        f"{PATH.name} not found.\n"
        "Pass the path to a workbook as the first argument. The default is the\n"
        "commercial reference-standard catalogue, which belongs to its publisher\n"
        "and is not redistributed with this repository."
    )


def dump(values_only: bool) -> None:
    wb = openpyxl.load_workbook(PATH, data_only=values_only, keep_vba=True)
    mode = "CACHED VALUES" if values_only else "FORMULAS"
    print(f"\n{'=' * 78}\n  {mode}\n{'=' * 78}")
    print(f"Sheets: {wb.sheetnames}\n")
    for ws in wb.worksheets:
        print(f"\n----- SHEET: {ws.title!r}  "
              f"(dims={ws.dimensions}, max_row={ws.max_row}, max_col={ws.max_column}) -----")
        for r, row in enumerate(ws.iter_rows(max_row=min(ws.max_row, MAX_ROWS),
                                             max_col=min(ws.max_column, MAX_COLS)), start=1):
            cells = []
            for c in row:
                if c.value is None:
                    continue
                text = str(c.value).replace("\n", " ⏎ ").strip()
                if len(text) > 70:
                    text = text[:67] + "..."
                cells.append(f"{c.coordinate}={text}")
            if cells:
                print(f"  r{r:<3} | " + " | ".join(cells))


if __name__ == "__main__":
    dump(values_only=True)
    dump(values_only=False)
