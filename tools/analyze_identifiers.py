"""
What identifiers do we actually have for the catalogue compounds?

Structure lookup strategy depends on this: CAS numbers resolve reliably, vendor
impurity names much less so.
"""

import re
from collections import Counter
from pathlib import Path

import openpyxl

WB = Path(__file__).parent.parent / "data" / "WEEK1_Nitroso_Impurities_MA.xlsm"
CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")


def clean(v) -> str:
    return "" if v is None else str(v).strip()


def main() -> None:
    ws = openpyxl.load_workbook(WB, data_only=True, keep_vba=True)["Working Sheet"]

    rows = []
    for row in ws.iter_rows(min_row=5, max_row=ws.max_row):
        c = {cell.column_letter: clean(cell.value) for cell in row}
        if not c.get("F"):
            continue
        rows.append({
            "parent": c.get("D", ""),
            "cas": c.get("E", ""),
            "name": c.get("F", ""),
            "cat": c.get("G", ""),
        })

    total = len(rows)
    valid_cas = [r for r in rows if CAS_RE.match(r["cas"])]
    no_cas = [r for r in rows if not CAS_RE.match(r["cas"])]
    has_parent = [r for r in rows if r["parent"] and r["parent"] not in {"-", "N/A"}]

    print(f"  Total impurities              : {total}")
    print(f"  With a valid CAS number       : {len(valid_cas):>4}  ({len(valid_cas)/total:.0%})")
    print(f"  Without                       : {len(no_cas):>4}  ({len(no_cas)/total:.0%})")
    print(f"  With a named parent drug      : {len(has_parent):>4}  ({len(has_parent)/total:.0%})")
    print()

    print("  What the non-CAS cells contain:")
    for value, n in Counter(r["cas"] or "(blank)" for r in no_cas).most_common(5):
        print(f"      {value:<12} {n}")
    print()

    print("  Sample names WITH CAS (should resolve easily):")
    for r in valid_cas[:5]:
        print(f"      {r['cas']:<14} {r['name'][:58]}")
    print()

    print("  Sample names WITHOUT CAS (vendor names - harder):")
    for r in no_cas[:8]:
        print(f"      parent={r['parent'][:18]:<18} {r['name'][:52]}")


if __name__ == "__main__":
    main()
