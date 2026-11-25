"""How much of the Veeprho nitrosamine catalogue has been CPCA-scored by hand?"""

from pathlib import Path

import openpyxl

WB = Path(__file__).parent.parent / "data" / "WEEK1_Nitroso_Impurities_MA.xlsm"

wb = openpyxl.load_workbook(WB, data_only=True, keep_vba=True)
ws = wb["Working Sheet"]

# Header at r4: D=Parent Drug, E=CAS, F=Impurity name, G=Catalogue No.
# CPCA answers at I (alpha-H?), J (>1 alpha-H?), K (tertiary alpha-C?)
total = scored_q1 = scored_all3 = 0
drugs, catalogue = set(), set()

for row in ws.iter_rows(min_row=5, max_row=ws.max_row):
    cells = {c.column_letter: c.value for c in row}
    name = cells.get("F")
    if not name:
        continue
    total += 1
    if cells.get("D"):
        drugs.add(str(cells["D"]).strip())
    if cells.get("G"):
        catalogue.add(str(cells["G"]).strip())
    if cells.get("I") not in (None, ""):
        scored_q1 += 1
        if cells.get("J") not in (None, "") and cells.get("K") not in (None, ""):
            scored_all3 += 1

print(f"  Impurities listed            : {total}")
print(f"  Distinct parent drugs        : {len(drugs)}")
print(f"  Distinct catalogue numbers   : {len(catalogue)}")
print()
print(f"  Q1 (alpha-H?) answered       : {scored_q1:>4}  ({scored_q1/total:.0%})")
print(f"  All 3 triage Qs answered     : {scored_all3:>4}  ({scored_all3/total:.0%})")
print(f"  Not yet touched              : {total - scored_q1:>4}  ({(total-scored_q1)/total:.0%})")
print()
print("  Sample parent drugs:", ", ".join(sorted(drugs)[:12]))
