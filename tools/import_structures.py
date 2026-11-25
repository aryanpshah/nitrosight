"""
Import vendor-supplied structures and re-score the catalogue.

Structures inferred from a parent drug agree with FDA only 47-92% of the time; verified
structures agree 100%. So the single highest-value input this project can receive is a
file of real structures, however it arrives: an export from the catalogue supplier, a
corrected version of the workbook, or a sheet a chemist filled in by hand.

This reads any csv / tsv / xlsx carrying SMILES or InChI (either will do), validates
every structure, merges them into the catalogue as high-confidence, re-scores, and
reports exactly what improved.

    python tools/import_structures.py structures.xlsx
    python tools/import_structures.py structures.csv --dry-run

Columns are detected by header name; nothing needs to be in a particular order. Rows are
matched to the catalogue by catalogue number where present, otherwise by impurity name.

One gotcha worth knowing: InChI strings contain commas, so in a CSV they must be quoted.
Unquoted, the row splits into the wrong columns and the structure is lost. The importer
warns when it sees a truncated InChI. Excel files have no such problem.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

from rdkit import Chem, RDLogger

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from cpca import score_smiles  # noqa: E402
from resolve import FDALimits  # noqa: E402

RDLogger.DisableLog("rdApp.*")

RESULTS = ROOT / "data" / "resolved_structures.json"

# Header keywords, in priority order.
FIELDS = {
    "smiles": ("smiles", "canonical smiles", "isomeric smiles", "structure smiles"),
    "inchi": ("inchi", "std inchi", "standard inchi", "inchi string"),
    "catalogue_no": ("catalogue", "catalog", "cat no", "cat.no", "cat_no",
                     "catalogue no", "product code", "sku"),
    "name": ("impurity name", "name", "compound", "product name", "impurity",
             "description"),
}


def normalise_name(text: str) -> str:
    text = text.lower()
    text = re.sub(r"\(.*?\)", " ", text)
    return re.sub(r"[^a-z0-9]+", "", text)


def read_table(path: Path) -> list[dict[str, str]]:
    suffix = path.suffix.lower()
    if suffix in (".xlsx", ".xlsm"):
        import openpyxl
        ws = openpyxl.load_workbook(path, data_only=True).worksheets[0]
        rows = [[("" if c.value is None else str(c.value).strip()) for c in row]
                for row in ws.iter_rows()]
    else:
        delim = "\t" if suffix in (".tsv", ".tab") else ","
        with path.open(encoding="utf-8-sig", newline="") as fh:
            rows = [[c.strip() for c in row] for row in csv.reader(fh, delimiter=delim)]

    if len(rows) < 2:
        return []

    header = [h.lower().strip() for h in rows[0]]
    mapping: dict[str, int] = {}
    for field, keywords in FIELDS.items():
        for i, h in enumerate(header):
            if any(k == h for k in keywords) or any(k in h for k in keywords):
                mapping.setdefault(field, i)
                break

    if "smiles" not in mapping and "inchi" not in mapping:
        raise SystemExit(
            f"No SMILES or InChI column found. Headers were: {header}\n"
            "Expected a column named something like 'SMILES' or 'InChI'.")

    print(f"  columns detected: "
          + ", ".join(f"{k}='{header[i]}'" for k, i in sorted(mapping.items())))

    out = []
    for row in rows[1:]:
        record = {k: (row[i] if i < len(row) else "") for k, i in mapping.items()}
        if any(record.values()):
            out.append(record)
    return out


def structure_from(record: dict[str, str]) -> tuple[str | None, str]:
    """Canonical SMILES plus a note on where it came from, cross-checking if both given."""
    from_smiles = from_inchi = None

    if record.get("smiles"):
        mol = Chem.MolFromSmiles(record["smiles"])
        if mol is not None:
            from_smiles = Chem.MolToSmiles(mol)
    if record.get("inchi"):
        inchi = record["inchi"].strip()
        if inchi.lower().startswith("inchi="):
            mol = Chem.MolFromInchi(inchi)
            if mol is not None:
                from_inchi = Chem.MolToSmiles(mol)

    if from_smiles and from_inchi:
        if from_smiles == from_inchi:
            return from_smiles, "SMILES and InChI agree"
        return None, "SMILES and InChI describe different molecules; not trusted"
    if from_inchi:
        return from_inchi, "from InChI"
    if from_smiles:
        return from_smiles, "from SMILES"

    raw_inchi = (record.get("inchi") or "").strip()
    if raw_inchi.lower().startswith("inchi=") and "/" in raw_inchi:
        return None, ("InChI present but would not parse; if this came from a CSV, "
                      "check it was quoted (InChI strings contain commas)")
    return None, "no valid structure in this row"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("path", help="csv / tsv / xlsx of vendor-supplied structures")
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would change without writing")
    args = ap.parse_args()

    records = read_table(Path(args.path))
    print(f"  {len(records)} rows read\n")

    rows = json.loads(RESULTS.read_text(encoding="utf-8"))
    by_cat = {(r.get("catalogue_no") or "").strip().upper(): r
              for r in rows if r.get("catalogue_no")}
    by_name = {normalise_name(r["impurity_name"]): r for r in rows}

    fda = FDALimits()
    stats = {"matched": 0, "unmatched": 0, "bad_structure": 0,
             "upgraded": 0, "changed_category": 0, "newly_scored": 0}
    changes = []

    for record in records:
        smiles, note = structure_from(record)
        if smiles is None:
            stats["bad_structure"] += 1
            label = record.get("catalogue_no") or record.get("name") or "(unnamed row)"
            print(f"    ! {label[:38]:<38} {note}")
            continue

        target = None
        cat_no = (record.get("catalogue_no") or "").strip().upper()
        if cat_no and cat_no in by_cat:
            target = by_cat[cat_no]
        elif record.get("name"):
            target = by_name.get(normalise_name(record["name"]))
        if target is None:
            stats["unmatched"] += 1
            continue

        stats["matched"] += 1
        before_cat = target.get("category")
        before_conf = target.get("confidence")

        result = score_smiles(smiles)
        target["smiles"] = smiles
        target["method"] = "vendor_supplied"
        target["confidence"] = "high"
        target["structure_note"] = note
        target.pop("candidates", None)
        target.pop("n_candidates", None)
        target.pop("issue", None)

        if result.applicable:
            site = result.governing_site
            target |= {
                "scored": True,
                "category": result.category,
                "ai_ng_day": result.ai_ng_day,
                "potency_score": site.potency_score if site else None,
                "factors": [{"name": f.name, "score": f.score,
                             "reference": f.reference, "uncertain": f.uncertain}
                            for f in (site.factors if site else [])],
                "triage_exit": site.triage_exit if site else None,
                "warnings": result.warnings,
            }
        else:
            target |= {"scored": False,
                       "issue": result.not_applicable_reason}

        pub = fda.lookup(target["impurity_name"])
        target["final_ai_ng_day"] = (pub.ai_ng_day if pub
                                     else target.get("ai_ng_day"))
        target["final_basis"] = f"fda_{pub.basis}" if pub else "cpca_prediction"

        if before_conf != "high":
            stats["upgraded"] += 1
        if before_cat is None and target.get("category"):
            stats["newly_scored"] += 1
        elif before_cat != target.get("category"):
            stats["changed_category"] += 1
            changes.append((target["impurity_name"], before_cat,
                            target.get("category"), before_conf))

    print("=" * 70)
    print("  IMPORT SUMMARY")
    print("=" * 70)
    for key, label in [("matched", "matched a catalogue row"),
                       ("unmatched", "no catalogue row matched"),
                       ("bad_structure", "row had no usable structure"),
                       ("upgraded", "confidence upgraded to high"),
                       ("newly_scored", "scored for the first time"),
                       ("changed_category", "potency category changed")]:
        print(f"    {label:<34} {stats[key]:>5}")
    print()

    if changes:
        print(f"  {len(changes)} category changes "
              f"(these are where the inferred structure was wrong):")
        for name, before, after, conf in changes[:25]:
            print(f"    {name[:46]:<46} {before} -> {after}   (was {conf})")
        if len(changes) > 25:
            print(f"    ... and {len(changes) - 25} more")
        print()

    if args.dry_run:
        print("  dry run, nothing written.")
    else:
        RESULTS.write_text(json.dumps(rows, indent=1), encoding="utf-8")
        print(f"  updated {RESULTS.relative_to(ROOT)}")
        print("  next: python tools/export_report.py --strict")
        print("        python tools/validate_vs_fda.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
