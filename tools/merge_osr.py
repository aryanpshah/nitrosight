"""
Merge image-derived structures into the catalogue, following what the evidence supports.

Adjudicating every disagreement against the workbook's own potency category gave a clear
pattern:

    our confidence    image right    inference right
    high                   0               3
    medium                 2               4
    low                    3               0

So the rule is not "prefer the picture" or "prefer the lookup". It is:

  * high confidence (CAS-verified)  -> keep ours. Recognition errs here, we do not.
  * medium                          -> keep ours, but flag any disagreement for review.
  * low / unresolved / no structure -> take the image. It wins every time we could judge.

One override sits above all of that. Many rows carry a potency category filled in by
someone else on the workbook's 'Nitrosamine impurity' sheet. Where that independent
category backs the image against us, the image is taken whatever our confidence said,
because two sources now agree and we are the outlier.

    python tools/merge_osr.py --dry-run
    python tools/merge_osr.py

ONE NAME, SEVERAL DIFFERENT PICTURES. Records are matched to the catalogue by name, and
two workbook rows can share a name while showing different pictures ("ABCRC-1
Dinitroso Impurity" appears on two sheets, once with the dinitroso compound and once
with a reused mono-nitroso picture). Taking whichever record happens to come first
would make the result depend on sheet order. Instead: if the readings differ, keep only
those carrying the number of N-nitroso groups the NAME implies (groups_from_name, the
same rule fetch_structures.py uses); if exactly one reading survives, use it; otherwise
take no picture, keep our structure and flag the row.

RUN IT ONCE, ON THE OUTPUT OF fetch_structures.py. Taking an image overwrites the row
and marks it "medium" confidence. Run again on that output, the "medium" makes this tool
protect the earlier image reading as if it were our own structure, so a corrected
picture can never replace a wrong one. The tool refuses to run on merged data.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from cpca import score_smiles  # noqa: E402
from resolve import FDALimits, groups_from_name  # noqa: E402

RDLogger.DisableLog("rdApp.*")

OSR = ROOT / "data" / "osr_results.json"
RESULTS = ROOT / "data" / "resolved_structures.json"

TAKE_IMAGE_WHEN = {"low", "conflict", "unresolved", None}
N_NITROSO = Chem.MolFromSmarts("[#7][#7]=[O]")


def nitroso_count(smiles: str) -> int | None:
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    return len(mol.GetSubstructMatches(N_NITROSO)) if mol else None


def pick_reading(name: str, records: list[dict]) -> tuple[dict | None, str]:
    """The one image reading to use for this catalogue name, and how it was chosen.

    Returns (record, "only_reading" | "name_nitroso_count") or (None, "ambiguous").
    """
    distinct = {r["osr_smiles"] for r in records}
    if len(distinct) == 1:
        return records[0], "only_reading"
    want = groups_from_name(name)
    fits = [r for r in records if nitroso_count(r["osr_smiles"]) == want]
    if len({r["osr_smiles"] for r in fits}) == 1:
        return fits[0], "name_nitroso_count"
    return None, "ambiguous"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    osr_rows = json.loads(OSR.read_text(encoding="utf-8"))
    by_name: dict[str, list[dict]] = {}
    for r in osr_rows:
        if r.get("osr_smiles"):
            by_name.setdefault(r["impurity_name"].strip().lower(), []).append(r)

    rows = json.loads(RESULTS.read_text(encoding="utf-8"))
    if any(r.get("method") == "image_recognition" for r in rows):
        raise SystemExit(
            "resolved_structures.json already contains image-derived rows, so it has been "
            "merged before. Run python tools/fetch_structures.py first. "
            "Merging (or comparing) on already-merged data treats an earlier "
            "image result as our own structure and protects it - which is how a "
            "corrected picture would fail to replace a wrong one.")
    fda = FDALimits()

    stats = Counter()
    changes, flagged = [], []

    for row in rows:
        records = by_name.get(row["impurity_name"].strip().lower())
        if not records:
            continue
        stats["has an image"] += 1
        osr, how = pick_reading(row["impurity_name"], records)
        if osr is None:
            row["image_disagrees"] = (
                f"{len(records)} workbook pictures carry this name and read as different "
                "structures; none matches the number of N-nitroso groups the name implies "
                "uniquely, so no picture was used. Worth a look.")
            flagged.append(row)
            stats["several different pictures, none chosen (flagged)"] += 1
            continue
        if how == "name_nitroso_count":
            stats["several different pictures, chosen by the name's nitroso count"] += 1

        confidence = row.get("confidence") if row.get("scored") else None
        take_image = confidence in TAKE_IMAGE_WHEN or not row.get("smiles")

        # Independent tiebreak: the workbook's own category column.
        sheet_cat = str(osr.get("sheet_potency_category") or "").strip()
        corroborated = False
        if sheet_cat.isdigit():
            ref = int(sheet_cat)
            if osr.get("osr_category") == ref and row.get("category") != ref:
                take_image = corroborated = True
                stats["image corroborated by the workbook"] += 1
            elif row.get("category") == ref and osr.get("osr_category") != ref:
                take_image = False   # we are backed up; recognition misread it
                stats["our structure corroborated by the workbook"] += 1

        if not take_image:
            # Keep our structure, but record that the image said something different.
            if osr["verdict"] == "differ":
                row["image_disagrees"] = (
                    f"Optical recognition of the vendor's drawing gives Category "
                    f"{osr.get('osr_category')} rather than {row.get('category')}. "
                    "Our structure is kept because it is CAS-verified. Worth a look.")
                flagged.append(row)
                stats["kept ours, image disagrees (flagged)"] += 1
            else:
                stats["kept ours, image agrees"] += 1
            continue

        before_cat = row.get("category")
        result = score_smiles(osr["osr_smiles"])

        row["smiles"] = osr["osr_smiles"]
        row["method"] = "image_recognition"
        row["confidence"] = "medium"      # a read drawing, never treated as verified
        row["structure_note"] = (
            f"Read from the workbook's structure image ({osr['image']}) by DECIMER "
            "optical recognition."
            + (" The workbook's own potency category backs this reading against our "
               "earlier structure." if corroborated
               else " Not a verified lookup; confirm before relying on it."))
        row.pop("candidates", None)
        row.pop("n_candidates", None)
        row.pop("issue", None)

        if result.applicable:
            site = result.governing_site
            row |= {
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
            stats["took the image"] += 1
            if before_cat is None:
                stats["  ... newly scoreable" ] += 1
            elif before_cat != result.category:
                stats["  ... category changed"] += 1
                changes.append((row["impurity_name"], before_cat, result.category))
        else:
            row |= {"scored": False, "issue": result.not_applicable_reason}
            stats["image read but CPCA does not apply"] += 1

        pub = fda.lookup(row["impurity_name"])
        row["final_ai_ng_day"] = pub.ai_ng_day if pub else row.get("ai_ng_day")
        row["final_basis"] = f"fda_{pub.basis}" if pub else "cpca_prediction"

    print("=" * 72)
    print("  MERGE SUMMARY")
    print("=" * 72)
    for key, n in stats.most_common():
        print(f"    {key:<44} {n:>4}")
    print()

    if changes:
        print(f"  {len(changes)} potency categories changed by the image:")
        for name, before, after in changes[:20]:
            arrow = "stricter" if (after or 9) < (before or 9) else "looser"
            print(f"    {name[:48]:<48} {before} -> {after}  ({arrow})")
        print()

    if flagged:
        print(f"  {len(flagged)} CAS-verified rows where the image disagrees "
              f"(kept ours, flagged for review):")
        for row in flagged[:10]:
            print(f"    {row['impurity_name'][:60]}")
        print()

    if args.dry_run:
        print("  dry run, nothing written.")
    else:
        RESULTS.write_text(json.dumps(rows, indent=1), encoding="utf-8")
        print(f"  updated {RESULTS.relative_to(ROOT)}")
        print("  next: python tools/export_tables.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
