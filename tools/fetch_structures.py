"""
Resolve every catalogue row to a structure, then score it.

Two tiers:

  1. CAS number (27% of rows) -> look it up directly on PubChem, then VERIFY the
     returned structure actually contains an N-nitroso group. A CAS that resolves to
     something without one is a data problem, not a result.

  2. No CAS (73%) -> look up the PARENT drug and build the NDSRI by nitrosating its
     amines. One nitrosatable site means one unambiguous product; several means several
     candidates, and we report the worst-case category.

Everything is cached, so re-running is cheap and an interrupted run resumes.

    python tools/fetch_structures.py            # whole catalogue
    python tools/fetch_structures.py --limit 40 # quick trial
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path

import openpyxl
from rdkit import Chem, RDLogger

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from cpca import score_smiles  # noqa: E402
from resolve import (  # noqa: E402
    FDALimits, PubChemClient, enumerate_ndsris, groups_from_name,
)

RDLogger.DisableLog("rdApp.*")

WB_PATH = ROOT / "data" / "WEEK1_Nitroso_Impurities_MA.xlsm"
OUT_PATH = ROOT / "data" / "resolved_structures.json"
CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")
NITROSO = Chem.MolFromSmarts("[NX3]-[NX2]=[OX1]")
SAVE_EVERY = 25


def clean(v) -> str:
    return "" if v is None else str(v).strip()


def has_nitroso(smiles: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None and mol.HasSubstructMatch(NITROSO)


def read_rows(limit: int | None) -> list[dict]:
    ws = openpyxl.load_workbook(WB_PATH, data_only=True, keep_vba=True)["Working Sheet"]
    rows = []
    for row in ws.iter_rows(min_row=5, max_row=ws.max_row):
        c = {cell.column_letter: clean(cell.value) for cell in row}
        if not c.get("F"):
            continue
        rows.append({
            "parent_drug": c.get("D", ""),
            "cas": c.get("E", ""),
            "impurity_name": c.get("F", ""),
            "catalogue_no": c.get("G", ""),
        })
        if limit and len(rows) >= limit:
            break
    return rows


def score_best(smiles_list: list[str]) -> dict:
    """Score candidates and keep the most potent (lowest category number)."""
    best, best_smiles = None, None
    for smi in smiles_list:
        r = score_smiles(smi)
        if not r.applicable:
            continue
        if best is None or r.category < best.category:
            best, best_smiles = r, smi
    if best is None:
        return {"scored": False}
    site = best.governing_site
    return {
        "scored": True,
        "scored_smiles": best_smiles,
        "category": best.category,
        "ai_ng_day": best.ai_ng_day,
        "potency_score": site.potency_score if site else None,
        "factors": [
            {"name": f.name, "score": f.score, "reference": f.reference,
             "uncertain": f.uncertain}
            for f in (site.factors if site else [])
        ],
        "triage_exit": site.triage_exit if site else None,
        "warnings": best.warnings,
    }


def resolve_row(row: dict, pc: PubChemClient) -> dict:
    out = dict(row)

    # Tier 1: CAS
    if CAS_RE.match(row["cas"]):
        hit = pc.lookup(row["cas"])
        if hit.found:
            if has_nitroso(hit.smiles):
                out |= {"smiles": hit.smiles, "method": "cas_lookup",
                        "confidence": "high", "pubchem_cid": hit.cid}
                out |= score_best([hit.smiles])
                return out
            out |= {"smiles": hit.smiles, "method": "cas_lookup",
                    "confidence": "conflict", "pubchem_cid": hit.cid,
                    "issue": "CAS resolved to a structure with no N-nitroso group; "
                             "the catalogue CAS may be for the non-nitrosated impurity",
                    "scored": False}
            return out

    # Tier 2: build from the parent drug
    parent = row["parent_drug"]
    if not parent or parent in {"-", "N/A"}:
        out |= {"method": "none", "confidence": "unresolved",
                "issue": "no CAS and no parent drug named", "scored": False}
        return out

    hit = pc.lookup(parent)
    if not hit.found:
        out |= {"method": "parent_lookup_failed", "confidence": "unresolved",
                "issue": f"parent drug '{parent}' not found on PubChem", "scored": False}
        return out

    n_groups = groups_from_name(row["impurity_name"])
    cands = enumerate_ndsris(hit.smiles, n_groups)
    if not cands:
        out |= {"method": "parent_nitrosation", "confidence": "unresolved",
                "parent_smiles": hit.smiles, "parent_cid": hit.cid,
                "issue": f"parent has no site able to carry {n_groups} N-nitroso "
                         "group(s) (no suitable secondary amine)", "scored": False}
        return out

    unambiguous = len(cands) == 1
    out |= {
        "method": "parent_nitrosation",
        "confidence": "medium" if unambiguous else "low",
        "parent_smiles": hit.smiles,
        "parent_cid": hit.cid,
        "n_nitroso_groups": n_groups,
        "n_candidates": len(cands),
        "routes": sorted({c.route for c in cands}),
        "candidates": [{"smiles": c.smiles, "site": c.site_description,
                        "route": c.route} for c in cands],
    }
    if not unambiguous:
        out["issue"] = (f"{len(cands)} possible nitrosation sites; cannot tell which one "
                        "the vendor's name refers to. Worst case reported.")
    out |= score_best([c.smiles for c in cands])
    out["smiles"] = out.get("scored_smiles")
    return out


def apply_published_limits(results: list[dict]) -> None:
    """Where FDA has already published a limit for this nitrosamine, that value governs.

    CPCA is for compounds that LACK an established limit - the guidance is explicit that
    compound-specific data takes precedence. Both numbers are kept so the prediction can
    still be audited against the published one.
    """
    fda = FDALimits()
    print(f"\n  FDA published limits: {len(fda)} named entries "
          f"(snapshot {fda.fetched})")
    if warning := fda.staleness_warning():
        print(f"  ! {warning}")

    overridden = agreed = 0
    for r in results:
        pub = fda.lookup(r["impurity_name"])
        if pub is None:
            # No published limit: the CPCA prediction stands.
            if r.get("scored"):
                r["final_ai_ng_day"] = r["ai_ng_day"]
                r["final_basis"] = "cpca_prediction"
            continue

        r["published_name"] = pub.name
        r["published_ai_ng_day"] = pub.ai_ng_day
        r["published_category"] = pub.category
        r["published_basis"] = pub.basis
        r["final_ai_ng_day"] = pub.ai_ng_day
        r["final_basis"] = f"fda_{pub.basis}"
        overridden += 1
        if r.get("scored") and pub.category is not None:
            if r["category"] == pub.category:
                agreed += 1
            else:
                r["override_note"] = (
                    f"CPCA predicted Category {r['category']} "
                    f"({r['ai_ng_day']:g} ng/day); FDA publishes Category "
                    f"{pub.category} ({pub.ai_ng_day:g} ng/day). Published value used.")
        if pub.is_interim:
            r.setdefault("issue", "")
            r["issue"] = (r["issue"] + " " if r["issue"] else "") + \
                f"FDA limit is INTERIM ({pub.duration or 'time-limited'})."

    print(f"  matched a published limit : {overridden}")
    print(f"  ... of those, CPCA agreed : {agreed}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, help="only process the first N rows")
    args = ap.parse_args()

    rows = read_rows(args.limit)
    pc = PubChemClient()
    print(f"  resolving {len(rows)} catalogue rows\n")

    results = []
    for i, row in enumerate(rows, start=1):
        results.append(resolve_row(row, pc))
        if i % SAVE_EVERY == 0:
            pc.save()
            print(f"    {i}/{len(rows)}  ({pc.live_requests} live requests so far)")

    pc.save()
    apply_published_limits(results)
    OUT_PATH.write_text(json.dumps(results, indent=1), encoding="utf-8")

    # summary
    conf = Counter(r.get("confidence", "?") for r in results)
    scored = [r for r in results if r.get("scored")]
    cats = Counter(r["category"] for r in scored)

    print(f"\n{'=' * 64}\n  RESOLUTION\n{'=' * 64}")
    for label in ("high", "medium", "low", "conflict", "unresolved"):
        if conf[label]:
            print(f"    {label:<12} {conf[label]:>5}  ({conf[label]/len(results):.0%})")
    print(f"    {'-' * 28}\n    {'scored':<12} {len(scored):>5}  "
          f"({len(scored)/len(results):.0%})")

    print(f"\n{'=' * 64}\n  POTENCY CATEGORIES\n{'=' * 64}")
    for cat in sorted(cats):
        ai = {1: 26.5, 2: 100, 3: 400, 4: 1500, 5: 1500}[cat]
        bar = "#" * max(1, round(40 * cats[cat] / max(cats.values())))
        print(f"    Category {cat} ({ai:>6} ng/day)  {cats[cat]:>5}  {bar}")

    strict = [r for r in scored if r["category"] <= 2]
    print(f"\n  {len(strict)} impurities fall in Category 1 or 2 "
          f"(limit 100 ng/day or tighter) - these are the ones to look at first.")
    print(f"\n  written to {OUT_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
