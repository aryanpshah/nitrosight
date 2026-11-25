"""Sanity-check every SMILES in the M7 addendum lookup table."""

import json
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).parent.parent
DATA = ROOT / "data" / "m7_addendum_ais.json"

# resolve/fda_limits.py reads the first of these when running from a clone and
# the second when installed. Nothing else keeps them in step, so check it here.
LIMITS = ROOT / "data" / "fda_published_limits.json"
LIMITS_PACKAGED = ROOT / "resolve" / "data" / "fda_published_limits.json"


def check_packaged_limits() -> int:
    """0 if the clone copy and the packaged copy agree, 1 otherwise."""
    if not LIMITS.exists():
        print("  MISSING  data/fda_published_limits.json")
        return 1
    if not LIMITS_PACKAGED.exists():
        print("  MISSING  resolve/data/fda_published_limits.json "
              "(tools/fetch_fda_limits.py writes it)")
        return 1
    if LIMITS.read_bytes() != LIMITS_PACKAGED.read_bytes():
        print("  DRIFT    the data/ and resolve/data/ copies of "
              "fda_published_limits.json differ; re-run tools/fetch_fda_limits.py")
        return 1
    print("  ok  packaged copy of fda_published_limits.json matches data/")
    return 0


def main() -> int:
    bad_limits = check_packaged_limits()
    compounds = json.loads(DATA.read_text(encoding="utf-8"))["compounds"]
    bad = 0
    for c in compounds:
        mol = Chem.MolFromSmiles(c["smiles"])
        if mol is None:
            print(f"  INVALID  {c['name']:<28} {c['smiles']}")
            bad += 1
            continue
        limits = " / ".join(f"{l['ug_day']}" + (f" ({l['route']})" if l["route"] != "all" else "")
                            for l in c["limits"])
        print(f"  ok  {c['name']:<28} {c['smiles']:<26} "
              f"MW={Descriptors.MolWt(mol):7.2f}  AI={limits} ug/day")
    print(f"\n{len(compounds) - bad}/{len(compounds)} SMILES parsed successfully")
    return 1 if (bad or bad_limits) else 0


if __name__ == "__main__":
    raise SystemExit(main())
