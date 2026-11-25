"""Sanity-check every SMILES in the M7 addendum lookup table."""

import json
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

DATA = Path(__file__).parent.parent / "data" / "m7_addendum_ais.json"


def main() -> int:
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
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
