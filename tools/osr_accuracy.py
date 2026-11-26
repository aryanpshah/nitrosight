"""
How accurately does the image path read a structure, and does it matter?

This is the validation table the image-to-category claim rests on. The proposal asks
for exact-structure accuracy, nitroso-site detection sensitivity, feature-level
agreement, category agreement, error types, and the share of cases needing manual
correction. All of that is computed here.

GROUND TRUTH. The reference structure is the one we resolved independently of the
image, by CAS lookup on PubChem. Only rows our resolver labelled "high" confidence are
used for the headline accuracy, because those are CAS-verified and were never derived
from the picture. Medium-confidence rows (structures built from the parent drug) are
reported separately and are NOT ground truth: there the image is often the better
witness, so scoring DECIMER against them would measure the wrong thing.

The distinction that matters for a user is not whether DECIMER drew every atom
correctly. It is whether a chemist who scanned a picture would end up with the right
limit. Both are reported, and they are very different numbers.

    python tools/osr_accuracy.py
    python tools/osr_accuracy.py --errors      # show every mismatch
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import inchi

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from cpca import score_smiles  # noqa: E402
from cpca.rules import NITROSAMINE_CORE  # noqa: E402

RDLogger.DisableLog("rdApp.*")

OSR = ROOT / "data" / "osr_results.json"
OUT = ROOT / "data" / "osr_accuracy.json"


def parse(smiles: str | None):
    if not smiles:
        return None
    return Chem.MolFromSmiles(smiles)


def keys(mol) -> tuple[str, str] | None:
    """(full InChIKey, skeleton block). The skeleton block ignores stereochemistry,
    isotopes and charge, so comparing it isolates connectivity errors."""
    if mol is None:
        return None
    try:
        k = inchi.MolToInchiKey(mol)
    except Exception:  # noqa: BLE001
        return None
    return (k, k.split("-")[0]) if k else None


def nitroso_count(mol) -> int:
    if mol is None:
        return 0
    seen = {m[0] for m in mol.GetSubstructMatches(NITROSAMINE_CORE)}
    return len(seen)


def heavy(mol) -> int:
    return mol.GetNumHeavyAtoms() if mol is not None else 0


def has_isotope(mol) -> bool:
    return mol is not None and any(a.GetIsotope() for a in mol.GetAtoms())


def category_of(smiles: str | None):
    if not smiles:
        return None
    r = score_smiles(smiles)
    return r.category if r.applicable else None


def classify(ref, got, ref_k, got_k) -> str:
    if got is None:
        return "unreadable: DECIMER output could not be parsed as a molecule"
    if got_k is None:
        return "unreadable: structure parsed but no InChI could be generated"
    if ref_k is None:
        return "reference structure could not be keyed"

    if ref_k[0] == got_k[0]:
        return "exact match"
    if ref_k[1] == got_k[1]:
        return "skeleton correct, stereochemistry or isotope label differs"

    ref_n, got_n = nitroso_count(ref), nitroso_count(got)
    if got_n == 0 and ref_n > 0:
        return "N-nitroso group missed entirely"
    if got_n != ref_n:
        return f"wrong number of N-nitroso groups ({got_n} read, {ref_n} expected)"
    if has_isotope(got) and not has_isotope(ref):
        return "spurious isotope label read from the drawing"

    d = heavy(got) - heavy(ref)
    if d > 0:
        return f"skeleton wrong, {d} extra heavy atom(s)"
    if d < 0:
        return f"skeleton wrong, {-d} heavy atom(s) missing"
    return "skeleton wrong, same atom count (bonds or elements misread)"


def difficulty(mol, ref_smiles: str) -> list[str]:
    """The hard-case axes the proposal asks to be reported separately."""
    tags = []
    if mol is None:
        return ["unparsable reference"]
    if nitroso_count(mol) > 1:
        tags.append("multiple N-nitroso sites")
    ri = mol.GetRingInfo()
    if ri.NumRings() >= 3:
        tags.append("three or more rings")
    if any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()):
        tags.append("stereocentres")
    if mol.GetNumHeavyAtoms() >= 30:
        tags.append("large (>=30 heavy atoms)")
    return tags or ["none"]


def pct(n: int, d: int) -> str:
    return f"{n / d:.1%}" if d else "n/a"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--errors", action="store_true", help="list every mismatch")
    args = ap.parse_args()

    rows = json.loads(OSR.read_text(encoding="utf-8"))

    by_name: dict[str, dict] = {}
    for r in rows:
        by_name.setdefault(r["impurity_name"].strip().lower(), r)
    records = list(by_name.values())

    graded, skipped = [], Counter()
    for r in records:
        conf = (r.get("our_confidence") or "").strip()
        ref_smiles = (r.get("our_smiles") or "").strip()
        got_smiles = (r.get("osr_smiles") or "").strip()

        if not ref_smiles:
            skipped["no reference structure to grade against"] += 1
            continue
        if not got_smiles:
            skipped["DECIMER returned nothing"] += 1
            continue

        ref, got = parse(ref_smiles), parse(got_smiles)
        ref_k, got_k = keys(ref), keys(got)
        if ref is None:
            skipped["reference structure unparsable"] += 1
            continue

        graded.append({
            "name": r["impurity_name"],
            "image": r.get("image"),
            "confidence": conf,
            "ref_smiles": ref_smiles,
            "osr_smiles": got_smiles,
            "error": classify(ref, got, ref_k, got_k),
            "exact": bool(ref_k and got_k and ref_k[0] == got_k[0]),
            "skeleton": bool(ref_k and got_k and ref_k[1] == got_k[1]),
            "ref_nitroso": nitroso_count(ref),
            "osr_nitroso": nitroso_count(got),
            "ref_category": category_of(ref_smiles),
            "osr_category": category_of(got_smiles),
            "difficulty": difficulty(ref, ref_smiles),
        })

    high = [g for g in graded if g["confidence"] == "high"]
    other = [g for g in graded if g["confidence"] != "high"]

    print("=" * 76)
    print("  OPTICAL STRUCTURE RECOGNITION: ACCURACY")
    print("=" * 76)
    print(f"    {len(records)} distinct compounds have a catalogue image")
    print(f"    {len(graded)} could be graded against an independently resolved structure")
    for reason, n in skipped.most_common():
        print(f"      {n:>4} skipped: {reason}")
    print()
    print(f"    {len(high)} graded against a CAS-VERIFIED reference (the headline set)")
    print(f"    {len(other)} graded against an inferred reference (reported, not headline)")

    def block(title: str, rows_: list[dict]) -> None:
        if not rows_:
            return
        n = len(rows_)
        exact = sum(1 for g in rows_ if g["exact"])
        skel = sum(1 for g in rows_ if g["skeleton"])
        found = sum(1 for g in rows_ if g["osr_nitroso"] >= 1 and g["ref_nitroso"] >= 1)
        need = sum(1 for g in rows_ if g["ref_nitroso"] >= 1)
        right_n = sum(1 for g in rows_ if g["osr_nitroso"] == g["ref_nitroso"])
        cat_ok = sum(1 for g in rows_ if g["osr_category"] == g["ref_category"])

        print(f"\n{'-' * 76}\n  {title}  (n = {n})\n{'-' * 76}")
        print(f"    exact structure (full InChIKey)        {exact:>4}/{n}  {pct(exact, n)}")
        print(f"    correct skeleton (connectivity only)   {skel:>4}/{n}  {pct(skel, n)}")
        print(f"    N-nitroso group detected at all        {found:>4}/{need}  {pct(found, need)}")
        print(f"    correct NUMBER of N-nitroso sites      {right_n:>4}/{n}  {pct(right_n, n)}")
        print(f"    SAME POTENCY CATEGORY as the reference {cat_ok:>4}/{n}  {pct(cat_ok, n)}")
        print(f"    would need manual correction           "
              f"{n - cat_ok:>4}/{n}  {pct(n - cat_ok, n)}")

    block("Against CAS-verified reference structures", high)
    block("Against inferred reference structures", other)
    block("All graded compounds", graded)

    print(f"\n{'-' * 76}\n  ERROR TYPES (all graded)\n{'-' * 76}")
    for err, n in Counter(g["error"] for g in graded).most_common():
        print(f"    {n:>4}  {err}")

    print(f"\n{'-' * 76}\n  BY DIFFICULTY (category agreement)\n{'-' * 76}")
    axes = Counter(t for g in graded for t in g["difficulty"])
    for tag, total in axes.most_common():
        subset = [g for g in graded if tag in g["difficulty"]]
        ok = sum(1 for g in subset if g["osr_category"] == g["ref_category"])
        print(f"    {tag:<34} {ok:>4}/{len(subset)}  {pct(ok, len(subset))}")

    if args.errors:
        print(f"\n{'-' * 76}\n  EVERY MISMATCH\n{'-' * 76}")
        for g in graded:
            if g["exact"]:
                continue
            nm = g["name"].encode("ascii", "replace").decode()
            print(f"\n  {nm[:66]}  [{g['confidence']}]")
            print(f"    {g['error']}")
            print(f"    category  reference={g['ref_category']}  image={g['osr_category']}")
            print(f"    ref : {g['ref_smiles'][:92]}")
            print(f"    osr : {g['osr_smiles'][:92]}")

    OUT.write_text(json.dumps(graded, indent=1), encoding="utf-8")
    print(f"\n  per-compound results -> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
