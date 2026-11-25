"""
Check that building NDSRIs from parent drugs produces chemically sensible structures.

The propranolol case is the anchor: it has exactly one secondary amine, so the generated
NDSRI must match the N-nitroso propranolol structure we already scored by hand.

Run:  python test_nitrosate.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from rdkit import Chem, RDLogger  # noqa: E402

from cpca import score_smiles  # noqa: E402
from resolve import enumerate_ndsris, groups_from_name  # noqa: E402

RDLogger.DisableLog("rdApp.*")

# parent drug name, parent SMILES, expected number of nitrosatable sites
PARENTS = [
    ("Propranolol", "CC(C)NCC(O)COc1cccc2ccccc12", 1),
    ("Atenolol", "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1", 1),
    ("Abacavir", "NC1=NC2=C(N=CN2[C@H]2C[C@H](CO)C=C2)C(NC2CC2)=N1", 1),
    ("Sertraline", "CNC1CCC(c2ccccc2)c2cc(Cl)c(Cl)cc21", 1),
    ("Varenicline (2 ring N)", "C1N2CC3=CC4=CC=CC=C4C=C3CC2CN1", 1),
]

CANON = Chem.CanonSmiles


def run() -> int:
    failures = 0

    print("=" * 72)
    print("  BUILDING NDSRIs FROM PARENT DRUGS")
    print("=" * 72)

    for name, parent, want_sites in PARENTS:
        cands = enumerate_ndsris(parent, 1)
        ok = len(cands) == want_sites
        if not ok:
            failures += 1
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
        print(f"         sites expected {want_sites}, generated {len(cands)}")
        for c in cands:
            r = score_smiles(c.smiles)
            verdict = (f"Category {r.category}, {r.ai_ng_day:g} ng/day"
                       if r.applicable else f"not applicable: {r.not_applicable_reason[:50]}")
            print(f"         -> {c.smiles}")
            print(f"            {c.site_description}  |  {verdict}")
        print()

    print("=" * 72)
    print("  ANCHOR: generated N-nitroso propranolol matches the hand-written structure")
    print("=" * 72)
    generated = enumerate_ndsris("CC(C)NCC(O)COc1cccc2ccccc12", 1)[0].smiles
    by_hand = "CC(C)N(N=O)CC(O)COc1cccc2ccccc12"
    match = CANON(generated) == CANON(by_hand)
    if not match:
        failures += 1
    print(f"  [{'PASS' if match else 'FAIL'}] generated : {CANON(generated)}")
    print(f"         by hand   : {CANON(by_hand)}")
    print()

    print("=" * 72)
    print("  DINITROSO ENUMERATION")
    print("=" * 72)
    # Piperazine has two secondary amines -> exactly one dinitroso product.
    n = groups_from_name("Some Dinitroso Impurity")
    cands = enumerate_ndsris("C1CNCCN1", n)
    ok = n == 2 and len(cands) == 1
    if not ok:
        failures += 1
    print(f"  [{'PASS' if ok else 'FAIL'}] piperazine, name implies {n} nitroso groups")
    for c in cands:
        r = score_smiles(c.smiles)
        print(f"         -> {c.smiles}  |  Category {r.category}, {r.ai_ng_day:g} ng/day")
    print()

    print("=" * 72)
    print("  TERTIARY DEMETHYLATION ROUTE")
    print("=" * 72)
    # Both of these appear in the catalogue as "N-Nitroso-desmethyl-X", so the
    # demethylation route has to produce them.
    DEMETHYL = [
        ("Sildenafil (N-methylpiperazine)",
         "CCCc1nn(C)c2c1nc(-c1cc(S(=O)(=O)N3CCN(C)CC3)ccc1OCC)[nH]c2=O"),
        ("Trimebutine (dimethylamino)",
         "CCC(c1ccccc1)(N(C)C)COC(=O)c1cc(OC)c(OC)c(OC)c1"),
    ]
    for name, parent in DEMETHYL:
        cands = enumerate_ndsris(parent, 1)
        got = [c for c in cands if c.route == "tertiary_demethylation"]
        ok = len(got) >= 1
        if not ok:
            failures += 1
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: "
              f"{len(got)} demethylation candidate(s)")
        for c in got[:2]:
            r = score_smiles(c.smiles)
            verdict = (f"Category {r.category}, {r.ai_ng_day:g} ng/day"
                       if r.applicable else "not applicable")
            print(f"         -> {c.smiles}")
            print(f"            {verdict}")
        print()

    print("=" * 72)
    print("  MUST PRODUCE NOTHING (no nitrosatable amine at all)")
    print("=" * 72)
    for name, smi in [
        ("Metformin (biguanide)", "CN(C)C(=N)NC(N)=N"),
        ("Paracetamol (amide N)", "CC(=O)Nc1ccc(O)cc1"),
        ("Benzene", "c1ccccc1"),
    ]:
        cands = enumerate_ndsris(smi, 1)
        ok = len(cands) == 0
        if not ok:
            failures += 1
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {len(cands)} candidates")
    print()

    print("=" * 72)
    print(f"  {'ALL PASSED' if failures == 0 else str(failures) + ' FAILURE(S)'}")
    print("=" * 72)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(run())
