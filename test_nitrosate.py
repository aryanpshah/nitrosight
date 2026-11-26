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
from resolve import (  # noqa: E402
    enumerate_ndsris, groups_from_name, implied_route, name_implied_index,
)

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

    print("  VENDOR NAME PICKS THE FORMATION ROUTE")
    print("=" * 72)
    # A tertiary amine cannot be nitrosated directly, so an impurity named
    # desmethyl-X must have lost a methyl first. The absence of the word is the
    # ordinary reading of an NDSRI name: direct nitrosation of the parent's own
    # secondary amine.
    ROUTE_NAMES = [
        ("N-Nitroso Desmethyl Tetracaine", "tertiary_demethylation"),
        ("N-Nitroso-N-Desmethyl Zolmitriptan", "tertiary_demethylation"),
        ("N-nitroso de-methylated ranitidine", "tertiary_demethylation"),
        ("N-Nitroso Tetracaine", "secondary_amine"),
        ("Abacavir Nitroso 1", "secondary_amine"),
        ("", None),
        ("   ", None),
    ]
    for name, want in ROUTE_NAMES:
        got = implied_route(name)
        ok = got == want
        if not ok:
            failures += 1
        print(f"  [{'PASS' if ok else 'FAIL'}] {name!r:<38} -> {got}")
    print()

    print("=" * 72)
    print("  NAME SINGLES OUT ONE CANDIDATE, OR NONE")
    print("=" * 72)
    # Exactly one matching route -> its index. Two matching routes (tigecycline
    # has two demethylatable nitrogens) -> None: the name does not say which.
    INDEX_CASES = [
        ("N-Nitroso Tetracaine", ["secondary_amine", "tertiary_demethylation"], 0),
        ("N-Nitroso Desmethyl Tetracaine", ["secondary_amine", "tertiary_demethylation"], 1),
        ("N-Nitroso Desmethyl Tigecycline 1",
         ["secondary_amine", "tertiary_demethylation", "tertiary_demethylation"], None),
        ("N-Nitroso Tigecycline",
         ["secondary_amine", "tertiary_demethylation", "tertiary_demethylation"], 0),
        ("N-Nitroso Something", ["tertiary_demethylation"], None),
        ("", ["secondary_amine", "tertiary_demethylation"], None),
    ]
    for name, routes, want in INDEX_CASES:
        got = name_implied_index(name, routes)
        ok = got == want
        if not ok:
            failures += 1
        print(f"  [{'PASS' if ok else 'FAIL'}] {name!r:<36} -> {got}")
    print()

    print("=" * 72)
    print("  TETRACAINE: THE NAME SEPARATES TWO CANDIDATES")
    print("=" * 72)
    # Tetracaine has a secondary butylamino on the ring AND a tertiary
    # dimethylaminoethyl ester tail, so it yields one candidate by each route.
    # Before the name was consulted both vendor names got the worst-scoring
    # candidate, which made them indistinguishable downstream even though FDA
    # publishes different categories for them.
    tetracaine = "CCCCNc1ccc(C(=O)OCCN(C)C)cc1"
    cands = enumerate_ndsris(tetracaine, 1)
    two = len(cands) == 2
    if not two:
        failures += 1
    print(f"  [{'PASS' if two else 'FAIL'}] tetracaine yields 2 candidates, got {len(cands)}")
    picked = {}
    for name in ("N-Nitroso Tetracaine", "N-Nitroso Desmethyl Tetracaine"):
        matching = [c for c in cands if c.route == implied_route(name)]
        ok = len(matching) == 1
        if not ok:
            failures += 1
        print(f"  [{'PASS' if ok else 'FAIL'}] {name!r} matches {len(matching)} candidate(s)")
        if ok:
            picked[name] = matching[0]
            print(f"         -> {matching[0].smiles}")
            print(f"            {matching[0].site_description}")
    distinct = len({c.smiles for c in picked.values()}) == len(picked) == 2
    if not distinct:
        failures += 1
    print(f"  [{'PASS' if distinct else 'FAIL'}] the two names get DIFFERENT structures")
    print()

    print("=" * 72)
    print(f"  {'ALL PASSED' if failures == 0 else str(failures) + ' FAILURE(S)'}")
    print("=" * 72)
    print("=" * 72)

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(run())

