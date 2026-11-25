"""
Validation suite for the CPCA scorer.

The two anchor cases are FDA's own worked examples (guidance Appendix A, pages 21-23):
if the scorer disagrees with those, it is wrong. The remaining cases are well-known
reference nitrosamines whose categories follow directly from the tables.

Run:  python test_cpca.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from cpca import format_report, score_smiles  # noqa: E402

# name, SMILES, expected category, expected potency score (None if triaged out), source
CASES = [
    (
        "N-nitroso-l-nebivolol",
        "OC(CN(N=O)CC(O)C1CCc2cc(F)ccc2O1)C1CCc2cc(F)ccc2O1",
        4, 4,
        "FDA guidance Example 1, p.21 (score 1+2+1=4)",
    ),
    (
        "N-nitroso-ketamine",
        "CN(N=O)C1(c2ccccc2Cl)CCCCC1=O",
        5, None,
        "FDA guidance Example 2, p.23 (tertiary alpha-carbon)",
    ),
    (
        "N-nitrosodiethylamine (NDEA)",
        "CCN(N=O)CC",
        1, 1,
        "alpha-H 2,2 -> 1; no other features. NDEA underpins the 26.5 ng/day limit.",
    ),
    (
        "N-nitrosomorpholine",
        "O=NN1CCOCC1",
        2, 2,
        "alpha-H 2,2 -> 1; morpholine ring +1",
    ),
    (
        "N-nitrosopiperidine",
        "O=NN1CCCCC1",
        3, 3,
        "alpha-H 2,2 -> 1; 6-membered ring +2",
    ),
    (
        "N-nitrosopyrrolidine",
        "O=NN1CCCC1",
        4, 4,
        "alpha-H 2,2 -> 1; pyrrolidine ring +3",
    ),
    (
        "N-nitrosothiomorpholine",
        "O=NN1CCSCC1",
        4, 4,
        "alpha-H 2,2 -> 1; 6-membered ring containing sulfur +3",
    ),
    (
        "N-nitrosoazepane (7-ring)",
        "O=NN1CCCCCC1",
        2, 2,
        "alpha-H 2,2 -> 1; 7-membered ring +1",
    ),
]

# Structures the CPCA must decline to score rather than guess at.
OUT_OF_SCOPE = [
    ("N-nitrosodimethylamine (NDMA)", "CN(N=O)C",
     "alpha-H pattern 3,3 is absent from Table A"),
    ("N-nitroso-N-methylacetamide", "CC(=O)N(C)N=O",
     "N-nitrosamide: alpha-carbon double bonded to a heteroatom"),
    ("Aspirin (no N-nitroso group)", "CC(=O)Oc1ccccc1C(=O)O",
     "not a nitrosamine"),
]


def run() -> int:
    failures = 0

    print("=" * 72)
    print("  SCORED CASES")
    print("=" * 72)
    for name, smiles, want_cat, want_score, note in CASES:
        r = score_smiles(smiles)
        site = r.governing_site
        got_score = site.potency_score if site else None
        ok = r.category == want_cat and got_score == want_score

        status = "PASS" if ok else "FAIL"
        if not ok:
            failures += 1
        print(f"  [{status}] {name}")
        print(f"         expected: category {want_cat}, score {want_score}")
        print(f"         got     : category {r.category}, score {got_score}, "
              f"AI {r.ai_ng_day} ng/day")
        print(f"         basis   : {note}")
        if not ok and site:
            for f in site.factors:
                print(f"           {f.score:+d}  {f.name}")
            if site.triage_exit:
                print(f"           triage: {site.triage_exit}")
            for n in site.notes:
                print(f"           note: {n}")
        print()

    print("=" * 72)
    print("  MUST DECLINE TO SCORE")
    print("=" * 72)
    for name, smiles, why in OUT_OF_SCOPE:
        r = score_smiles(smiles)
        site = r.governing_site
        declined = (not r.applicable) or site is None or site.potency_score is None
        status = "PASS" if declined else "FAIL"
        if not declined:
            failures += 1
        print(f"  [{status}] {name}")
        print(f"         reason expected: {why}")
        detail = r.not_applicable_reason or "; ".join(
            r.warnings + [n for s in r.sites for n in s.notes])
        print(f"         got            : {detail or 'SCORED (should not have)'}")
        print()

    total = len(CASES) + len(OUT_OF_SCOPE)
    print("=" * 72)
    print(f"  {total - failures}/{total} passed")
    print("=" * 72)

    if failures == 0:
        print()
        print("Sample report:")
        print()
        print(format_report(score_smiles(CASES[0][1]), CASES[0][0]))

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(run())
