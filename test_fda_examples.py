"""
Check our CPCA implementation against FDA's own reference tool.

FDA publishes Featurize-Nitrosamines (github.com/FDA/featurize-nitrosamines), the
reference implementation behind the CPCA guidance. Its README documents four worked
examples with the complete feature-level breakdown, so they can be checked without
building the Java project.

This is a stronger test than category agreement alone: it compares every individual
feature we score against FDA's own, which is exactly the "feature-level agreement"
metric a validation paper needs.

Run:  python test_fda_examples.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from cpca import score_smiles  # noqa: E402

# From the FDA README's documented output, verbatim.
# (smiles, alpha_hydrogens, alpha_score, potency_score, category, {feature: score})
CASES = [
    ("C=CCN(CC=C)N=O", "2,2", 1, 1, 1, {}),
    ("CN(N=O)C1=CC=NC=C1", "0,3", 2, 2, 2, {}),
    ("CN(N=O)C1=CC=CN=C1", "0,3", 2, 2, 2, {}),
    ("OC(=O)C1CCCCN1N=O", "1,2", 3, 8, 4, {
        "Carboxylic acid group anywhere on molecule": 3,
        "N-nitroso group in a 6-membered ring": 2,
    }),
]


def ours(smiles: str) -> dict:
    result = score_smiles(smiles)
    if not result.applicable:
        return {"error": result.not_applicable_reason}
    site = result.governing_site
    a, b = site.alpha_h_counts
    return {
        "alpha_hydrogens": f"{a},{b}",
        "alpha_score": site.alpha_h_score,
        "potency_score": site.potency_score,
        "category": result.category,
        "features": {f.name: f.score for f in site.factors
                     if f.name != "alpha-Hydrogen Score"},
    }


def run() -> int:
    failures = 0
    print("=" * 76)
    print("  OUR IMPLEMENTATION vs FDA's REFERENCE TOOL (documented examples)")
    print("=" * 76)
    print()

    for smiles, a_h, a_score, score, cat, features in CASES:
        got = ours(smiles)
        if "error" in got:
            print(f"  [FAIL] {smiles}\n         {got['error']}\n")
            failures += 1
            continue

        checks = [
            ("alpha-hydrogens", a_h, got["alpha_hydrogens"]),
            ("alpha-hydrogen score", a_score, got["alpha_score"]),
            ("potency score", score, got["potency_score"]),
            ("potency category", cat, got["category"]),
        ]
        for name, want in features.items():
            checks.append((name, want, got["features"].get(name)))

        bad = [c for c in checks if c[1] != c[2]]
        failures += bool(bad)
        print(f"  [{'PASS' if not bad else 'FAIL'}] {smiles}")
        for name, want, mine in checks:
            mark = " " if want == mine else "  <-- MISMATCH"
            print(f"         {name:<36} FDA={str(want):<8} ours={str(mine):<8}{mark}")
        extra = set(got["features"]) - set(features)
        if extra:
            print(f"         extra features we scored: "
                  + ", ".join(f"{k} ({got['features'][k]:+d})" for k in sorted(extra)))
        print()

    print("=" * 76)
    print(f"  {len(CASES) - failures}/{len(CASES)} FDA reference examples reproduced")
    print("=" * 76)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(run())
