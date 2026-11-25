"""
Why are some structures ambiguous, and how much does it actually matter?

When a parent drug has more than one nitrogen that can be nitrosated, the vendor's name
("Nitroso Impurity 1", "Nitroso Impurity 2") does not say which one. We generate every
possibility. Sometimes they all land in the same potency category, in which case the
ambiguity is harmless. Sometimes they do not, and then we genuinely cannot say what the
limit is without knowing the real structure.

This separates those two cases and prints examples.

    python tools/explain_ambiguity.py
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

from rdkit import RDLogger

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from cpca import score_smiles  # noqa: E402

RDLogger.DisableLog("rdApp.*")

RESULTS = ROOT / "data" / "resolved_structures.json"
AI = {1: 26.5, 2: 100.0, 3: 400.0, 4: 1500.0, 5: 1500.0}


def main() -> None:
    rows = json.loads(RESULTS.read_text(encoding="utf-8"))
    seen: dict[str, dict] = {}
    for r in rows:
        seen.setdefault((r.get("catalogue_no") or "").strip()
                        or r["impurity_name"].lower(), r)
    catalogue = list(seen.values())

    ambiguous = [r for r in catalogue if (r.get("n_candidates") or 0) > 1]

    harmless, harmful = [], []
    for r in ambiguous:
        cats = []
        for c in r.get("candidates", []):
            res = score_smiles(c["smiles"])
            if res.applicable:
                cats.append((res.category, c))
        if not cats:
            continue
        distinct = {c for c, _ in cats}
        (harmless if len(distinct) == 1 else harmful).append((r, cats, distinct))

    total = len(catalogue)
    print("=" * 78)
    print("  HOW BIG IS THE AMBIGUITY PROBLEM?")
    print("=" * 78)
    print(f"  catalogue compounds                          : {total}")
    print(f"  built from a parent with >1 possible site    : {len(ambiguous)}")
    print(f"     ... all candidates agree on the category  : {len(harmless)}"
          f"   (ambiguity does not matter)")
    print(f"     ... candidates DISAGREE on the category   : {len(harmful)}"
          f"   (we genuinely cannot say)")
    print()

    if harmful:
        spread = Counter()
        for _, cats, distinct in harmful:
            lo, hi = min(distinct), max(distinct)
            spread[f"Category {lo} vs {hi}"] += 1
        print("  size of the disagreement:")
        for label, n in spread.most_common(8):
            print(f"    {label:<24} {n:>4}")
        print()

    print("=" * 78)
    print("  WORKED EXAMPLES: same vendor name, different possible molecules")
    print("=" * 78)

    # Widest AI spread first: those are the ones worth asking about.
    harmful.sort(key=lambda x: AI[max(x[2])] / AI[min(x[2])], reverse=True)

    for r, cats, distinct in harmful[:6]:
        lo, hi = min(distinct), max(distinct)
        print()
        print(f"  {r['impurity_name']}")
        print(f"     parent drug     : {r['parent_drug']}")
        print(f"     catalogue no.   : {r.get('catalogue_no', '(none)')}")
        print(f"     CAS in workbook : {r.get('cas') or '(none)'}")
        print(f"     possible sites  : {r['n_candidates']}")
        print(f"     the answer could be anywhere from "
              f"{AI[lo]:g} ng/day (Cat {lo}) to {AI[hi]:g} ng/day (Cat {hi})"
              f"  -> a {AI[hi]/AI[lo]:.0f}x difference")
        print()
        for cat, cand in sorted(cats, key=lambda x: x[0])[:4]:
            print(f"       Cat {cat}  {AI[cat]:>6g} ng/day   at {cand['site']}")
            print(f"                              {cand['smiles'][:66]}")
        print(f"     we currently report the worst case: Category {lo}")

    print()
    print("=" * 78)
    print("  WHERE AMBIGUITY IS HARMLESS (for contrast)")
    print("=" * 78)
    for r, cats, distinct in harmless[:4]:
        cat = next(iter(distinct))
        print(f"  {r['impurity_name'][:52]:<52} {r['n_candidates']} sites, "
              f"all Category {cat}")


if __name__ == "__main__":
    main()
