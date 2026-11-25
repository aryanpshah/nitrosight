"""
Is our agreement with FDA inflated by having looked at the data?

The proposal asks that development and testing compounds be kept apart. That request
needs care here, because the usual train/test logic does not map cleanly onto a
rule-based model. Nothing in this scorer is a fitted parameter. Every rule was
transcribed from the guidance text and from FDA's published Java source. There is no
coefficient that could have absorbed information from the catalogue.

What COULD have leaked is softer and worth measuring: catalogue disagreements were used
to decide where to re-read the specification, and in four places the guidance and
FDA's own code conflict, so a choice had to be made. If those choices were made to
maximise agreement with the evaluation set, the evaluation is optimistic.

Three checks, reported together:

  1. TEMPORAL HOLD-OUT. FDA's revision table dates each compound's addition to the
     published limits page. Compounds added most recently are scored separately. They
     are the closest thing available to data that could not have shaped the rules.

  2. RULE-CHOICE SENSITIVITY. The four discretionary choices are exposed as rulesets
     (cpca/policy.py). If agreement with FDA's published categories is unchanged when
     those choices are flipped, then the choices demonstrably did not tune that metric,
     and the concern in the paragraph above does not apply to it.

  3. INSPECTED vs NOT INSPECTED. Compounds that appear in our benchmark disagreement
     file were looked at individually during development. Everything else was not.
     Both subsets are scored.

    python tools/holdout.py
"""

from __future__ import annotations

import json
import re
import sys
from collections import Counter
from datetime import date
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from cpca import PRESET_NAMES, score_smiles  # noqa: E402
from resolve import FDALimits  # noqa: E402

RESULTS = ROOT / "data" / "resolved_structures.json"
LIMITS = ROOT / "data" / "fda_published_limits.json"
BENCH = ROOT / "data" / "fda_benchmark_comparison.json"
OUT = ROOT / "data" / "holdout_analysis.json"

# Anything FDA posted on or after this date is treated as the temporal hold-out.
# Chosen as the start of the most recent full year present in the revision table.
CUTOFF = date(2025, 1, 1)


def parse_date(text: str) -> date | None:
    m = re.search(r"(\d{1,2})/(\d{1,2})/(\d{4})", text or "")
    if not m:
        return None
    mm, dd, yyyy = (int(g) for g in m.groups())
    try:
        return date(yyyy, mm, dd)
    except ValueError:
        return None


def norm(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (name or "").lower())


def pct(n: int, d: int) -> str:
    return f"{n / d:.0%}" if d else "n/a"


def main() -> int:
    limits_raw = json.loads(LIMITS.read_text(encoding="utf-8"))
    revisions = limits_raw.get("revisions", [])

    # First date each compound name appears in FDA's revision log.
    first_seen: dict[str, date] = {}
    for r in revisions:
        d = parse_date(r.get("date_posted", ""))
        key = norm(r.get("name", ""))
        if d and key and (key not in first_seen or d < first_seen[key]):
            first_seen[key] = d

    fda = FDALimits()
    rows = json.loads(RESULTS.read_text(encoding="utf-8"))
    seen: dict[str, dict] = {}
    for r in rows:
        if r.get("smiles"):
            seen.setdefault((r.get("catalogue_no") or "").strip()
                            or r["impurity_name"].lower(), r)

    # Every compound where FDA publishes a potency category we can be graded against.
    graded = []
    for r in seen.values():
        pub = fda.lookup(r["impurity_name"])
        if pub is None or pub.category is None:
            continue
        added = first_seen.get(norm(pub.name)) or first_seen.get(norm(r["impurity_name"]))
        graded.append({
            "name": r["impurity_name"],
            "smiles": r["smiles"],
            "published_category": pub.category,
            "date_added": added.isoformat() if added else None,
            "confidence": r.get("confidence"),
        })

    for g in graded:
        for preset in PRESET_NAMES:
            res = score_smiles(g["smiles"], preset)
            g[f"cat_{preset}"] = res.category if res.applicable else None
        g["agrees"] = g["cat_balanced"] == g["published_category"]

    # Which compounds did we inspect individually during development?
    inspected: set[str] = set()
    if BENCH.exists():
        for rec in json.loads(BENCH.read_text(encoding="utf-8")):
            inspected.add(norm(rec.get("name", "")))
    for g in graded:
        g["inspected_during_development"] = norm(g["name"]) in inspected

    n = len(graded)
    agree = sum(1 for g in graded if g["agrees"])

    print("=" * 78)
    print("  HOLD-OUT AND LEAKAGE ANALYSIS")
    print("=" * 78)
    print(f"    {n} compounds have an FDA-published potency category and a structure")
    print(f"    overall agreement: {agree}/{n}  ({pct(agree, n)})")

    # 1. temporal
    dated = [g for g in graded if g["date_added"]]
    undated = n - len(dated)
    recent = [g for g in dated if date.fromisoformat(g["date_added"]) >= CUTOFF]
    older = [g for g in dated if date.fromisoformat(g["date_added"]) < CUTOFF]

    print(f"\n{'-' * 78}\n  1. TEMPORAL HOLD-OUT (FDA revision-table date added)\n{'-' * 78}")
    print(f"    {len(dated)} of {n} could be dated from FDA's revision table "
          f"({undated} undated)")
    if len(dated) < 20:
        print("    NOT INFORMATIVE. FDA's revision table only dates compounds added")
        print("    after the original August 2023 posting, and this catalogue overlaps")
        print("    almost entirely with that original list. Too few dated compounds")
        print("    remain to support a temporal split. Reported for completeness only;")
        print("    check 3 below is the meaningful one.")
    for label, subset in (("added before " + CUTOFF.isoformat(), older),
                          ("added on or after " + CUTOFF.isoformat() + "  <- hold-out",
                           recent)):
        ok = sum(1 for g in subset if g["agrees"])
        print(f"    {label:<48} {ok:>3}/{len(subset):<4} {pct(ok, len(subset))}")

    by_year = Counter(g["date_added"][:4] for g in dated)
    print("\n    agreement by year added:")
    for year in sorted(by_year):
        subset = [g for g in dated if g["date_added"].startswith(year)]
        ok = sum(1 for g in subset if g["agrees"])
        print(f"      {year}   {ok:>3}/{len(subset):<4} {pct(ok, len(subset))}")

    # 2. rule-choice sensitivity
    print(f"\n{'-' * 78}\n  2. RULE-CHOICE SENSITIVITY\n{'-' * 78}")
    print("    If the four discretionary guidance-vs-code choices had been made to fit")
    print("    this evaluation set, flipping them would move the number.\n")
    for preset in PRESET_NAMES:
        ok = sum(1 for g in graded if g[f"cat_{preset}"] == g["published_category"])
        print(f"      ruleset {preset:<10} {ok:>3}/{n}  {pct(ok, n)}")
    spread = {sum(1 for g in graded if g[f"cat_{p}"] == g["published_category"])
              for p in PRESET_NAMES}
    print(f"\n    spread across rulesets: {max(spread) - min(spread)} compound(s).")
    if max(spread) - min(spread) == 0:
        print("    The discretionary choices do not affect this metric at all, so this")
        print("    evaluation cannot have been tuned by making them.")

    # 3. inspected vs not
    print(f"\n{'-' * 78}\n  3. INSPECTED DURING DEVELOPMENT vs NOT\n{'-' * 78}")
    for label, flag in (("looked at individually", True), ("never inspected", False)):
        subset = [g for g in graded if g["inspected_during_development"] is flag]
        ok = sum(1 for g in subset if g["agrees"])
        print(f"    {label:<28} {ok:>3}/{len(subset):<4} {pct(ok, len(subset))}")

    OUT.write_text(json.dumps(graded, indent=1), encoding="utf-8")
    print(f"\n  per-compound results -> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
