"""
Command-line front end for the CPCA scorer.

    python score.py "CCN(N=O)CC"
    python score.py "CCN(N=O)CC" --name "NDEA"
    python score.py "CCN(N=O)CC" --json
    python score.py "CCN(N=O)CC" --rules reference
    python score.py "CCN(N=O)CC" --compare-rules
"""

import argparse
import dataclasses
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from cpca import PRESET_NAMES, PRESETS, format_report, score_smiles  # noqa: E402
from cpca.policy import DEFAULT  # noqa: E402


def _compare(smiles: str, name: str | None) -> int:
    """Score under every ruleset side by side. Where FDA's guidance and FDA's own
    reference code disagree, this shows whether that disagreement reaches this
    compound."""
    width = max(len(n) for n in PRESET_NAMES)
    rows = []
    for preset in PRESET_NAMES:
        r = score_smiles(smiles, preset)
        if r.applicable:
            site = r.governing_site
            verdict = (f"Category {r.category}, AI = {r.ai_ng_day:g} ng/day"
                       f"   (potency score {site.potency_score})"
                       if site and site.potency_score is not None
                       else f"Category {r.category}, AI = {r.ai_ng_day:g} ng/day")
        else:
            verdict = f"not scored: {r.not_applicable_reason}"
        rows.append((preset, verdict))

    print("=" * 72)
    print(f"  RULESET COMPARISON{'  -  ' + name if name else ''}")
    print("=" * 72)
    print(f"SMILES: {smiles}\n")
    for preset, verdict in rows:
        marker = " *" if preset == DEFAULT.name else "  "
        print(f"{marker} {preset:<{width}}  {verdict}")
    print(f"\n  * default. {PRESETS[DEFAULT.name].description}")

    if len({v for _, v in rows}) == 1:
        print("\n  All rulesets agree on this compound.")
    else:
        print("\n  Rulesets DISAGREE on this compound. See "
              "docs/reference-implementation-divergences.md")
    print("=" * 72)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Score a nitrosamine impurity against FDA's CPCA "
                    "(Carcinogenic Potency Categorization Approach).")
    ap.add_argument("smiles", help="SMILES string for the impurity")
    ap.add_argument("--name", help="compound name, shown in the report header")
    ap.add_argument("--json", action="store_true",
                    help="emit the full result as JSON instead of a text report")
    ap.add_argument("--rules", choices=PRESET_NAMES, default=DEFAULT.name,
                    help="what to follow where FDA's guidance and FDA's reference "
                         "code disagree (default: %(default)s). "
                         "'guidance' follows the published PDF throughout; "
                         "'reference' reproduces FDA's Java tool throughout")
    ap.add_argument("--compare-rules", action="store_true",
                    help="score under every ruleset side by side instead of "
                         "producing one report")
    args = ap.parse_args()

    if args.compare_rules:
        return _compare(args.smiles, args.name)

    result = score_smiles(args.smiles, args.rules)

    if args.json:
        print(json.dumps(dataclasses.asdict(result), indent=2))
    else:
        print(format_report(result, args.name))

    return 0 if result.applicable else 2


if __name__ == "__main__":
    raise SystemExit(main())
