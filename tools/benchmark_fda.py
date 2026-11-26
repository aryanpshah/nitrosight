"""
Head-to-head benchmark against FDA's reference CPCA implementation.

FDA publishes Featurize-Nitrosamines, the Java tool behind the CPCA guidance. It emits
not just a category but every individual feature it scored, so this compares the two
implementations feature by feature rather than only on the final answer. That is the
feature-level agreement a validation paper needs, and it is a far stronger claim than
agreeing on the category alone.

    python tools/benchmark_fda.py --export     # write the input TSV
    python tools/benchmark_fda.py --run        # run the FDA jar
    python tools/benchmark_fda.py --compare    # compare the two
    python tools/benchmark_fda.py              # all three

A molecule carrying several N-nitroso groups produces one row per site from the FDA
tool. We report the most potent site for a molecule, so FDA's rows are collapsed the
same way before anything is compared.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from cpca import PRESET_NAMES, score_smiles  # noqa: E402
from cpca.policy import DEFAULT  # noqa: E402

RESULTS = ROOT / "data" / "resolved_structures.json"
TSV_IN = ROOT / "data" / "fda_benchmark_input.tsv"
TSV_OUT = ROOT / "data" / "fda_benchmark_output.tsv"
REPORT = ROOT / "data" / "fda_benchmark_comparison.json"
JAR_DIR = ROOT / "external" / "featurize-nitrosamines" / "target"

RING_COLUMN = ("NNO in 5- or 6-membered ring (excluding pyrrolidine, 6-membered "
               "S-containing ring and morpholine) score")

# our factor name -> FDA score column
FEATURE_MAP = {
    "Carboxylic acid group anywhere on molecule":
        "Carboxylic acid group anywhere on molecule score",
    "N-nitroso group in a pyrrolidine ring": "NNO in pyrrolidine ring score",
    "N-nitroso group in a 6-membered ring containing sulfur":
        "NNO in 6-membered ring with S score",
    "N-nitroso group in a 5-membered ring": RING_COLUMN,
    "N-nitroso group in a 6-membered ring": RING_COLUMN,
    "N-nitroso group in a morpholine ring": "NNO in morpholine ring score",
    "N-nitroso group in a 7-membered ring": "NNO in a 7-membered ring score",
    "Chains of >=5 consecutive non-hydrogen atoms on both sides of acyclic "
    "N-nitroso group": "Chains of >=5-non-H atoms on both sides of NNO score",
    "Electron-withdrawing group bonded to alpha-carbon on one side":
        "EWG on alpha-carbon on ONLY one side of NNO (excluding carboxylic acid, "
        "aryl and ketone) score",
    "Electron-withdrawing group bonded to alpha-carbon on both sides":
        "EWG on alpha-carbon on BOTH sides of NNO (excluding carboxylic acid, aryl "
        "and ketone) score",
    "Hydroxyl group bonded to beta-carbon on one side":
        "Beta-hydroxyl on ONLY one side score",
    "Hydroxyl group bonded to beta-carbon on both sides":
        "Beta-hydroxyl on BOTH sides score",
    "Aryl group bonded to alpha-carbon (benzylic substituent)":
        "Aryl bonded to alpha-carbon score",
    "Methyl group bonded to beta-carbon": "Methyl group on beta-carbon score",
}


def catalogue() -> list[dict]:
    rows = json.loads(RESULTS.read_text(encoding="utf-8"))
    seen: dict[str, dict] = {}
    for r in rows:
        if r.get("smiles"):
            key = (r.get("catalogue_no") or "").strip() or r["impurity_name"].lower()
            seen.setdefault(key, r)
    return list(seen.values())


def export() -> int:
    rows = catalogue()
    with TSV_IN.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["ID", "SMILES"])
        for i, r in enumerate(rows):
            writer.writerow([f"R{i:04d}", r["smiles"]])
    print(f"  wrote {len(rows)} structures to {TSV_IN.relative_to(ROOT)}")
    return len(rows)


def find_jar() -> Path:
    jars = (sorted(JAR_DIR.glob("*jar-with-dependencies.jar"))
            or sorted(JAR_DIR.glob("*.jar")))
    if not jars:
        raise SystemExit(
            f"No jar in {JAR_DIR}. Build it first:\n"
            "  cd external/featurize-nitrosamines\n"
            "  ../apache-maven-3.9.16/bin/mvn clean package -DskipTests")
    return jars[-1]


def run_fda() -> None:
    jar = find_jar()
    print(f"  running {jar.name}")
    proc = subprocess.run(
        ["java", "-jar", str(jar), "-i", str(TSV_IN), "-o", str(TSV_OUT)],
        capture_output=True, text=True, encoding="utf-8", errors="replace")
    if proc.returncode != 0 or not TSV_OUT.exists():
        print((proc.stdout or "")[-2000:])
        print((proc.stderr or "")[-2000:])
        raise SystemExit(f"FDA tool exited {proc.returncode}")
    print(f"  wrote {TSV_OUT.relative_to(ROOT)}")


def load_fda() -> dict[str, dict]:
    """Most potent site per structure, which is how we report a molecule."""
    best: dict[str, dict] = {}
    with TSV_OUT.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            name = (row.get("Structure_Name") or "").strip()
            if not name:
                continue
            try:
                cat = int(row["Potency Category"])
            except (KeyError, ValueError, TypeError):
                continue
            if name not in best or cat < int(best[name]["Potency Category"]):
                best[name] = row
    return best


def our_features(site) -> dict[str, int]:
    """Our factor scores collapsed onto FDA's column names."""
    out: dict[str, int] = {}
    for factor in site.factors:
        column = FEATURE_MAP.get(factor.name)
        if column:
            out[column] = out.get(column, 0) + factor.score
    return out


def compare(ruleset: str = DEFAULT.name, quiet: bool = False) -> tuple[int, int]:
    rows = catalogue()
    fda = load_fda()
    if not quiet:
        print(f"  FDA scored {len(fda)} of {len(rows)} structures")
        print(f"  our ruleset: {ruleset}\n")

    verdicts = Counter()
    feature_mismatch = Counter()
    records = []

    for i, r in enumerate(rows):
        key = f"R{i:04d}"
        theirs = fda.get(key)
        mine = score_smiles(r["smiles"], ruleset)

        if theirs is None:
            verdicts["FDA produced no row" if mine.applicable else "both declined"] += 1
            continue
        if not mine.applicable:
            verdicts["we declined, FDA scored"] += 1
            records.append({
                "id": key, "name": r["impurity_name"], "smiles": r["smiles"],
                "verdict": "we declined", "our_category": None,
                "fda_category": int(theirs["Potency Category"]),
                "our_reason": mine.not_applicable_reason})
            continue

        site = mine.governing_site
        their_cat = int(theirs["Potency Category"])
        agree = their_cat == mine.category
        verdicts["category agrees" if agree else "category differs"] += 1

        ours = our_features(site)
        diffs = []
        # Sorted, not just deduplicated. Python randomises string hashing per
        # process, so iterating the bare set writes feature_diffs in a different
        # order on every run. That makes the committed JSON churn between runs
        # and makes a diff look like the results moved when nothing did.
        for column in sorted(set(FEATURE_MAP.values())):
            want = ours.get(column, 0)
            try:
                got = int(theirs.get(column, 0) or 0)
            except (TypeError, ValueError):
                continue
            if want != got:
                diffs.append({"feature": column, "fda": got, "ours": want})
                feature_mismatch[column] += 1

        their_ah = (theirs.get("Alpha-Hydrogens") or "").strip()
        our_ah = (f"{site.alpha_h_counts[0]},{site.alpha_h_counts[1]}"
                  if site.alpha_h_counts else "")
        if their_ah and our_ah and their_ah != our_ah:
            feature_mismatch["Alpha-Hydrogens"] += 1

        if not agree or diffs:
            records.append({
                "id": key, "name": r["impurity_name"], "smiles": r["smiles"],
                "verdict": "category differs" if not agree else "features differ",
                "fda_category": their_cat, "our_category": mine.category,
                "fda_score": theirs.get("Potency Score"),
                "our_score": site.potency_score,
                "fda_alpha_h": their_ah, "our_alpha_h": our_ah,
                "our_triaged": bool(site.triage_exit),
                "our_triage_reason": site.triage_exit,
                "feature_diffs": diffs,
                "our_confidence": r.get("confidence")})

    agree = verdicts["category agrees"]
    total = agree + verdicts["category differs"]

    if quiet:
        return agree, total

    print("=" * 74)
    print("  CATEGORY AGREEMENT WITH FDA's REFERENCE IMPLEMENTATION")
    print("=" * 74)
    for label, n in verdicts.most_common():
        print(f"    {n:>5}  {label}")
    if total:
        print(f"\n    {agree}/{total} ({agree / total:.1%}) agree on potency category")

    if feature_mismatch:
        print(f"\n{'=' * 74}\n  FEATURE-LEVEL DISAGREEMENTS\n{'=' * 74}")
        for column, n in feature_mismatch.most_common():
            print(f"    {n:>5}  {column[:64]}")
    else:
        print("\n  no feature-level disagreements")

    REPORT.write_text(json.dumps(records, indent=1), encoding="utf-8")
    print(f"\n  {len(records)} differing structures -> {REPORT.relative_to(ROOT)}")
    return agree, total


def compare_all_rulesets() -> int:
    """Agreement under each ruleset. The 'reference' preset should score highest by
    construction, since it is the one that sets out to reproduce FDA's tool. Anything
    else would mean we had not understood the divergences."""
    print("=" * 74)
    print("  AGREEMENT WITH FDA's REFERENCE TOOL, BY RULESET")
    print("=" * 74)
    for name in PRESET_NAMES:
        agree, total = compare(name, quiet=True)
        marker = " *" if name == DEFAULT.name else "  "
        pct = f"{agree / total:.1%}" if total else "n/a"
        print(f"{marker} {name:<10}  {agree:>3}/{total}  ({pct})")
    print("\n  * default ruleset")
    print("  'reference' is expected to score highest: it is the preset that")
    print("  deliberately reproduces FDA's tool, divergences from the published")
    print("  guidance included. See docs/reference-implementation-divergences.md")
    print("=" * 74)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--export", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--compare", action="store_true")
    ap.add_argument("--rules", choices=PRESET_NAMES, default=DEFAULT.name,
                    help="which of our rulesets to compare against FDA's tool "
                         "(default: %(default)s)")
    ap.add_argument("--all-rulesets", action="store_true",
                    help="report agreement under every ruleset")
    args = ap.parse_args()

    if args.all_rulesets:
        return compare_all_rulesets()

    if not any([args.export, args.run, args.compare]):
        args.export = args.run = args.compare = True

    if args.export:
        export()
    if args.run:
        run_fda()
    if args.compare:
        compare(args.rules)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
