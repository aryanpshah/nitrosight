"""
Three-way benchmark: this scorer, FDA's reference implementation, and the Novartis tool.

The proposal asks that our CPCA scores be cross-checked against other published
implementations rather than against FDA's alone. This runs the comparison that can
actually be automated.

  FDA featurize-nitrosamines   Java, batch TSV. Run via tools/benchmark_fda.py.
  Novartis NDSRIs_in_silico    Python. Its README documents combine_all_rules_together
                               in CPCA_rules.py as the batch entry point, so it is
                               called directly here.
  NitroLimit                   PyQt5 desktop GUI, one SMILES at a time, no batch or
                               command-line entry point, so it cannot be included in an
                               automated comparison. Recorded as such.
  USP Excel calculator         A spreadsheet in which the structure box is for reference
                               only and scoring is done by hand. Not automatable.

The Novartis tool is vendored under external/NDSRIs_in_silico_tool (Apache-2.0). Its
CPCA_rules module imports only RDKit and base64 and performs no I/O, which was checked
before running it.

    python tools/benchmark_tools.py
    python tools/benchmark_tools.py --errors
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

NOVARTIS = ROOT / "external" / "NDSRIs_in_silico_tool"

from rdkit import RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")

from cpca import score_smiles  # noqa: E402
from resolve import FDALimits  # noqa: E402

RESULTS = ROOT / "data" / "resolved_structures.json"
FDA_OUT = ROOT / "data" / "fda_benchmark_output.tsv"
OUT = ROOT / "data" / "tool_comparison.json"


def load_novartis():
    """Import Novartis CPCA_rules, working around a defect in the published file.

    As published, CPCA_rules.py opens with a thirteen-line licence header that has no
    comment markers, so the file is not valid Python and cannot be imported at all.
    template.py has the same defect; streamlit.py in the same repository comments its
    header correctly, so this looks like an oversight rather than a convention.

    Rather than edit their source, the header is commented out in memory before the
    module is compiled. Their file on disk is left exactly as cloned, and the patch is
    reported so it appears in the record of any benchmark run.
    """
    path = NOVARTIS / "CPCA_rules.py"
    if not path.exists():
        return None

    source = path.read_text(encoding="utf-8")
    lines = source.splitlines()

    # First line that could legally begin a Python module.
    start = 0
    for i, line in enumerate(lines):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if s.startswith(("import ", "from ", '"""', "'''")):
            start = i
            break
    else:
        start = 0

    patched = 0
    if start:
        for i in range(start):
            if lines[i].strip() and not lines[i].lstrip().startswith("#"):
                lines[i] = "# " + lines[i]
                patched += 1

    if patched:
        print(f"  note: commented out {patched} line(s) of uncommented licence header "
              f"in CPCA_rules.py so it could be imported.")
        print("        Their file on disk is unchanged. This is a defect in the "
              "published repository,")
        print("        not a modification to the rules being benchmarked.")

    import types  # noqa: PLC0415
    module = types.ModuleType("novartis_cpca_rules")
    module.__file__ = str(path)
    sys.path.insert(0, str(NOVARTIS))
    try:
        exec(compile("\n".join(lines), str(path), "exec"), module.__dict__)  # noqa: S102
    except Exception as exc:  # noqa: BLE001
        print(f"  ! Novartis tool could not be loaded: {type(exc).__name__}: {exc}")
        return None
    if not hasattr(module, "combine_all_rules_together"):
        print("  ! Novartis module loaded but combine_all_rules_together is missing")
        return None
    return module


def novartis_category(mod, smiles: str):
    """(category, note). None category means the tool declined or failed.

    CPCA_rules.py carries leftover debug print statements that emit atom indices and
    element symbols for every molecule. They are swallowed here so the comparison
    output stays readable; nothing else about the call is changed.
    """
    import contextlib  # noqa: PLC0415
    import io  # noqa: PLC0415
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            message, score = mod.combine_all_rules_together(smiles)
    except Exception as exc:  # noqa: BLE001
        return None, f"crashed: {type(exc).__name__}"
    if isinstance(message, list) and message and "No nitrosamine pattern" in str(message[0]):
        return None, "no nitrosamine pattern found"
    try:
        return mod.score_to_category(score), ""
    except Exception as exc:  # noqa: BLE001
        return None, f"category mapping failed: {type(exc).__name__}"


def load_fda() -> dict[str, int]:
    """Most potent site per structure, keyed by the R#### id used on export."""
    import csv  # noqa: PLC0415
    best: dict[str, int] = {}
    if not FDA_OUT.exists():
        return best
    with FDA_OUT.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            name = (row.get("Structure_Name") or "").strip()
            try:
                cat = int(row["Potency Category"])
            except (KeyError, ValueError, TypeError):
                continue
            if name and (name not in best or cat < best[name]):
                best[name] = cat
    return best


def catalogue() -> list[dict]:
    rows = json.loads(RESULTS.read_text(encoding="utf-8"))
    seen: dict[str, dict] = {}
    for r in rows:
        if r.get("smiles"):
            seen.setdefault((r.get("catalogue_no") or "").strip()
                            or r["impurity_name"].lower(), r)
    return list(seen.values())


def pct(n: int, d: int) -> str:
    return f"{n / d:.1%}" if d else "n/a"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--errors", action="store_true")
    args = ap.parse_args()

    mod = load_novartis()
    if mod is None:
        print("  Novartis tool unavailable. Clone it into external/ first:")
        print("    git clone https://github.com/Novartis/NDSRIs_in_silico_tool")
        return 1

    rows = catalogue()
    fda = load_fda()
    fdalimits = FDALimits()

    records = []
    for i, r in enumerate(rows):
        smiles = r["smiles"]
        ours = score_smiles(smiles)
        nov_cat, nov_note = novartis_category(mod, smiles)
        pub = fdalimits.lookup(r["impurity_name"])
        records.append({
            "id": f"R{i:04d}",
            "name": r["impurity_name"],
            "smiles": smiles,
            "ours": ours.category if ours.applicable else None,
            "ours_reason": "" if ours.applicable else ours.not_applicable_reason,
            "fda": fda.get(f"R{i:04d}"),
            "novartis": nov_cat,
            "novartis_note": nov_note,
            "published": pub.category if pub else None,
        })

    print("=" * 78)
    print("  THREE-WAY TOOL COMPARISON")
    print("=" * 78)
    print(f"    {len(records)} catalogue structures\n")

    for label, key in (("this scorer", "ours"), ("FDA featurize-nitrosamines", "fda"),
                       ("Novartis NDSRIs_in_silico", "novartis")):
        scored = sum(1 for r in records if r[key] is not None)
        print(f"    {label:<30} produced a category for {scored}/{len(records)}")

    def agreement(a: str, b: str) -> tuple[int, int]:
        both = [r for r in records if r[a] is not None and r[b] is not None]
        return sum(1 for r in both if r[a] == r[b]), len(both)

    print(f"\n{'-' * 78}\n  PAIRWISE AGREEMENT (where both produced a category)\n{'-' * 78}")
    for a, b, label in (("ours", "fda", "this scorer vs FDA"),
                        ("ours", "novartis", "this scorer vs Novartis"),
                        ("fda", "novartis", "FDA vs Novartis")):
        ok, n = agreement(a, b)
        print(f"    {label:<32} {ok:>4}/{n:<5} {pct(ok, n)}")

    three = [r for r in records
             if r["ours"] is not None and r["fda"] is not None
             and r["novartis"] is not None]
    unanimous = sum(1 for r in three
                    if r["ours"] == r["fda"] == r["novartis"])
    print(f"\n    all three agree                  {unanimous:>4}/{len(three):<5} "
          f"{pct(unanimous, len(three))}")

    print(f"\n{'-' * 78}\n  AGAINST FDA's PUBLISHED CATEGORIES\n{'-' * 78}")
    for key, label in (("ours", "this scorer"), ("fda", "FDA tool"),
                       ("novartis", "Novartis tool")):
        both = [r for r in records if r["published"] is not None and r[key] is not None]
        ok = sum(1 for r in both if r[key] == r["published"])
        print(f"    {label:<32} {ok:>4}/{len(both):<5} {pct(ok, len(both))}")

    print(f"\n{'-' * 78}\n  WHY THE NOVARTIS TOOL DECLINED\n{'-' * 78}")
    for note, n in Counter(r["novartis_note"] for r in records
                           if r["novartis"] is None).most_common(6):
        print(f"    {n:>4}  {note or '(no note)'}")

    if args.errors:
        print(f"\n{'-' * 78}\n  THREE-WAY DISAGREEMENTS\n{'-' * 78}")
        shown = 0
        for r in three:
            if r["ours"] == r["fda"] == r["novartis"]:
                continue
            shown += 1
            if shown > 25:
                break
            nm = r["name"].encode("ascii", "replace").decode()
            print(f"\n  {nm[:62]}")
            print(f"    ours={r['ours']}  fda={r['fda']}  novartis={r['novartis']}"
                  f"  published={r['published']}")
            print(f"    {r['smiles'][:92]}")

    OUT.write_text(json.dumps(records, indent=1), encoding="utf-8")
    print(f"\n  per-compound results -> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
