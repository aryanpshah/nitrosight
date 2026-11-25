"""
When the image and our inference disagree, who is right?

Three independent opinions exist for many of these compounds:

  1. our structure        (CAS lookup, or built from the parent drug)
  2. the image            (DECIMER optical recognition)
  3. the workbook itself  (the 'Nitrosamine impurity' sheet carries a potency category
                           that was filled in by someone else entirely)

Where the workbook agrees with one side, that side wins. This turns a pile of
disagreements into a measured error rate for each method, rather than a shrug.

    python tools/adjudicate_osr.py
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).parent.parent
OSR = ROOT / "data" / "osr_results.json"
OUT = ROOT / "data" / "osr_adjudicated.json"


def main() -> int:
    rows = json.loads(OSR.read_text(encoding="utf-8"))

    # One record per compound; the same compound can appear on two sheets.
    by_name: dict[str, dict] = {}
    for r in rows:
        by_name.setdefault(r["impurity_name"].strip().lower(), r)
    records = list(by_name.values())

    disagreements = [r for r in records if r["verdict"] == "differ"]
    print(f"  distinct compounds with an image : {len(records)}")
    print(f"  image and inference disagree     : {len(disagreements)}")
    print()

    verdict = Counter()
    adjudicated = []

    for r in disagreements:
        sheet = str(r.get("sheet_potency_category") or "").strip()
        ours = r.get("our_category")
        osr = r.get("osr_category")

        if not sheet.isdigit():
            call = "no independent category to judge with"
        else:
            ref = int(sheet)
            if osr == ref and ours != ref:
                call = "image was right, our inference was wrong"
            elif ours == ref and osr != ref:
                call = "our inference was right, the image was misread"
            elif ours == ref and osr == ref:
                call = "both match the workbook (scoring quirk, not a structure issue)"
            else:
                call = "neither matches the workbook"
        verdict[call] += 1
        adjudicated.append(r | {"adjudication": call, "workbook_category": sheet})

    print("=" * 74)
    print("  WHO WAS RIGHT, JUDGED AGAINST THE WORKBOOK'S OWN CATEGORY")
    print("=" * 74)
    for call, n in verdict.most_common():
        print(f"    {n:>4}  {call}")
    print()

    judged = sum(n for call, n in verdict.items()
                 if not call.startswith("no independent"))
    if judged:
        image_right = verdict["image was right, our inference was wrong"]
        ours_right = verdict["our inference was right, the image was misread"]
        print(f"  of the {judged} we could judge:")
        print(f"    image wins     {image_right:>4}  ({image_right/judged:.0%})")
        print(f"    inference wins {ours_right:>4}  ({ours_right/judged:.0%})")
        print()

    # Does the outcome depend on how confident we already were?
    print("=" * 74)
    print("  BY OUR EXISTING CONFIDENCE")
    print("=" * 74)
    for conf in ("high", "medium", "low", None):
        sub = [a for a in adjudicated if a["our_confidence"] == conf
               and not a["adjudication"].startswith("no independent")]
        if not sub:
            continue
        img = sum(1 for a in sub
                  if a["adjudication"].startswith("image was right"))
        our = sum(1 for a in sub
                  if a["adjudication"].startswith("our inference was right"))
        label = conf or "unscored"
        print(f"    {label:<10} image right {img:>3}   inference right {our:>3}   "
              f"(of {len(sub)} judged)")
    print()

    wrong_inference = [a for a in adjudicated
                       if a["adjudication"].startswith("image was right")]
    if wrong_inference:
        print("=" * 74)
        print(f"  {len(wrong_inference)} COMPOUNDS WHERE THE IMAGE CORRECTS US")
        print("=" * 74)
        for a in wrong_inference[:15]:
            print(f"    {a['impurity_name'][:48]:<48} "
                  f"ours Cat {a['our_category']} -> image Cat {a['osr_category']} "
                  f"(workbook: {a['workbook_category']})")

    misread = [a for a in adjudicated
               if a["adjudication"].startswith("our inference was right")]
    if misread:
        print()
        print("=" * 74)
        print(f"  {len(misread)} COMPOUNDS WHERE THE IMAGE WAS MISREAD")
        print("=" * 74)
        for a in misread[:10]:
            print(f"    {a['impurity_name'][:48]:<48} "
                  f"image Cat {a['osr_category']} but workbook says "
                  f"{a['workbook_category']}")

    OUT.write_text(json.dumps(adjudicated, indent=1), encoding="utf-8")
    print(f"\n  wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
