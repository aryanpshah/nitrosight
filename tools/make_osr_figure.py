"""
Draw the image-recognition accuracy figure from data/osr_accuracy.json.

Every bar is computed from the graded records, so the figure cannot drift away
from the table in the paper. Run it after tools/osr_accuracy.py.

    python tools/make_osr_figure.py                  -> figures/fig-osr-accuracy.png
    python tools/make_osr_figure.py --out some.png
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).parent.parent
ACCURACY = ROOT / "data" / "osr_accuracy.json"

# The headline set is the CAS-verified one: references that were never derived
# from the image being graded. See tools/osr_accuracy.py.
HEADLINE_CONFIDENCE = "high"

DARK = "#22406e"
LIGHT = "#95a4c4"


def metrics(rows: list[dict]) -> dict[str, tuple[float, int, int]]:
    """Each bar as (percentage, numerator, denominator)."""
    def pct(num: int, den: int) -> tuple[float, int, int]:
        return (100.0 * num / den if den else 0.0, num, den)

    # Detection is only meaningful where the reference actually has a nitroso
    # group, so it carries its own denominator. The table says so in a footnote.
    detectable = [r for r in rows if (r["ref_nitroso"] or 0) > 0]

    return {
        "Exact structure": pct(sum(1 for r in rows if r["exact"]), len(rows)),
        "Correct skeleton": pct(sum(1 for r in rows if r["skeleton"]), len(rows)),
        "Correct potency\ncategory": pct(
            sum(1 for r in rows if r["osr_category"] == r["ref_category"]), len(rows)),
        "N-nitroso group\nfound": pct(
            sum(1 for r in detectable if (r["osr_nitroso"] or 0) > 0), len(detectable)),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "figures" / "fig-osr-accuracy.png"))
    args = ap.parse_args()

    rows = json.loads(ACCURACY.read_text(encoding="utf-8"))
    headline = [r for r in rows if r["confidence"] == HEADLINE_CONFIDENCE]

    all_m = metrics(rows)
    cas_m = metrics(headline)
    labels = list(all_m)
    y = range(len(labels))
    height = 0.38

    fig, ax = plt.subplots(figsize=(7.2, 3.9), dpi=300)

    cas_bars = ax.barh([i + height / 2 for i in y],
                       [cas_m[k][0] for k in labels],
                       height=height, color=DARK,
                       label=f"CAS-verified reference (n={len(headline)})")
    all_bars = ax.barh([i - height / 2 for i in y],
                       [all_m[k][0] for k in labels],
                       height=height, color=LIGHT,
                       label=f"All graded (n={len(rows)})")

    for bars, src in ((cas_bars, cas_m), (all_bars, all_m)):
        for bar, key in zip(bars, labels):
            value = src[key][0]
            ax.text(value + 1.2, bar.get_y() + bar.get_height() / 2,
                    f"{value:.1f}%", va="center", ha="left", fontsize=9)

    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel("Percentage of compounds", fontsize=10)
    ax.set_xlim(0, 118)
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.xaxis.grid(True, color="#d9d9d9", linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="y", length=0)

    # Above the axes, so it can never sit on top of a bar or a value label.
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2,
              frameon=False, fontsize=9)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"  wrote {out}")
    for key in labels:
        v, num, den = all_m[key]
        cv, cnum, cden = cas_m[key]
        flat = key.replace("\n", " ")
        print(f"    {flat:24s} CAS {cv:5.1f}% ({cnum}/{cden})   "
              f"all {v:5.1f}% ({num}/{den})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
