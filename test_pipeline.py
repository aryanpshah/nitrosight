"""
End-to-end checks on the shared pipeline: resolve -> score -> published override -> ppm.

Run:  python test_pipeline.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from pipeline import classify_input, format_scan_report, scan  # noqa: E402

CLASSIFY = [
    ("CCN(N=O)CC", "smiles"),
    ("62-53-3", "cas"),
    ("N-nitrosomorpholine", "name"),
    ("O=NN1CCOCC1", "smiles"),
]


def run() -> int:
    failures = 0

    print("=" * 72)
    print("  INPUT CLASSIFICATION")
    print("=" * 72)
    for text, want in CLASSIFY:
        got = classify_input(text)
        ok = got == want
        failures += not ok
        print(f"  [{'PASS' if ok else 'FAIL'}] {text[:34]:<34} -> {got} (want {want})")
    print()

    print("=" * 72)
    print("  SCORING FROM SMILES (no network)")
    print("=" * 72)
    for smiles, want_cat in [("CCN(N=O)CC", 1), ("O=NN1CCOCC1", 2),
                             ("O=NN1CCCCC1", 3)]:
        r = scan(smiles, kind="smiles")
        ok = r.cpca is not None and r.cpca.category == want_cat
        failures += not ok
        print(f"  [{'PASS' if ok else 'FAIL'}] {smiles:<18} "
              f"CPCA cat {r.cpca.category if r.cpca else '-'} (want {want_cat})  "
              f"governing {r.final_ai_ng_day} ng/day via {r.final_basis}")
    print()

    print("=" * 72)
    print("  PPM CONTROL LIMIT  (ppm = AI ng/day / MDD mg)")
    print("=" * 72)
    r = scan("CCN(N=O)CC", kind="smiles")
    ppm = r.control_limit_ppm(500)
    want = 26.5 / 500
    ok = ppm is not None and abs(ppm - want) < 1e-9
    failures += not ok
    print(f"  [{'PASS' if ok else 'FAIL'}] NDEA at 500 mg/day -> {ppm:.5g} ppm "
          f"(want {want:.5g})")
    print()

    print("=" * 72)
    print("  PUBLISHED-LIMIT OVERRIDE")
    print("=" * 72)
    # NDMA has a published limit (96 ng/day) and is outside CPCA's Table A.
    r = scan("N-nitrosodimethylamine", kind="name")
    has_pub = r.published is not None
    failures += not has_pub
    print(f"  [{'PASS' if has_pub else 'FAIL'}] NDMA by name -> "
          f"{r.final_ai_ng_day} ng/day via {r.final_basis}")
    if r.cpca:
        print(f"         CPCA alone: applicable={r.cpca.applicable} "
              f"({r.cpca.not_applicable_reason[:60]})")
    print()

    print("=" * 72)
    print("  REPORT RENDERS")
    print("=" * 72)
    text = format_scan_report(scan("CCN(N=O)CC", kind="smiles"), 500)
    ok = "GOVERNING LIMIT" in text and "CPCA BREAKDOWN" in text and "ppm" in text
    failures += not ok
    print(f"  [{'PASS' if ok else 'FAIL'}] report contains the expected sections "
          f"({len(text)} chars)")
    print()

    print("=" * 72)
    print("  WORKBOOK PICTURES: ROW = WHERE THE PICTURE'S CENTRE SITS")
    print("=" * 72)
    # The extractor once used the row holding a picture's top-left corner. Pictures
    # in the vendor workbook start a few points above their own row, so twenty were
    # credited to the previous entry (N-Nitroso Midodrine Impurity got the
    # migalastat picture of the next row). Rows below are 20 pt high.
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "extract_workbook_images",
        Path(__file__).parent / "tools" / "extract_workbook_images.py")
    ewi = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ewi)
    row_tops = [0.0, 0.0, 20.0, 40.0, 60.0, 80.0, 100.0]   # rows 1-5
    PLACEMENT = [
        # (description, picture top, picture bottom, top-left row, expected row)
        ("wholly inside row 3", 42.0, 58.0, 3, 3),
        ("2.5% in row 2, rest below (the midodrine case)", 39.0, 79.0, 2, 3),
        ("mostly in row 2, spilling into row 3", 21.0, 44.0, 2, 2),
    ]
    for label, y0, y1, top_left, want in PLACEMENT:
        got = ewi.centre_row(row_tops, y0, y1)
        ok = got == want
        failures += not ok
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}: row {got}"
              f" (top-left corner says {top_left})")
    print()

    print("=" * 72)
    print(f"  {'ALL PASSED' if failures == 0 else str(failures) + ' FAILURE(S)'}")
    print("=" * 72)

    if failures == 0:
        print()
        print(text)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(run())
