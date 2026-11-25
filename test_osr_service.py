"""
Check the resident recognition service: start once, answer several images, shut down.

The point of keeping it resident is that only the first image pays for loading the
model. This measures that.

Run:  python test_osr_service.py
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from cpca import score_smiles  # noqa: E402
from osr import OsrService, is_available  # noqa: E402

IMAGES = Path(__file__).parent / "data" / "structure_images"

SAMPLES = [
    "unnamed__Bosutinib-Nitroso-Impurity.png",
    "unnamed__Eletriptan-Nitroso-Impurity.png",
    "unnamed__Frovatriptan-Nitroso-Impurity-1.png",
]


def run() -> int:
    if not is_available():
        print("  recognition environment not installed; skipping")
        return 0

    service = OsrService()
    failures = 0

    try:
        print("  starting the service (loads TensorFlow and the model)...")
        started = time.monotonic()
        service.start()
        if not service.wait_until_ready():
            print("  FAIL service never became ready")
            return 1
        print(f"  ready in {time.monotonic() - started:.0f}s\n")

        for name in SAMPLES:
            path = IMAGES / name
            if not path.exists():
                print(f"  [skip] {name} not found")
                continue

            t0 = time.monotonic()
            result = service.recognise(path)
            elapsed = time.monotonic() - t0

            if not result.ok:
                print(f"  [FAIL] {name}: {result.error}")
                failures += 1
                continue

            scored = score_smiles(result.smiles)
            verdict = (f"Category {scored.category}, {scored.ai_ng_day:g} ng/day"
                       if scored.applicable
                       else f"CPCA n/a ({scored.not_applicable_reason[:40]})")
            print(f"  [ok] {name[:46]:<46} {elapsed:>5.1f}s   {verdict}")
            print(f"       {result.smiles[:88]}")
    finally:
        service.stop()
        print("\n  service stopped")

    print(f"\n  {'PASSED' if failures == 0 else str(failures) + ' FAILURE(S)'}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(run())
