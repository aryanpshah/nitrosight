"""
The whole image path, end to end, on one picture.

    python tools/demo_image_to_score.py data/structure_images/<file>.png
"""

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from cpca import format_report, score_smiles  # noqa: E402
from osr import SETUP_HINT, venv_python  # noqa: E402

WORKER_PY = venv_python()
WORKER = ROOT / "tools" / "osr_worker.py"


def main() -> int:
    if not WORKER_PY.exists():
        print(SETUP_HINT)
        return 1

    image = Path(sys.argv[1]).resolve()
    work = ROOT / "data" / "osr_work"
    work.mkdir(parents=True, exist_ok=True)
    in_file, out_file = work / "demo_in.json", work / "demo_out.json"

    print(f"  STEP 1  picture            {image.name}")

    in_file.write_text(json.dumps([str(image)]), encoding="utf-8")
    # stderr is kept so a failed DECIMER import or model download is diagnosable
    # rather than surfacing as an unexplained missing output file.
    proc = subprocess.run([str(WORKER_PY), str(WORKER),
                           "--in", str(in_file), "--out", str(out_file)],
                          stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)

    if not out_file.exists():
        print("  recognition failed: the worker produced no output")
        if proc.stderr:
            print(proc.stderr.decode("utf-8", "replace").strip())
        return 1

    payload = json.loads(out_file.read_text(encoding="utf-8"))
    smiles = next(iter(payload.values())).get("smiles")
    if not smiles:
        print("  recognition failed")
        return 1

    print(f"  STEP 2  read to SMILES     {smiles}")

    result = score_smiles(smiles)
    if result.applicable:
        print(f"  STEP 3  scored             Category {result.category}, "
              f"{result.ai_ng_day:g} ng/day")
    else:
        print(f"  STEP 3  CPCA does not apply: {result.not_applicable_reason}")
    print()
    print(format_report(result, image.stem))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
