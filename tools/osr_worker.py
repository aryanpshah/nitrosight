"""
Runs INSIDE .venv-osr (Python 3.13 + TensorFlow). Not importable from the main project.

Two modes, both file-based rather than piped: TensorFlow and DECIMER write progress
output to the console, which would corrupt JSON sent over stdout.

Batch, for processing a list in one go:

    .venv-osr/Scripts/python tools/osr_worker.py --in list.json --out results.json
    ... --preprocess     # strip the pale watermark first

Serve, for interactive use. Loading TensorFlow and the model takes close to a minute, so
a fresh process per image would make the desktop app unusable. In serve mode the model
is loaded once and the process waits for requests:

    .venv-osr/Scripts/python tools/osr_worker.py --serve --dir <work dir>

It watches <dir> for req_<id>.json holding {"image": path}, and answers with
res_<id>.json holding {"smiles": ...} or {"error": ...}. Writing stop.flag ends it.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")   # keep stdout clean for JSON


def strip_watermark(path: Path) -> Path:
    """Catalogue images carry a pale tinted watermark across the structure. The drawing is
    dark or strongly coloured; the watermark is not. Dropping low-contrast pixels leaves
    the structure and removes most of the overlay."""
    from PIL import Image
    import numpy as np

    img = Image.open(path).convert("RGB")
    arr = np.asarray(img).astype("int16")

    lightness = arr.max(axis=2)
    saturation = arr.max(axis=2) - arr.min(axis=2)

    # Keep a pixel if it is dark (ink) or strongly coloured (heteroatom labels).
    keep = (lightness < 170) | (saturation > 90)
    out = np.full(arr.shape, 255, dtype="uint8")
    out[keep] = arr[keep].astype("uint8")

    tmp = Path(tempfile.gettempdir()) / f"osr_clean_{path.stem}.png"
    Image.fromarray(out).save(tmp)
    return tmp


def serve(work_dir: Path, preprocess: bool) -> int:
    """Keep the model loaded and answer one request at a time."""
    import time

    work_dir.mkdir(parents=True, exist_ok=True)
    stop_flag = work_dir / "stop.flag"
    stop_flag.unlink(missing_ok=True)

    from DECIMER import predict_SMILES

    (work_dir / "ready.flag").write_text("ready", encoding="utf-8")
    print("model loaded, waiting for requests", flush=True)

    while not stop_flag.exists():
        requests = sorted(work_dir.glob("req_*.json"))
        if not requests:
            time.sleep(0.15)
            continue

        for req in requests:
            request_id = req.stem[len("req_"):]
            try:
                image = Path(json.loads(req.read_text(encoding="utf-8"))["image"])
                target = strip_watermark(image) if preprocess else image
                payload = {"smiles": predict_SMILES(str(target))}
            except Exception as e:  # noqa: BLE001
                payload = {"error": f"{type(e).__name__}: {e}"}
            (work_dir / f"res_{request_id}.json").write_text(
                json.dumps(payload), encoding="utf-8")
            req.unlink(missing_ok=True)
            print(f"answered {request_id}", flush=True)

    (work_dir / "ready.flag").unlink(missing_ok=True)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="infile",
                    help="JSON file holding a list of image paths")
    ap.add_argument("--out", dest="outfile",
                    help="where to write the JSON results")
    ap.add_argument("--serve", action="store_true",
                    help="stay resident and answer requests one at a time")
    ap.add_argument("--dir", dest="workdir",
                    help="request/response directory, with --serve")
    ap.add_argument("--preprocess", action="store_true",
                    help="remove the pale watermark before recognition")
    args = ap.parse_args()

    if args.serve:
        if not args.workdir:
            ap.error("--serve requires --dir")
        return serve(Path(args.workdir), args.preprocess)
    if not args.infile or not args.outfile:
        ap.error("batch mode requires --in and --out")

    paths = [Path(p) for p in json.loads(Path(args.infile).read_text(encoding="utf-8"))]
    out_path = Path(args.outfile)

    from DECIMER import predict_SMILES  # heavy import, after argument parsing

    results: dict[str, dict] = {}
    for i, path in enumerate(paths, start=1):
        key = str(path)
        try:
            target = strip_watermark(path) if args.preprocess else path
            results[key] = {"smiles": predict_SMILES(str(target))}
        except Exception as e:  # noqa: BLE001
            results[key] = {"error": f"{type(e).__name__}: {e}"}
        print(f"  [{i}/{len(paths)}] {path.name}", flush=True)
        # Written every iteration so a long run is resumable if it is interrupted.
        out_path.write_text(json.dumps(results, indent=1), encoding="utf-8")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
