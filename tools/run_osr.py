"""
Read the workbook's structure images, and compare what they say against what we inferred.

This is the honest test of the image path. Optical structure recognition on watermarked
vendor images is not going to be perfect, so the result is not "trust the picture" but a
three-way split:

  AGREE      the image and our inferred structure are the same molecule
             -> confidence goes to high, no human needed
  DISAGREE   they are different molecules
             -> flag for a chemist. These are exactly the rows where our name-based
                inference has been getting it wrong.
  NEW        we had no structure at all
             -> the image is the only evidence we have

    python tools/run_osr.py --limit 20            # try a sample first
    python tools/run_osr.py --limit 20 --preprocess
    python tools/run_osr.py                        # all of them
    python tools/run_osr.py --stamp-legacy         # one-time cache migration, see below

CACHING IS BY CONTENT, NOT BY FILE NAME. A cached recognition is reused only when the
image behind the file name still has the same content hash. It used to be keyed on the
path alone, so when extract_workbook_images.py was fixed to put the right picture
behind a file name, the cache would have handed back the old picture's reading. Old
cache entries carry no hash; `--stamp-legacy` records one for each from the image
currently on disk, and must be run BEFORE re-extracting, while those images are still
the ones the entries were made from. An entry with no hash is treated as a miss.

INPUT MUST BE UN-MERGED. "Our" structure is read from resolved_structures.json, so it
must be the output of fetch_structures.py, before merge_osr.py has replaced anything
with an image reading. Otherwise the comparison is between the picture and an earlier
reading of a picture. The tool refuses to run on merged data.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

from cpca import score_smiles  # noqa: E402

RDLogger.DisableLog("rdApp.*")

IMAGES = ROOT / "data" / "structure_images"
INDEX = ROOT / "data" / "structure_images_index.json"
RESULTS = ROOT / "data" / "resolved_structures.json"
OUT = ROOT / "data" / "osr_results.json"
CACHE = ROOT / "data" / "osr_cache.json"
from osr import venv_python  # noqa: E402

WORKER_PY = venv_python()
WORKER = ROOT / "tools" / "osr_worker.py"


def _parent(smiles: str | None) -> Chem.Mol | None:
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    return max(frags, key=lambda m: m.GetNumHeavyAtoms()) if len(frags) > 1 else mol


def canon(smiles: str | None) -> str | None:
    mol = _parent(smiles)
    return Chem.MolToSmiles(mol) if mol else None


def skeleton(smiles: str | None) -> str | None:
    """Connectivity only: no stereochemistry, no isotope labels.

    Recognition routinely reads a drawn 'H' label as deuterium and drops or invents
    wedge bonds. Neither affects a CPCA score, which depends on connectivity and
    hydrogen counts, so a skeleton match is a real match for our purposes.
    """
    mol = _parent(smiles)
    if mol is None:
        return None
    mol = Chem.Mol(mol)
    Chem.RemoveStereochemistry(mol)
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
    return Chem.MolToSmiles(mol)


def file_md5(path: Path) -> str:
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def run_worker(paths: list[Path], preprocess: bool) -> dict[str, dict]:
    """Invoke the recogniser in its own environment, handing work over via files.

    TensorFlow and DECIMER both write to the console, so piping JSON over stdout is not
    safe. Their output is left to stream through to the terminal instead, which also
    makes a long run watchable.
    """
    if not WORKER_PY.exists():
        raise SystemExit(
            f"OSR environment not found at {WORKER_PY}.\n"
            "Create it with:  py -3.13 -m venv .venv-osr\n"
            "                 .venv-osr\\Scripts\\python -m pip install decimer")

    work_dir = ROOT / "data" / "osr_work"
    work_dir.mkdir(parents=True, exist_ok=True)
    in_file = work_dir / "input.json"
    out_file = work_dir / "output.json"
    in_file.write_text(json.dumps([str(p) for p in paths]), encoding="utf-8")
    if out_file.exists():
        out_file.unlink()

    cmd = [str(WORKER_PY), str(WORKER), "--in", str(in_file), "--out", str(out_file)]
    if preprocess:
        cmd.append("--preprocess")

    proc = subprocess.run(cmd)
    if not out_file.exists():
        raise SystemExit(
            f"OSR worker produced no results (exit code {proc.returncode}). "
            "Scroll up for its output.")
    results = json.loads(out_file.read_text(encoding="utf-8"))
    if proc.returncode != 0 and len(results) < len(paths):
        print(f"  ! worker exited {proc.returncode} after {len(results)}/{len(paths)}; "
              "keeping what it produced")
    return results


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, help="only process the first N images")
    ap.add_argument("--preprocess", action="store_true",
                    help="strip the pale watermark before recognition")
    ap.add_argument("--no-cache", action="store_true")
    ap.add_argument("--stamp-legacy", action="store_true",
                    help="record a content hash for cache entries that lack one, from "
                         "the image now on disk, then stop. Run BEFORE re-extracting.")
    args = ap.parse_args()

    if args.stamp_legacy:
        cache = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}
        stamped = missing = 0
        for key, value in cache.items():
            _tag, _, path = key.partition(":")
            if isinstance(value, dict) and "md5" not in value:
                if Path(path).exists():
                    value["md5"] = file_md5(Path(path))
                    stamped += 1
                else:
                    missing += 1
        CACHE.write_text(json.dumps(cache, indent=1), encoding="utf-8")
        print(f"  stamped {stamped} cache entries with the hash of the image now on disk"
              f" ({missing} entries point at a file that no longer exists)")
        return 0

    index = json.loads(INDEX.read_text(encoding="utf-8"))
    rows = json.loads(RESULTS.read_text(encoding="utf-8"))
    if any(r.get("method") == "image_recognition" for r in rows):
        raise SystemExit(
            "resolved_structures.json already contains image-derived rows, so it has been "
            "merged before. Run python tools/fetch_structures.py first. "
            "Merging (or comparing) on already-merged data treats an earlier "
            "image result as our own structure and protects it - which is how a "
            "corrected picture would fail to replace a wrong one.")
    by_name = {r["impurity_name"].strip().lower(): r for r in rows}

    entries = [e for e in index if (IMAGES / e["file"]).exists()]
    if args.limit:
        entries = entries[:args.limit]

    cache: dict[str, dict] = {}
    if CACHE.exists() and not args.no_cache:
        cache = json.loads(CACHE.read_text(encoding="utf-8"))

    tag = "pre" if args.preprocess else "raw"
    # Content-addressed: a reading made from identical bytes is reused whatever the
    # file is called. Two index names can be the SAME file on Windows and macOS
    # ("N-Nitroso-labetalol.png" / "N-Nitroso-Labetalol.png"), and Windows paths
    # compare equal regardless of case, so only one of them is sent to the worker;
    # looking up by name alone left the other with its old reading.
    def by_content() -> dict[str, dict]:
        return {v["md5"]: v for k, v in cache.items()
                if k.startswith(f"{tag}:") and isinstance(v, dict) and v.get("md5")}

    def reading(path: Path) -> dict | None:
        md5 = file_md5(path)
        value = cache.get(f"{tag}:{path}")
        if isinstance(value, dict) and value.get("md5") == md5:
            return value
        return by_content().get(md5)

    def cached(path: Path) -> bool:
        return reading(path) is not None

    todo = list(dict.fromkeys(IMAGES / e["file"] for e in entries
                              if not cached(IMAGES / e["file"])))

    print(f"  {len(entries)} images, {len(todo)} to recognise "
          f"({'watermark stripped' if args.preprocess else 'raw'})")
    if todo:
        print("  running DECIMER, this is slow on CPU...\n")
        fresh = run_worker(todo, args.preprocess)
        for path, value in fresh.items():
            cache[f"{tag}:{path}"] = {**value, "md5": file_md5(Path(path))}
        CACHE.write_text(json.dumps(cache, indent=1), encoding="utf-8")

    verdicts, results = Counter(), []
    for e in entries:
        path = IMAGES / e["file"]
        payload = reading(path) or {}
        if payload:
            cache[f"{tag}:{path}"] = payload      # record it under this name too
        osr_smiles = canon(payload.get("smiles"))

        row = by_name.get(e["impurity_name"].strip().lower())
        ours = canon(row.get("smiles")) if row else None

        osr_cat = our_cat = None
        if osr_smiles:
            scored = score_smiles(osr_smiles)
            osr_cat = scored.category if scored.applicable else None
        if row and row.get("scored"):
            our_cat = row.get("category")

        if osr_smiles is None:
            verdict = "osr_failed"
        elif ours is None:
            verdict = "new"
        elif osr_smiles == ours:
            verdict = "agree_exact"
        elif skeleton(payload.get("smiles")) == skeleton(row.get("smiles")):
            verdict = "agree_skeleton"
        elif osr_cat is not None and osr_cat == our_cat:
            verdict = "agree_category"
        else:
            verdict = "differ"
        verdicts[verdict] += 1

        record = {
            "impurity_name": e["impurity_name"],
            "parent_drug": e.get("parent_drug", ""),
            "sheet_potency_category": e.get("sheet_potency_category", ""),
            "image": e["file"],
            "verdict": verdict,
            "osr_smiles": osr_smiles,
            "our_smiles": ours,
            "our_confidence": row.get("confidence") if row else None,
            "our_category": row.get("category") if row else None,
            "osr_error": payload.get("error"),
            "osr_category": osr_cat,
        }
        results.append(record)

    OUT.write_text(json.dumps(results, indent=1), encoding="utf-8")
    CACHE.write_text(json.dumps(cache, indent=1), encoding="utf-8")

    total = len(results)
    AGREEING = ("agree_exact", "agree_skeleton", "agree_category")
    print(f"\n{'=' * 70}\n  IMAGE vs INFERRED  ({tag})\n{'=' * 70}")
    for key, label in [
            ("agree_exact", "identical structure"),
            ("agree_skeleton", "same structure ignoring stereo / isotopes"),
            ("agree_category", "different structure, same potency category"),
            ("differ", "DIFFERENT potency category (needs a chemist)"),
            ("new", "we had no structure at all"),
            ("osr_failed", "recognition failed or gave invalid SMILES")]:
        n = verdicts[key]
        print(f"    {label:<48} {n:>4}  ({n/max(1,total):.0%})")

    comparable = sum(verdicts[k] for k in AGREEING) + verdicts["differ"]
    if comparable:
        agreeing = sum(verdicts[k] for k in AGREEING)
        print(f"\n    of the {comparable} we could compare, the image and our inference")
        print(f"    reach the same answer {agreeing} times "
              f"({agreeing/comparable:.0%})")

    # Does agreement track the confidence we already assigned?
    print(f"\n{'=' * 70}\n  AGREEMENT BY OUR EXISTING CONFIDENCE\n{'=' * 70}")
    for conf in ("high", "medium", "low"):
        sub = [r for r in results
               if r["our_confidence"] == conf
               and r["verdict"] in AGREEING + ("differ",)]
        if sub:
            ok = sum(1 for r in sub if r["verdict"] in AGREEING)
            print(f"    {conf:<8} {ok:>3}/{len(sub):<3}  ({ok/len(sub):.0%})")

    disagreements = [r for r in results if r["verdict"] == "differ"]
    if disagreements:
        print(f"\n{'=' * 70}\n  DISAGREEMENTS (first 10)\n{'=' * 70}")
        for r in disagreements[:10]:
            print(f"\n  {r['impurity_name'][:58]}   [{r['our_confidence']}]")
            print(f"     from image : Cat {r.get('osr_category')}  {r['osr_smiles'][:62]}")
            print(f"     inferred   : Cat {r.get('our_category')}  {(r['our_smiles'] or '')[:62]}")
            if r["sheet_potency_category"]:
                print(f"     workbook says Category {r['sheet_potency_category']}")

    print(f"\n  wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
