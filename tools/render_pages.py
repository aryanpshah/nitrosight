"""Render specific PDF pages to PNG so figures (which are images, not text) can be read."""

import sys
from pathlib import Path

import pymupdf

HERE = Path(__file__).parent.parent
OUT = HERE / "figures"
OUT.mkdir(exist_ok=True)


def main() -> None:
    pdf_name = sys.argv[1]
    pages = [int(p) for p in sys.argv[2].split(",")]
    zoom = float(sys.argv[3]) if len(sys.argv) > 3 else 2.5

    doc = pymupdf.open(HERE / "guidelines" / pdf_name)
    for p in pages:
        page = doc[p - 1]  # 1-indexed on the command line
        pix = page.get_pixmap(matrix=pymupdf.Matrix(zoom, zoom))
        dest = OUT / f"{Path(pdf_name).stem}_p{p:02d}.png"
        pix.save(dest)
        print(f"  page {p:>2} -> {dest.name}  ({pix.width}x{pix.height})")


if __name__ == "__main__":
    main()
