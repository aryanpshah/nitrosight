"""Extract text from the ICH M7 guideline PDFs so we can read the scoring rules."""

import sys
from pathlib import Path

from pypdf import PdfReader

GUIDELINES = Path(__file__).parent.parent / "guidelines"


def extract(pdf_path: Path) -> Path:
    reader = PdfReader(pdf_path)
    out = pdf_path.with_suffix(".txt")
    chunks = []
    for i, page in enumerate(reader.pages, start=1):
        chunks.append(f"\n\n===== PAGE {i} =====\n")
        chunks.append(page.extract_text() or "")
    out.write_text("".join(chunks), encoding="utf-8")
    print(f"{pdf_path.name}: {len(reader.pages)} pages -> {out.name}")
    return out


if __name__ == "__main__":
    targets = sys.argv[1:] or [p.name for p in sorted(GUIDELINES.glob("*.pdf"))]
    for name in targets:
        extract(GUIDELINES / name)
