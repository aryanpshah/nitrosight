"""
Launch the desktop app.

    python run_gui.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from gui.main import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
