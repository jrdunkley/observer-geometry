from __future__ import annotations

import sys
from pathlib import Path


SRC = Path(__file__).resolve().parents[1] / "src"
src_str = str(SRC)
if src_str not in sys.path:
    sys.path.insert(0, src_str)
