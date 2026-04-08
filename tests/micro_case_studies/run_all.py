from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CASES = (
    "2604_06078_observability",
    "2604_06108_amp_protocol",
    "2604_06064_spmi_degeneracy",
    "2604_06099_robustness_views",
)


def main() -> None:
    for case in CASES:
        script = ROOT / case / "run_main.py"
        subprocess.run([sys.executable, str(script)], check=True)


if __name__ == "__main__":
    main()
