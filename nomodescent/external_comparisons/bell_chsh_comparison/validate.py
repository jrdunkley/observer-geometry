from __future__ import annotations

import json
from pathlib import Path

from external_comparisons.bell_chsh_comparison.run_main import main

ROOT = Path(__file__).resolve().parent


def validate() -> dict[str, object]:
    main()
    summary = json.loads((ROOT / "outputs" / "summary.json").read_text(encoding="utf-8"))
    if summary["external_method_result"] != "non_obstructed":
        raise AssertionError("arcsine CHSH side should be non-obstructive on the count table")
    if summary["raw_qd_classification"] != "incompatible_by_linear_inconsistency":
        raise AssertionError("raw QD path should detect exact incompatibility")
    if summary["false_collapse_classification"] != "false_collapse_detected":
        raise AssertionError("raw Bell comparison should expose false collapse")
    return summary


if __name__ == "__main__":
    print(json.dumps(validate(), indent=2, sort_keys=True))
