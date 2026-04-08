from __future__ import annotations

import json
from pathlib import Path

from external_comparisons.free_gaussian_rg_comparison.run_main import main

ROOT = Path(__file__).resolve().parent


def validate() -> dict[str, object]:
    main()
    summary = json.loads((ROOT / "outputs" / "summary.json").read_text(encoding="utf-8"))
    if summary["schur_residual"] > 1e-12:
        raise AssertionError("QD retained-mode elimination should match Schur complement at machine precision")
    if summary["tower_classification"] != "exact_tower_agreement":
        raise AssertionError("staged elimination should satisfy the tower law")
    if summary["common_refinement_classification"] != "exact_common_refinement":
        raise AssertionError("mode2 and block4 should admit an exact common refinement")
    return summary


if __name__ == "__main__":
    print(json.dumps(validate(), indent=2, sort_keys=True))
