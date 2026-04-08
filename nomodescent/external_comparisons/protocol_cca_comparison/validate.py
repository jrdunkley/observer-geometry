from __future__ import annotations

import json
from pathlib import Path

from external_comparisons.protocol_cca_comparison.run_main import main

ROOT = Path(__file__).resolve().parent


def validate() -> dict[str, object]:
    main()
    summary = json.loads((ROOT / "outputs" / "summary.json").read_text(encoding="utf-8"))
    if summary["qd_relation_classification"] != "non_nested_observers":
        raise AssertionError("QD should classify the sepal and petal panels as non-nested")
    if summary["qd_completion_classification"] != "underdetermined_affine_family":
        raise AssertionError("QD should preserve honest underdetermination on the two-panel problem")
    if summary["qd_refinement_winner"] != "full_iris":
        raise AssertionError("full_iris should be the minimal common refinement in the finite candidate family")
    return summary


if __name__ == "__main__":
    print(json.dumps(validate(), indent=2, sort_keys=True))
