from __future__ import annotations

import json
from pathlib import Path

from micro_real_bundles.iris_protocol_mismatch.run_main import main

ROOT = Path(__file__).resolve().parent


def validate() -> dict[str, object]:
    main()
    summary = json.loads((ROOT / "summary.json").read_text(encoding="utf-8"))
    if summary["relation_classification"] != "non_nested_observers":
        raise AssertionError("iris panel relation should be non-nested")
    if summary["completion_classification"] != "underdetermined_affine_family":
        raise AssertionError("iris panel completion should remain honestly underdetermined")
    if summary["refinement_winner"] != "full_iris":
        raise AssertionError("full 4-feature iris panel should be the minimal common refinement")
    return summary


if __name__ == "__main__":
    print(json.dumps(validate(), indent=2, sort_keys=True))
