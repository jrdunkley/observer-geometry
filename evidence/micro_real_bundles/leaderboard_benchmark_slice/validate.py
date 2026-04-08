from __future__ import annotations

import json
from pathlib import Path

from micro_real_bundles.leaderboard_benchmark_slice.run_main import main

ROOT = Path(__file__).resolve().parent


def validate() -> dict[str, object]:
    main()
    summary = json.loads((ROOT / "summary.json").read_text(encoding="utf-8"))
    if summary["initial_assembly_classification"] != "underdetermined_evidence":
        raise AssertionError("benchmark bundle should remain unresolved before explicit observer selection")
    if summary["relation_classification"] != "non_nested_observers":
        raise AssertionError("selected benchmark suites should be non-nested observers")
    if summary["completion_classification"] != "underdetermined_affine_family":
        raise AssertionError("benchmark common completion should remain honestly underdetermined")
    return summary


if __name__ == "__main__":
    print(json.dumps(validate(), indent=2, sort_keys=True))
