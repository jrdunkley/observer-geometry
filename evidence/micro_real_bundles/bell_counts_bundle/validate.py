from __future__ import annotations

import json
from pathlib import Path

from micro_real_bundles.bell_counts_bundle.run_main import main

ROOT = Path(__file__).resolve().parent


def validate() -> dict[str, object]:
    main()
    summary = json.loads((ROOT / "summary.json").read_text(encoding="utf-8"))
    if summary["raw_downstream_classification"] != "incompatible_by_linear_inconsistency":
        raise AssertionError("raw Bell count bundle should expose finite-sample linear inconsistency")
    if summary["projected_downstream_classification"] != "approximate_common_descent":
        raise AssertionError("projected Bell count bundle should admit audited approximate common descent")
    if summary["projected_best_psd_margin"] <= 0.0:
        raise AssertionError("projected Bell bundle should have positive PSD margin")
    return summary


if __name__ == "__main__":
    print(json.dumps(validate(), indent=2, sort_keys=True))
