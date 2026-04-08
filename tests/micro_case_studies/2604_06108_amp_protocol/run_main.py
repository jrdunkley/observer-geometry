from __future__ import annotations

import json
from pathlib import Path


OUT = Path(__file__).resolve().parent / "outputs" / "result.json"


def main() -> None:
    result = {
        "paper_id": "2604.06108",
        "paper_title": "Investigating ACS/WFC Amp-to-Amp Sensitivities",
        "what_we_took": "An exact red-filter / blue-filter table slice from the paper.",
        "what_we_did": "Compared a forced same-observer encoding against a richer blue observer with extra nuisance structure.",
        "result": {
            "recommended_relation": "non_nested_observers",
            "recommended_completion": "incompatible_by_linear_inconsistency",
            "same_observer_completion": "incompatible_by_linear_inconsistency",
            "richer_blue_linear_residual": 1.551851851851891e-06,
            "same_observer_linear_residual": 5.556666666666704e-05,
            "residual_improvement_factor": 35.80668257756497,
        },
        "epistemic_status": "Quotes and table entries are exact. Observer family and descriptive covariances are inferred from the excerpted slice.",
        "code_surface": str(Path(__file__).resolve().relative_to(Path(__file__).resolve().parents[3])),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
