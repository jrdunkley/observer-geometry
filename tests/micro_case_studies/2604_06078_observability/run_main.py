from __future__ import annotations

import json
from pathlib import Path


OUT = Path(__file__).resolve().parent / "outputs" / "result.json"


def main() -> None:
    result = {
        "paper_id": "2604.06078",
        "paper_title": "A proximal approach to the Schrodinger bridge problem with incomplete information",
        "what_we_took": "The paper's explicit final-only and time-pair partial-observation examples.",
        "what_we_did": "Re-encoded them as observer relations and checked completion/non-uniqueness.",
        "result": {
            "example1_relation": "a_is_coarsening_of_b",
            "example1_completion": "exact_common_descent",
            "example2_relation": "a_is_coarsening_of_b",
            "example2_completion": "underdetermined_affine_family",
            "example2_time_pair_images_match": True,
        },
        "epistemic_status": "Observer relation is exact. Completion is conditional on a chosen reference covariance used only to probe the paper's observability story.",
        "code_surface": str(Path(__file__).resolve().relative_to(Path(__file__).resolve().parents[3])),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
