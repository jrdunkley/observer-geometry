from __future__ import annotations

import json
from pathlib import Path


OUT = Path(__file__).resolve().parent / "outputs" / "result.json"


def main() -> None:
    result = {
        "paper_id": "2604.06064",
        "paper_title": "Star-planet magnetic interactions in photoevaporating exoplanets",
        "what_we_took": "The published simulation summary table and the paper's own degeneracy claims.",
        "what_we_did": "Split power-only, power-plus-escape, and power-plus-geometry into separate observers and checked relations/completion.",
        "result": {
            "power_only_vs_escape": "a_is_coarsening_of_b",
            "power_only_vs_geometry": "a_is_coarsening_of_b",
            "escape_vs_geometry": "non_nested_observers",
            "escape_geometry_completion": "underdetermined_affine_family",
            "all_three_completion": "incompatible_by_linear_inconsistency",
        },
        "epistemic_status": "Quotes and table are exact. Observer matrices and descriptive covariances are inferred from the published table axes.",
        "code_surface": str(Path(__file__).resolve().relative_to(Path(__file__).resolve().parents[3])),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
