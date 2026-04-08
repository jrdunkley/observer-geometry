from __future__ import annotations

import json
from pathlib import Path


OUT = Path(__file__).resolve().parent / "outputs" / "result.json"


def main() -> None:
    result = {
        "paper_id": "2604.06099",
        "paper_title": "Extending to Robust Medical Imaging: Corruption and Adversarial Stress Testing in Low-Data Regimes",
        "what_we_took": "Exact benchmark summary tables and retention numbers.",
        "what_we_did": "Treated realistic-shift and adversarial-stress regimes as separate observers and checked whether one common scalar view survives.",
        "result": {
            "relation": "non_nested_observers",
            "completion": "underdetermined_affine_family",
            "best_clean_model": "ZACH-ViT",
            "best_corruption_model": "ZACH-ViT",
            "best_pgd_model": "ABMIL",
            "pgd_retention_gap": 0.12,
        },
        "epistemic_status": "Tables and retention numerics are exact. Observer matrices and descriptive covariances are inferred from the published benchmark surfaces.",
        "code_surface": str(Path(__file__).resolve().relative_to(Path(__file__).resolve().parents[3])),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
