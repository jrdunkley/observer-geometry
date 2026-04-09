from __future__ import annotations

from pathlib import Path

import numpy as np

import nomogeo
from nomogeo import connection_current, fixed_observer_coordinates, observer_transition

ROOT = Path(__file__).resolve().parents[1]


def run_install_surface_smoke() -> dict[str, object]:
    package_file = Path(nomogeo.__file__).resolve()
    expected_root = (ROOT / "src").resolve()
    imported_from_workspace = package_file.is_relative_to(expected_root)

    H = np.diag([2.0, 3.0, 5.0, 7.0])
    C = np.array([[1.0, 0.0, 1.0, 0.0], [0.0, 1.0, 0.0, 1.0]])
    Delta = np.diag([0.5, -0.25, 0.125, 0.0])

    chart = fixed_observer_coordinates(H, C)
    current = connection_current(H, C, Delta)
    transition = observer_transition(H, C, C)

    return {
        "package_file": str(package_file),
        "imported_from_workspace": imported_from_workspace,
        "phi_trace": float(np.trace(chart.phi)),
        "current_norm": float(np.linalg.norm(current.current, ord="fro")),
        "forcing_residual": float(np.linalg.norm(current.forcing - current.q, ord="fro")),
        "identity_transition_residual": float(transition.residual),
    }


def main() -> None:
    result = run_install_surface_smoke()
    print(result)
    if not result["imported_from_workspace"]:
        raise SystemExit("nomogeo did not import from this workspace")
    if result["forcing_residual"] > 1e-10:
        raise SystemExit("connection current forcing identity failed in install smoke")
    if result["identity_transition_residual"] > 1e-10:
        raise SystemExit("observer transition identity failed in install smoke")


if __name__ == "__main__":
    main()
