from __future__ import annotations

import numpy as np

from nomogeo import dv_bridge, local_visible_calculus, visible_precision


def main() -> None:
    H0 = np.array([[4.0, 0.8, 0.0], [0.8, 3.0, 0.4], [0.0, 0.4, 2.2]])
    Jhat = np.array([[0.0, 1.0, -0.4], [-1.0, 0.0, 0.7], [0.4, -0.7, 0.0]])
    C = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])

    bridge = dv_bridge(H0, Jhat)
    local = local_visible_calculus(H0, C, bridge.delta_dv)

    eps = 1e-2
    phi_shift = visible_precision(H0 + (eps * eps) * bridge.delta_dv, C) - visible_precision(H0, C)

    print("Delta_DV:")
    print(bridge.delta_dv)
    print()
    print("quadratic visible onset estimate:")
    print(phi_shift / (eps * eps))
    print()
    print("predicted V:")
    print(local.V)


if __name__ == "__main__":
    main()
