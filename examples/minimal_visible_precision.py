from __future__ import annotations

import numpy as np

from nomogeo import canonical_lift, visible_precision


def main() -> None:
    H = np.array([[4.0, 1.0, 0.2], [1.0, 3.0, 0.5], [0.2, 0.5, 2.5]])
    C = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])

    phi = visible_precision(H, C)
    lift = canonical_lift(H, C)

    print("Phi_C(H):")
    print(phi)
    print()
    print("L_{C,H}:")
    print(lift)


if __name__ == "__main__":
    main()
