from __future__ import annotations

import numpy as np

from nomogeo import local_visible_calculus


def main() -> None:
    H = np.array([[5.0, 1.0, 0.2], [1.0, 3.5, 0.4], [0.2, 0.4, 2.8]])
    C = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    Delta = np.array([[0.6, 0.2, 0.1], [0.2, 0.3, 0.15], [0.1, 0.15, 0.2]])

    result = local_visible_calculus(H, C, Delta)
    print("V:")
    print(result.V)
    print()
    print("Q:")
    print(result.Q)
    print()
    print("determinant-curvature split:", result.det_split)


if __name__ == "__main__":
    main()
