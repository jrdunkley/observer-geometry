from __future__ import annotations

import numpy as np

from nomogeo import gaussian_data_processing_contraction, visible_precision


def main() -> None:
    H_a = np.array([[1.0, 0.0], [0.0, 4.0]])
    H_b = np.array([[1.0, 0.0], [0.0, 1.0]])
    coarse = np.array([[1.0, 0.0]])
    rich = (1.0 / np.sqrt(2.0)) * np.array([[1.0, 1.0]])

    print("coarse visible precisions:")
    print(visible_precision(H_a, coarse), visible_precision(H_b, coarse))
    print()
    print("rich visible precisions:")
    print(visible_precision(H_a, rich), visible_precision(H_b, rich))
    print()

    H1 = np.diag([1.0, 2.0, 3.0, 4.0])
    H2 = np.diag([1.5, 2.2, 2.7, 4.5])
    C1 = np.eye(4)
    D = np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 1.0, 0.0]])
    contraction = gaussian_data_processing_contraction(H1, H2, C1, D)

    print("forward KL fine/coarse:", contraction.forward_kl_fine, contraction.forward_kl_coarse)
    print("Hellinger^2 fine/coarse:", contraction.hellinger_sq_fine, contraction.hellinger_sq_coarse)


if __name__ == "__main__":
    main()
