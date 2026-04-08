from __future__ import annotations

import numpy as np

from nomogeo import clock, transport_hidden_load, visible_from_hidden_load


def main() -> None:
    T = np.diag([2.0, 1.0, 0.0])
    lambda_a = np.diag([0.2, 0.5])
    lambda_b = np.diag([0.1, 0.3])
    total = transport_hidden_load(lambda_a, lambda_b)

    print("X(lambda_a):")
    print(visible_from_hidden_load(T, lambda_a, lambda_representation="reduced"))
    print()
    print("transported load:")
    print(total)
    print()
    print("clock additivity:", clock(total), "=", clock(lambda_a) + clock(lambda_b))


if __name__ == "__main__":
    main()
