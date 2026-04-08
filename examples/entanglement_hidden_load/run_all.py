from __future__ import annotations

from examples.entanglement_hidden_load.run_main import main as run_main
from examples.entanglement_hidden_load.validate import main as run_validate


def main() -> None:
    run_main()
    run_validate()


if __name__ == "__main__":
    main()
