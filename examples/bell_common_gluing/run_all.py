from __future__ import annotations

from examples.bell_common_gluing.run_main import main as run_main
from examples.bell_common_gluing.validate import main as run_validate


def main() -> None:
    run_main()
    run_validate()


if __name__ == "__main__":
    main()
