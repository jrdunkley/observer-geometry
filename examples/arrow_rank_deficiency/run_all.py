from __future__ import annotations

from examples.arrow_rank_deficiency.run_main import main as run_main
from examples.arrow_rank_deficiency.validate import main as run_validate


def main() -> None:
    run_main()
    run_validate()


if __name__ == "__main__":
    main()
