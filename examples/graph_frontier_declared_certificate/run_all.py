from __future__ import annotations

from examples.graph_frontier_declared_certificate.run_main import main as run_main
from examples.graph_frontier_declared_certificate.validate import main as run_validate


def main() -> None:
    run_main()
    run_validate()


if __name__ == "__main__":
    main()
