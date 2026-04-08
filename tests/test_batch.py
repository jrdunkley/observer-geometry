from __future__ import annotations

from nomogeo import batch_map

from .helpers import square_value


def test_batch_map_preserves_order_in_process_mode() -> None:
    values = [5, 1, 4, 2]
    result = batch_map(square_value, values, max_workers=2, backend="process")
    assert result == [25, 1, 16, 4]
