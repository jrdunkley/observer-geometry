from __future__ import annotations

import itertools
from collections.abc import Callable, Sequence
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from typing import Any

from .exceptions import BatchTaskError, InputValidationError


def _apply_task(func: Callable[..., Any], task: Any) -> Any:
    if isinstance(task, dict):
        return func(**task)
    if isinstance(task, tuple):
        return func(*task)
    return func(task)


def _apply_task_indexed(index: int, func: Callable[..., Any], task: Any) -> Any:
    try:
        return _apply_task(func, task)
    except Exception as exc:
        raise BatchTaskError(index, task, exc) from exc


def batch_map(
    func: Callable[..., Any],
    tasks: Sequence[Any],
    *,
    max_workers: int | None = None,
    backend: str = "process",
) -> list[Any]:
    """Apply a pure function to independent tasks with deterministic ordering."""
    if backend not in {"serial", "thread", "process"}:
        raise InputValidationError("backend must be 'serial', 'thread', or 'process'")

    indexed_tasks = list(enumerate(tasks))
    if backend == "serial" or max_workers == 1 or len(indexed_tasks) <= 1:
        return [_apply_task_indexed(index, func, task) for index, task in indexed_tasks]

    executor_cls = ThreadPoolExecutor if backend == "thread" else ProcessPoolExecutor
    try:
        with executor_cls(max_workers=max_workers) as executor:
            iterator = executor.map(
                _apply_task_indexed,
                (index for index, _task in indexed_tasks),
                itertools.repeat(func),
                (task for _index, task in indexed_tasks),
            )
            return list(iterator)
    except PermissionError:
        if backend != "process":
            raise
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            iterator = executor.map(
                _apply_task_indexed,
                (index for index, _task in indexed_tasks),
                itertools.repeat(func),
                (task for _index, task in indexed_tasks),
            )
            return list(iterator)
