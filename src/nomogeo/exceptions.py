class NomogeoError(Exception):
    """Base exception for the package."""


class InputValidationError(NomogeoError):
    """Raised when an input violates a mathematical precondition."""


class SupportError(InputValidationError):
    """Raised when support-aware hidden-load assumptions fail."""


class BatchTaskError(NomogeoError):
    """Raised when a batch task fails."""

    def __init__(self, index: int, task: object, original: Exception) -> None:
        self.index = index
        self.task = task
        self.original = original
        message = f"batch task {index} failed for task={task!r}: {original}"
        super().__init__(message)

