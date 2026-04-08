from __future__ import annotations


class EvidenceInputError(ValueError):
    """Raised when an evidence bundle or item is malformed."""


class EvidenceAssemblyError(ValueError):
    """Raised when evidence cannot be assembled honestly into a problem spec."""
