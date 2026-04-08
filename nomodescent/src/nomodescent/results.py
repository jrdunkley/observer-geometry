from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .audit import AuditReport
from .specs import ObserverSpec


@dataclass(frozen=True)
class ObstructionCertificate:
    kind: str
    exact: bool
    summary: str
    details: dict[str, float | int | str | list[float] | dict[str, object]]


@dataclass(frozen=True)
class DescentResult:
    classification: str
    exact: bool
    factor_map: np.ndarray | None
    residuals: dict[str, float]
    certificates: tuple[ObstructionCertificate, ...]
    common_covariance: np.ndarray | None
    common_precision: np.ndarray | None
    audit: AuditReport
    details: dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class RefinementSearchResult:
    winner: str | None
    objective: str
    scores: dict[str, float]
    certificates: tuple[ObstructionCertificate, ...]
    audit: AuditReport
    details: dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class CommonRefinementResult:
    classification: str
    exact: bool
    common_refinement: ObserverSpec | None
    factor_maps: dict[str, np.ndarray]
    residuals: dict[str, float]
    certificates: tuple[ObstructionCertificate, ...]
    audit: AuditReport
    details: dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class FalseCollapseResult:
    classification: str
    exact: bool
    coarse_summary_label: str
    coarse_summary_max_gap: float
    coarse_summary_agreement: bool
    coarse_certified_compatible: bool
    fine_result: DescentResult
    audit: AuditReport
    details: dict[str, object] = field(default_factory=dict)
