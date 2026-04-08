from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]


@dataclass(frozen=True)
class LinearAlgebraMetadata:
    atol: float
    rtol: float
    rank_tol: float
    method: str
    ambient_dim: int
    support_rank: int
    visible_dim: int | None = None
    support_restricted: bool = False
    notes: tuple[str, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class VisibleGeometryResult:
    phi: Array
    lift: Array
    projector: Array
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class LocalCalculusResult:
    phi: Array
    lift: Array
    projector: Array
    V: Array
    Q: Array
    det_split: float
    active_support: Array
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class HiddenLoadResult:
    lambda_: Array
    reduced_lambda: Array
    X: Array
    ceiling: Array
    support_basis: Array
    rank: int
    clock: float
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class RealisationResult:
    matrix: Array
    factor: Array
    support_basis: Array
    rank: int
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class DVBridgeResult:
    h_dv: Array
    delta_dv: Array
    gram_factor: Array
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class GaussianContractionResult:
    forward_kl_fine: float
    forward_kl_coarse: float
    reverse_kl_fine: float
    reverse_kl_coarse: float
    hellinger_sq_fine: float
    hellinger_sq_coarse: float
    bhattacharyya_fine: float
    bhattacharyya_coarse: float
    metadata: LinearAlgebraMetadata

