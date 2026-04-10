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
    condition_number: float | None = None
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


@dataclass(frozen=True)
class ClosureScoresResult:
    leakage: float
    visible_score: float
    eta: float
    total_curvature: float
    projector: Array
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class ClosureAdaptedObserverResult:
    B: Array
    C: Array
    projector: Array
    scores: ClosureScoresResult
    common_basis: Array
    spectral_energies: Array
    selected_indices: tuple[int, ...]
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class LeakageChannelsResult:
    whitened_perturbation: Array
    visible_operator: Array
    leakage_gram: Array
    singular_values: Array
    visible_channel_basis: Array
    hidden_channel_basis: Array
    coupling: Array
    projector: Array
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class ObserverComparisonResult:
    left_scores: ClosureScoresResult
    right_scores: ClosureScoresResult
    leakage_delta: float
    visible_score_delta: float
    eta_delta: float
    total_curvature_delta: float
    left_dominates: bool
    right_dominates: bool
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class FixedObserverCoordinatesResult:
    phi: Array
    hidden_block: Array
    coupling: Array
    adapted_basis: Array
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class ObserverTransitionResult:
    left: FixedObserverCoordinatesResult
    right: FixedObserverCoordinatesResult
    transform: Array
    a: Array
    b: Array
    c: Array
    d: Array
    c_hat: Array
    d_hat: Array
    right_phi_from_left: Array
    right_hidden_from_left: Array
    right_coupling_from_left: Array
    residual: float
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class ConnectionCurrentResult:
    phi: Array
    hidden_block: Array
    coupling: Array
    coupling_velocity: Array
    current: Array
    forcing: Array
    q: Array
    adapted_basis: Array
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class SupportStratumTransportResult:
    lambda_: Array
    pi: Array
    generator: Array
    lambda_rhs: Array
    pi_rhs: Array
    clock_rate: float
    a_min: float | None
    a_max: float | None
    generator_psd: bool
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class SupportRestartResult:
    event_kind: str
    lambda_before: Array
    lambda_after: Array
    old_basis: Array
    new_basis: Array
    basis_map: Array
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class KernelJetResult:
    order: int | None
    kernel_basis: Array
    gap_basis: Array
    gap_operator: Array
    effective_coefficients: tuple[Array, ...]
    leading_effective: Array
    leading_eigenvalues: Array
    birth_count_forward: int
    death_count_forward: int
    zero_count: int
    event_kind: str
    semisimple: bool
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class SemisimpleEventBlockResult:
    order: int
    side_sign: int
    dimension: int
    leading_matrix: Array
    pole_coefficient: float
    death_like: bool
    birth_like: bool
    clock_log_coefficient: float
    desingularisation_power: float
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class LocalCoupledBirthResult:
    phi: Array
    lift: Array
    hidden_basis: Array
    hidden_metric: Array
    V: Array
    B: Array
    beta: Array
    Q: Array
    observer_tensor: Array
    W: Array
    active_support_basis: Array
    a_cpl: Array
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class SampledIntervalLeakageResult:
    leakage: float
    visible_score: float
    stationarity_residual: float
    projector: Array
    sample_count: int
    weights: Array
    sampled_exact_closure: bool
    metadata: LinearAlgebraMetadata


@dataclass(frozen=True)
class IntervalHessianResult:
    hessian_operator: Array
    spectral_gap: float | None
    rigidity_lower_bound: float | None
    locally_rigid: bool
    metadata: LinearAlgebraMetadata
