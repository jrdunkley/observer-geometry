# Admissibility Contract v0

Status: route contract, not a promoted attachment schema.

This contract banks the current operational lesson:

```text
formula surface -> route -> declared finite positive object -> allowed surface or refusal
```

A formula, sensor law, calibration curve, or wave equation is not admitted because it is familiar. It is admitted only after it declares what finite object the module can see and what it still cannot see.

The higher-level scientific synthesis of this rule is banked in `declared_local_positive_contact_v0.md` and `declared_local_positive_contact_v0.json`.

## Routes

| Route | Maturity | Emitted Object | Allowed Surfaces |
| --- | --- | --- | --- |
| `exact_quadratic_storage` | `carded` | finite storage Hessian or precision object | `visible_precision` |
| `local_quadratic_storage` | `carded` | local Hessian at declared operating point | `visible_precision` after operating point declaration |
| `declared_ceiling_hidden_load` | `carded` | storage Hessian plus visible-block ceiling/reference | `visible_precision`, `hidden_load` after ceiling dominance check |
| `measurement_equation_local_linearisation` | `proofed` | local covariance or precision from sensitivity/Jacobian | local positive curvature record only |
| `observation_equation_fisher` | `proofed` | scalar or matrix Fisher precision | local positive curvature record only |
| `finite_modal_eigenbasis` | `proofed` | finite modal mass and stiffness matrices | finite modal positive curvature or observer-map record |
| `calibration_analysis_function` | `proofed_minimal` | analysis function with local uncertainty and validity range | analysis-function record; local calibrated-measurand precision only after uncertainty propagation |
| `process_control_gate` | `proofed_gate_only` | measurement-process validity decision | gate or refusal only |

## Mandatory Refusal Pattern

Every route has mandatory refusals. The shared refusal is:

```text
hidden_load without a declared ceiling/reference is refused.
```

Additional important refusals:

- measurement equations refuse missing uncertainty models and failed first-order linearisation;
- observation equations refuse missing noise, zero sensitivity, rank deficiency, and residual-gate failure;
- finite modal routes refuse raw infinite wave formulas and missing truncation;
- calibration analysis refuses extrapolation and inverse use without validity range;
- process-control gates refuse theorem-local kernel calls, missing-baseline gates, and out-of-control measurement processes.

## Promotion Rule

Do not promote a broad measurement attachment schema until these all pass together:

- carded storage routes;
- measurement-equation local linearisation proofs;
- observation-equation Fisher proofs;
- finite-modal eigenbasis proofs;
- minimal calibration-analysis proof;
- process-control gate proof.

The contract is machine-readable in `admissibility_contract_v0.json` and validated by `validate_operationalise_metadata.py`.
