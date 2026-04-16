# Attachment Layer v0 Schema

Status: operational source-of-truth schema for physical attachment cards.

The goal is to convert formula-sheet language into declared observer-geometry input. A card is successful when it states what the module can see, what it cannot see, and which `nomogeo` calls are licensed by the declared object.

## Required Fields

- `card_id`: stable snake-case identifier.
- `schema_version`: currently `attachment_card.v0`.
- `system`: short physical system name.
- `attachment_type`: one of the controlled values below.
- `state_variables`: latent physical coordinates with names and units.
- `sensor_surface`: split into `ideal` and `practical`.
- `control_variables`: controllable physical parameters, or explicit `none_declared`.
- `regime_assumptions`: conditions required before the module object is valid.
- `equation_hook`: formula-sheet equation or compact physical equation.
- `units_and_dimensions`: units for variables, parameters, and matrices.
- `coordinate_note`: short warning about coordinate choices, especially when matrix blocks use mixed physical units.
- `module_object`: declared finite object supplied to the module.
- `derivation_status`: exactness or approximation status of the object.
- `observer_map`: declared observer `C`, or `null` for non-attachable cards.
- `ceiling_or_reference`: declared ceiling/reference for hidden-load calls, or explicit `none`.
- `allowed_nomogeo_calls`: exact list of allowed calls for this card.
- `nomogeo_readout`: requested readout fields.
- `positive_scope`: one sentence describing what this card can expose.
- `boundary`: one sentence describing what this card cannot expose.
- `failure_mode`: concrete conditions that invalidate the attachment.
- `next_measurement`: next useful measurement or control action.
- `default_parameters`: runnable default values.
- `declared_objects`: runnable matrices such as `H`, `C`, and optional `T`.

## Attachment Types

- `exact_quadratic`
- `exact_quadratic_with_declared_ceiling`
- `exact_special_sector`
- `local_quadratic_approximation`
- `declared_fluctuation_model`
- `weighted_family`
- `support_path`
- `not_yet_attachable`

`sensor_comparison` is intentionally not an attachment type in v0. It is a task mode built on top of cards with concrete attachment types.

## Runner Rule

A card may only execute module calls licensed by `allowed_nomogeo_calls` and supported by its declared object fields.

For example:

- `visible_precision` requires declared `H` and `C`.
- `hidden_load` requires declared `H`, `C`, and `T`; `T` must be a ceiling/reference for the computed visible precision, and the runner checks `T - Phi >= 0` before calling the module surface.
- `local_visible_calculus` requires declared `H`, `C`, and `Delta`.
- `not_yet_attachable` must refuse theorem-local kernel calls.

This rule is the core safeguard of the attachment layer.
