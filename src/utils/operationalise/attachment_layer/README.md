# Attachment Layer v0

This folder is the first operational attachment grammar for physical systems.

It does not add theory to `nomogeo`. It converts simple physical formula-sheet objects into declared finite inputs that the module may or may not be allowed to touch.

## Contents

- `schema.md`: required card fields and allowed attachment types.
- `cards/*.json`: source-of-truth attachment cards.
- `run_card.py`: strict runner for cards.
- `outputs/*.json`: generated readouts from the runner.
- `outputs/*.md`: generated human summaries from the JSON cards and readouts.
- `ap_hook_digest.md`: first extraction of promising hooks from AP/A-level material.
- `core_coverage_map.md`: breadth map across the supplied formula-sheet domains.
- `source_exhaustion_log.md`: source-by-source record of which PDF material has been routed.
- `candidate_backlog.json`: machine-readable backlog of candidate, deferred, refusal, and carded hooks.
- `two_track_programme.md`: programme memo for diagnosing non-integrating hooks and signing off good cards.
- `card_signoff_v0.md`: sign-off status for the current card set and analogue sweep result.
- `analogue_sweep.py`: deterministic coupled spring / coupled LC stress check.
- `non_integration_diagnostics.json`: machine-readable route classification for every non-carded hook.
- `non_integration_diagnostics.md`: human summary of why non-carded hooks are not v0 static cards.
- `validate_operationalise_metadata.py`: checks backlog, diagnostics, sign-off, and sweep consistency.
- `failure_probe_suite.py`: representative numerical probes for failure classes without promoting them to cards.
- `next_schema_decisions.md`: ranked next schema decisions derived from the failure probes.
- `local_positive_curvature_admissibility.md`: theoretical admissibility memo for declared local positive curvature.
- `declared_local_positive_contact_v0.md` / `declared_local_positive_contact_v0.json`: scientific synthesis of the current contact doctrine, route ontology, gold quartet, and mandatory refusals.
- `main_module_contact_bridge_v0.md` / `main_module_contact_bridge_v0.json`: wrap-ready bridge back to the main module, including handoff types, boundaries, and next-phase recommendation.
- `gold_quartet_demonstration_v0.md` / `gold_quartet_demonstration_v0.json`: compact demonstration bundle checked against generated card and empirical-probe outputs.
- `gold_quartet_demonstration.py`: generator for the compact gold-quartet bundle from the current card and empirical-probe outputs.
- `local_curvature_proofs.py`: RC decay and Beer-Lambert proof-of-admissibility generator.
- `local_curvature_stress_tests.py`: reparameterization, noise-scaling, and negative-control checks for the local curvature proofs.
- `measurement_inference_source_digest.md`: first digest of the GUM/NIST measurement-model, process-control, process-modeling, and eigenbasis sources.
- `measurement_compilation_pipeline.md`: non-promoted compiler design note for measurement equations, observation equations, calibration/analysis functions, process gates, and finite modal routes.
- `metric_modal_route_plan.md`: scientific plan for separating metric curvature, stiffness curvature, and finite modal readouts.
- `fifteen_point_operationalisation_plan.md`: banked five-step spine expanded into a 15-point execution plan.
- `admissibility_contract_v0.md` / `admissibility_contract_v0.json`: route contract for admissible objects and mandatory refusals.
- `measurement_route_backlog.json`: machine-readable backlog for the measurement/inference routes.
- `measurement_equation_linearisation_proofs.py`: generated proof packet for local linearisation of measurement equations, including a negative control.
- `observation_equation_fisher_proofs.py`: generated proof packet for Fisher precision from observation equations, including rank and zero-sensitivity refusals.
- `finite_modal_eigenbasis_proofs.py`: generated proof packet for finite wave/modal attachment, including raw-wave and zero-stiffness refusals.
- `calibration_analysis_function_proofs.py`: generated proof packet for calibrated analysis-function inversion, including extrapolation and rank refusals.
- `process_control_gate_proofs.py`: generated proof packet for check-standard process gates, including drift and missing-baseline refusals.
- `empirical_dataset_protocol_v0.md`: acquisition protocol for small empirical datasets that can test the measurement and observation routes.
- `dataset_templates/*.json` and `dataset_templates/*.csv`: source templates for RC decay traces and Beer-Lambert calibration datasets.
- `online_dataset_acquisition_log.md`: source and admission log for external/manual empirical datasets.
- `beer_lambert_uvvis_probe.py`: UV-Vis Beer-Lambert empirical probe over the admitted Edinburgh Rhodamine dataset.
- `rc_decay_q7_probe.py`: RC decay empirical probe over the manually supplied Q7 capacitor voltage dataset.

## Run

From this folder:

```powershell
python run_card.py --all --markdown
```

If `nomogeo` is not installed in the active environment, pass the source path:

```powershell
python run_card.py --all --markdown --nomogeo-src C:\path\to\observer_geometry\src
```

or set:

```powershell
$env:NOMOGEO_SRC = "C:\path\to\observer_geometry\src"
```

Run a single card:

```powershell
python run_card.py --card mass_spring_coupled_observe_one --markdown
```

Validate card and runner invariants:

```powershell
python validate_cards.py --nomogeo-src C:\path\to\observer_geometry\src
```

Run the coupled spring / coupled LC sign-off sweep:

```powershell
python analogue_sweep.py --nomogeo-src C:\path\to\observer_geometry\src
```

Validate operational metadata:

```powershell
python validate_operationalise_metadata.py
```

Run representative non-integration probes:

```powershell
python failure_probe_suite.py
```

Generate local positive curvature admissibility proofs:

```powershell
python local_curvature_proofs.py
```

Stress-test the local curvature proofs:

```powershell
python local_curvature_stress_tests.py
```

Generate measurement-equation local-linearisation proofs:

```powershell
python measurement_equation_linearisation_proofs.py
```

Generate observation-equation Fisher proofs:

```powershell
python observation_equation_fisher_proofs.py
```

Generate finite-modal eigenbasis proofs:

```powershell
python finite_modal_eigenbasis_proofs.py
```

Generate calibration-analysis function proofs:

```powershell
python calibration_analysis_function_proofs.py
```

Generate process-control gate proofs:

```powershell
python process_control_gate_proofs.py
```

Run the admitted empirical probes:

```powershell
python beer_lambert_uvvis_probe.py
python rc_decay_q7_probe.py
```

Generate the compact gold-quartet demonstration:

```powershell
python gold_quartet_demonstration.py
```

## Rule

A card may only execute calls listed in `allowed_nomogeo_calls` and supported by its declared object fields.

For hidden-load cards, the runner checks that the declared ceiling dominates the computed visible precision before calling the module surface.

The first expected success is boring: a repeatable, validated conversion from formula-sheet language into declared observer-geometry input.
