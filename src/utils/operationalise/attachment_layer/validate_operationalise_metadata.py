from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parent


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _fail(message: str) -> None:
    raise AssertionError(message)


def _ids(items: list[dict[str, Any]]) -> set[str]:
    return {str(item["id"]) for item in items}


def _assert_close(actual: Any, expected: Any, label: str, *, relative_tolerance: float = 1.0e-10) -> None:
    actual_value = float(actual)
    expected_value = float(expected)
    scale = max(1.0, abs(expected_value))
    if abs(actual_value - expected_value) > relative_tolerance * scale:
        _fail(f"{label}: expected {expected_value}, got {actual_value}")


def validate_backlog(backlog: dict[str, Any]) -> None:
    if backlog.get("schema_version") != "attachment_candidate_backlog.v0":
        _fail("candidate_backlog.json schema_version mismatch")
    items = backlog.get("items")
    if not isinstance(items, list) or not items:
        _fail("candidate_backlog.json must contain non-empty items")
    seen: set[str] = set()
    for item in items:
        item_id = item.get("id")
        if not isinstance(item_id, str) or not item_id:
            _fail("backlog item missing id")
        if item_id in seen:
            _fail(f"duplicate backlog id: {item_id}")
        seen.add(item_id)
        for field in ("domain", "status", "source_surfaces", "attachment_test", "module_object_candidate", "sensor_surface", "control_variables", "next_action"):
            if field not in item:
                _fail(f"{item_id}: missing backlog field {field}")


def validate_diagnostics(backlog: dict[str, Any], diagnostics: dict[str, Any]) -> None:
    if diagnostics.get("schema_version") != "non_integration_diagnostics.v0":
        _fail("non_integration_diagnostics.json schema_version mismatch")
    backlog_items = backlog["items"]
    expected = {
        item["id"]
        for item in backlog_items
        if item["status"] not in {"carded", "carded_refusal"}
    }
    diagnostic_items = diagnostics.get("items")
    if not isinstance(diagnostic_items, list):
        _fail("diagnostics items must be a list")
    actual = _ids(diagnostic_items)
    if expected != actual:
        _fail(f"diagnostic ids do not match non-carded backlog ids; missing={sorted(expected - actual)} extra={sorted(actual - expected)}")

    route_classes = set(diagnostics.get("route_classes", {}))
    if not route_classes:
        _fail("diagnostics route_classes must be non-empty")
    for item in diagnostic_items:
        item_id = item["id"]
        if item.get("route_class") not in route_classes:
            _fail(f"{item_id}: unknown route_class {item.get('route_class')!r}")
        for field in ("likely_module_route", "why_not_static_card", "missing_declarations", "next_proof_obligation", "gap_assessment"):
            if field not in item:
                _fail(f"{item_id}: missing diagnostic field {field}")
        if not isinstance(item["likely_module_route"], list) or not item["likely_module_route"]:
            _fail(f"{item_id}: likely_module_route must be non-empty")
        if not isinstance(item["missing_declarations"], list):
            _fail(f"{item_id}: missing_declarations must be a list")


def validate_sweep(sweep: dict[str, Any]) -> None:
    if sweep.get("schema_version") != "attachment_sweep.v0":
        _fail("sweep schema_version mismatch")
    rows = sweep.get("rows")
    if not isinstance(rows, list) or not rows:
        _fail("sweep rows must be non-empty")
    fields = sweep.get("monotone_checks")
    if not isinstance(fields, list) or not fields:
        _fail("sweep monotone_checks must be non-empty")
    by_family: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        family = row.get("family")
        if not isinstance(family, str) or not family:
            _fail("sweep row missing family")
        by_family.setdefault(family, []).append(row)
        if row.get("executed_calls") != ["visible_precision", "hidden_load"]:
            _fail(f"{family}: sweep row did not execute visible_precision and hidden_load")
    for family, family_rows in by_family.items():
        family_rows.sort(key=lambda row: float(row["coupling"]))
        for field in fields:
            values = [float(row[field]) for row in family_rows]
            if any(b <= a for a, b in zip(values, values[1:])):
                _fail(f"{family}: sweep field {field} is not strictly increasing: {values}")


def validate_failure_probes(backlog: dict[str, Any], diagnostics: dict[str, Any], probes: dict[str, Any]) -> None:
    if probes.get("schema_version") != "failure_probe_suite.v0":
        _fail("failure_probe_suite.json schema_version mismatch")
    probe_items = probes.get("probes")
    if not isinstance(probe_items, list) or not probe_items:
        _fail("failure probes must be a non-empty list")

    diagnostic_by_id = {item["id"]: item for item in diagnostics["items"]}
    backlog_by_id = {item["id"]: item for item in backlog["items"]}
    route_classes = set(diagnostics.get("route_classes", {}))
    for probe in probe_items:
        probe_id = probe.get("id")
        route_class = probe.get("route_class")
        if route_class not in route_classes:
            _fail(f"{probe_id}: unknown probe route_class")
        if probe_id in diagnostic_by_id:
            if route_class != diagnostic_by_id[probe_id]["route_class"]:
                _fail(f"{probe_id}: probe route_class does not match diagnostics")
        elif backlog_by_id.get(probe_id, {}).get("status") not in {"carded", "carded_refusal"}:
            _fail(f"failure probe {probe_id!r} is neither diagnosed nor a carded item")
        for field in (
            "surface_formula",
            "refusal_without_extra_declarations",
            "minimal_extra_declarations",
            "declared_object_after_extra_declarations",
            "numeric_readout",
            "lesson",
        ):
            if field not in probe:
                _fail(f"{probe_id}: missing failure probe field {field}")
        if not isinstance(probe["minimal_extra_declarations"], list) or not probe["minimal_extra_declarations"]:
            _fail(f"{probe_id}: minimal_extra_declarations must be non-empty")
        if not isinstance(probe["numeric_readout"], dict) or not probe["numeric_readout"]:
            _fail(f"{probe_id}: numeric_readout must be non-empty")


def validate_local_curvature_proofs(proofs: dict[str, Any]) -> None:
    if proofs.get("schema_version") != "local_positive_curvature_proofs.v0":
        _fail("local_positive_curvature_proofs.json schema_version mismatch")
    proof_items = proofs.get("proofs")
    if not isinstance(proof_items, list) or len(proof_items) != 2:
        _fail("local positive curvature proofs must contain exactly two proofs")
    expected_ids = {"rc_decay_local_curvature", "beer_lambert_local_curvature"}
    actual_ids = {proof.get("id") for proof in proof_items}
    if actual_ids != expected_ids:
        _fail(f"local curvature proof ids mismatch: {actual_ids}")
    for proof in proof_items:
        proof_id = proof["id"]
        if proof.get("allowed_route") != "local_positive_curvature":
            _fail(f"{proof_id}: allowed_route must be local_positive_curvature")
        positive = proof.get("positive_object")
        if not isinstance(positive, dict):
            _fail(f"{proof_id}: positive_object must be present")
        value = float(positive.get("value", 0.0))
        if value <= 0.0:
            _fail(f"{proof_id}: positive_object value must be positive")
        if "static storage Hessian" not in proof.get("refused_routes", []):
            _fail(f"{proof_id}: must refuse static storage Hessian route")
        if "hidden_load without a ceiling" not in proof.get("refused_routes", []):
            _fail(f"{proof_id}: must refuse hidden_load without a ceiling")


def validate_local_curvature_stress_tests(stress: dict[str, Any]) -> None:
    if stress.get("schema_version") != "local_curvature_stress_tests.v0":
        _fail("local_curvature_stress_tests.json schema_version mismatch")
    tests = stress.get("tests")
    if not isinstance(tests, list) or len(tests) != 2:
        _fail("local curvature stress tests must contain exactly two tests")
    expected_ids = {"rc_decay_local_curvature_stress", "beer_lambert_local_curvature_stress"}
    actual_ids = {test.get("id") for test in tests}
    if actual_ids != expected_ids:
        _fail(f"local curvature stress test ids mismatch: {actual_ids}")
    for test in tests:
        test_id = test["id"]
        if test.get("status") != "passed":
            _fail(f"{test_id}: stress status must be passed")
        values = [float(value) for value in test.get("fisher_by_noise", [])]
        if len(values) < 2:
            _fail(f"{test_id}: fisher_by_noise must have at least two values")
        if any(b >= a for a, b in zip(values, values[1:])):
            _fail(f"{test_id}: Fisher values should strictly decrease as noise increases")
        controls = test.get("negative_controls", {})
        if float(controls.get("zero_sensitivity_fisher", -1.0)) != 0.0:
            _fail(f"{test_id}: zero sensitivity control should have zero Fisher")
        for control_name in ("zero_noise_refused", "hidden_load_route_refused_without_ceiling", "static_storage_route_refused"):
            if controls.get(control_name) is not True:
                _fail(f"{test_id}: control {control_name} must be true")


def validate_signoff(signoff_path: Path) -> None:
    text = signoff_path.read_text(encoding="utf-8")
    for required in ("good_v0_baseline", "good_v0_storage", "good_v0_metric", "good_v0_finite_modal", "good_v0_modal_observer", "good_v0_modal_reference", "good_v0_bridge", "good_v0_regime", "good_v0_refusal"):
        if required not in text:
            _fail(f"card_signoff_v0.md missing {required}")
    for card_id in (
        "mass_spring_single",
        "lc_resonator_single",
        "linear_elastic_bar",
        "capacitor_charge_coordinate",
        "kinetic_and_rotational_energy_metrics",
        "string_fixed_end_three_mode",
        "string_point_sensor_modal_observer",
        "string_modal_hidden_load_with_declared_ceiling",
        "acoustic_tube_three_mode_pressure_observer",
        "mass_spring_coupled_observe_one",
        "coupled_lc_resonators_observe_one_node",
        "rlc_resonator_linear_storage",
        "coupled_rlc_resonators_observe_one_node",
        "pendulum_small_angle",
        "pendulum_full",
        "ideal_gas_law_boundary",
    ):
        if card_id not in text:
            _fail(f"card_signoff_v0.md missing {card_id}")


def validate_measurement_route_backlog(backlog: dict[str, Any]) -> None:
    if backlog.get("schema_version") != "measurement_route_backlog.v0":
        _fail("measurement_route_backlog.json schema_version mismatch")
    items = backlog.get("items")
    if not isinstance(items, list) or not items:
        _fail("measurement_route_backlog.json must contain non-empty items")
    allowed_routes = {
        "measurement_equation_local_linearisation",
        "observation_equation_fisher",
        "calibration_analysis_function",
        "process_control_gate",
        "finite_modal_eigenbasis",
    }
    seen: set[str] = set()
    for item in items:
        item_id = item.get("id")
        if not isinstance(item_id, str) or not item_id:
            _fail("measurement route item missing id")
        if item_id in seen:
            _fail(f"duplicate measurement route id: {item_id}")
        seen.add(item_id)
        if item.get("route") not in allowed_routes:
            _fail(f"{item_id}: unknown measurement route {item.get('route')!r}")
        for field in (
            "status",
            "source_basis",
            "declared_object_after_admission",
            "required_declarations",
            "refusal_conditions",
        ):
            if field not in item:
                _fail(f"{item_id}: missing measurement route field {field}")
        if not isinstance(item["source_basis"], list) or not item["source_basis"]:
            _fail(f"{item_id}: source_basis must be non-empty")
        if not isinstance(item["required_declarations"], list) or not item["required_declarations"]:
            _fail(f"{item_id}: required_declarations must be non-empty")
        if not isinstance(item["refusal_conditions"], list) or not item["refusal_conditions"]:
            _fail(f"{item_id}: refusal_conditions must be non-empty")


def validate_measurement_equation_proofs(proofs: dict[str, Any]) -> None:
    if proofs.get("schema_version") != "measurement_equation_linearisation_proofs.v0":
        _fail("measurement_equation_linearisation_proofs.json schema_version mismatch")
    proof_items = proofs.get("proofs")
    if not isinstance(proof_items, list) or len(proof_items) != 3:
        _fail("measurement equation linearisation proofs must contain exactly three proofs")
    expected_ids = {
        "beer_lambert_concentration_measurement_equation",
        "cadmium_calibration_standard_measurement_equation",
        "bromine_abundance_linearisation_negative_control",
    }
    actual_ids = {proof.get("id") for proof in proof_items}
    if actual_ids != expected_ids:
        _fail(f"measurement equation proof ids mismatch: {actual_ids}")
    for proof in proof_items:
        proof_id = proof["id"]
        if proof.get("route") != "measurement_equation_local_linearisation":
            _fail(f"{proof_id}: route must be measurement_equation_local_linearisation")
        if proof_id.endswith("negative_control"):
            if proof.get("status") != "first_order_refused":
                _fail(f"{proof_id}: negative control must refuse first-order route")
            if float(proof.get("first_order_variance", -1.0)) != 0.0:
                _fail(f"{proof_id}: first-order variance should be zero in the negative control")
            if float(proof.get("second_order_standard_uncertainty_for_normal_local_model", 0.0)) <= 0.0:
                _fail(f"{proof_id}: second-order standard uncertainty should be positive")
            if "first_order local precision" not in proof.get("refused_routes", []):
                _fail(f"{proof_id}: must explicitly refuse first_order local precision")
        else:
            if proof.get("status") != "admissible_local_precision":
                _fail(f"{proof_id}: expected admissible local precision")
            if float(proof.get("local_variance", 0.0)) <= 0.0:
                _fail(f"{proof_id}: local_variance must be positive")
            if float(proof.get("local_precision", 0.0)) <= 0.0:
                _fail(f"{proof_id}: local_precision must be positive")
            if float(proof.get("max_sensitivity_error", 1.0)) > 1.0e-5:
                _fail(f"{proof_id}: sensitivity check too large")
            if "hidden_load without a ceiling" not in proof.get("refused_routes", []):
                _fail(f"{proof_id}: must refuse hidden_load without a ceiling")


def validate_observation_equation_proofs(proofs: dict[str, Any]) -> None:
    if proofs.get("schema_version") != "observation_equation_fisher_proofs.v0":
        _fail("observation_equation_fisher_proofs.json schema_version mismatch")
    proof_items = proofs.get("proofs")
    if not isinstance(proof_items, list) or len(proof_items) != 3:
        _fail("observation equation Fisher proofs must contain exactly three proofs")
    expected_ids = {
        "rc_decay_tau_observation_equation",
        "repeated_current_mean_observation_equation",
        "load_cell_linear_calibration_observation_equation",
    }
    actual_ids = {proof.get("id") for proof in proof_items}
    if actual_ids != expected_ids:
        _fail(f"observation equation proof ids mismatch: {actual_ids}")
    for proof in proof_items:
        proof_id = proof["id"]
        if proof.get("route") != "observation_equation_fisher":
            _fail(f"{proof_id}: route must be observation_equation_fisher")
        if proof.get("allowed_route") != "local_positive_curvature":
            _fail(f"{proof_id}: allowed_route must be local_positive_curvature")
        if proof.get("residual_gate", {}).get("passed_simple_gate") is not True:
            _fail(f"{proof_id}: residual gate did not pass")
        if "hidden_load without a ceiling" not in proof.get("refused_routes", []):
            _fail(f"{proof_id}: must refuse hidden_load without a ceiling")
        fisher = proof.get("fisher_object")
        if not isinstance(fisher, dict):
            _fail(f"{proof_id}: missing fisher_object")
        if "value" in fisher:
            if float(fisher.get("value", 0.0)) <= 0.0:
                _fail(f"{proof_id}: scalar Fisher value must be positive")
            if float(fisher.get("local_variance", 0.0)) <= 0.0:
                _fail(f"{proof_id}: scalar local variance must be positive")
        else:
            spd = fisher.get("spd_check", {})
            eigenvalues = [float(item) for item in spd.get("eigenvalues", [])]
            if len(eigenvalues) != 2 or any(value <= 0.0 for value in eigenvalues):
                _fail(f"{proof_id}: Fisher matrix must have two positive eigenvalues")
            if float(spd.get("condition_number", 0.0)) <= 1.0:
                _fail(f"{proof_id}: condition_number should be present and greater than one")

    controls = proofs.get("negative_controls")
    if not isinstance(controls, list) or len(controls) != 2:
        _fail("observation equation Fisher proofs must contain exactly two negative controls")
    control_by_id = {control.get("id"): control for control in controls}
    zero = control_by_id.get("rc_decay_zero_sensitivity_grid_refusal")
    if zero is None or zero.get("status") != "refused" or float(zero.get("computed_fisher", -1.0)) != 0.0:
        _fail("RC zero sensitivity negative control did not refuse with zero Fisher")
    rank = control_by_id.get("load_cell_rank_deficient_design_refusal")
    if rank is None or rank.get("status") != "refused":
        _fail("load-cell rank deficient negative control missing or not refused")
    rank_eigs = [float(item) for item in rank.get("fisher_eigenvalues", [])]
    if len(rank_eigs) != 2 or min(abs(value) for value in rank_eigs) > 1.0e-8:
        _fail("load-cell rank deficient negative control should have a near-zero Fisher eigenvalue")


def validate_finite_modal_proofs(proofs: dict[str, Any]) -> None:
    if proofs.get("schema_version") != "finite_modal_eigenbasis_proofs.v0":
        _fail("finite_modal_eigenbasis_proofs.json schema_version mismatch")
    proof_items = proofs.get("proofs")
    if not isinstance(proof_items, list) or len(proof_items) != 3:
        _fail("finite modal eigenbasis proofs must contain exactly three proofs")
    proof_by_id = {proof.get("id"): proof for proof in proof_items}
    string = proof_by_id.get("string_fixed_end_three_mode")
    if string is None or string.get("status") != "admissible_finite_modal_family":
        _fail("string_fixed_end_three_mode proof missing or not admissible")
    module_objects = string.get("module_objects", {})
    for check_name in ("mass_metric_spd", "stiffness_spd"):
        eigenvalues = [float(item) for item in module_objects.get(check_name, {}).get("eigenvalues", [])]
        if len(eigenvalues) != 3 or any(value <= 0.0 for value in eigenvalues):
            _fail(f"string_fixed_end_three_mode: {check_name} must have three positive eigenvalues")
    frequencies = [float(item) for item in string.get("modal_readout", {}).get("frequencies_hz", [])]
    if len(frequencies) != 3 or any(b <= a for a, b in zip(frequencies, frequencies[1:])):
        _fail("string_fixed_end_three_mode: modal frequencies should be strictly increasing")
    if "hidden_load without a ceiling" not in string.get("refused_routes", []):
        _fail("string_fixed_end_three_mode must refuse hidden_load without a ceiling")

    projection = proof_by_id.get("string_point_sensor_two_mode_projection")
    if projection is None or projection.get("status") != "observer_map_recorded_not_hidden_load":
        _fail("string point sensor projection proof missing or wrong status")
    if projection.get("rank") != 1:
        _fail("string point sensor projection rank should be one")
    if "hidden_load without a ceiling" not in projection.get("refused_routes", []):
        _fail("string point sensor projection must refuse hidden_load without a ceiling")

    acoustic = proof_by_id.get("acoustic_tube_three_mode_pressure_observer")
    if acoustic is None or acoustic.get("status") != "admissible_finite_modal_pressure_observer":
        _fail("acoustic tube pressure observer proof missing or wrong status")
    acoustic_objects = acoustic.get("module_objects", {})
    for check_name in ("mass_metric_spd", "stiffness_spd"):
        eigenvalues = [float(item) for item in acoustic_objects.get(check_name, {}).get("eigenvalues", [])]
        if len(eigenvalues) != 3 or any(value <= 0.0 for value in eigenvalues):
            _fail(f"acoustic_tube_three_mode_pressure_observer: {check_name} must have three positive eigenvalues")
    pressure_precision = np.asarray(acoustic_objects.get("visible_pressure_precision", []), dtype=float)
    if pressure_precision.shape != (1, 1) or float(pressure_precision[0, 0]) <= 0.0:
        _fail("acoustic_tube_three_mode_pressure_observer: pressure visible precision must be positive scalar")
    if "raw acoustic impedance" not in acoustic.get("refused_routes", []):
        _fail("acoustic_tube_three_mode_pressure_observer must refuse raw acoustic impedance")
    if "hidden_load without a ceiling" not in acoustic.get("refused_routes", []):
        _fail("acoustic_tube_three_mode_pressure_observer must refuse hidden_load without a ceiling")

    controls = proofs.get("negative_controls")
    if not isinstance(controls, list) or len(controls) != 2:
        _fail("finite modal eigenbasis proofs must contain exactly two negative controls")
    control_by_id = {control.get("id"): control for control in controls}
    if control_by_id.get("raw_standing_wave_no_truncation_refusal", {}).get("status") != "refused":
        _fail("raw standing wave no-truncation control must refuse")
    zero = control_by_id.get("string_zero_tension_stiffness_refusal")
    if zero is None or zero.get("status") != "refused":
        _fail("zero-tension string control must refuse")
    zero_eigs = [float(item) for item in zero.get("stiffness_eigenvalues", [])]
    if len(zero_eigs) != 3 or any(value != 0.0 for value in zero_eigs):
        _fail("zero-tension string control should have three zero stiffness eigenvalues")


def validate_calibration_analysis_proofs(proofs: dict[str, Any]) -> None:
    if proofs.get("schema_version") != "calibration_analysis_function_proofs.v0":
        _fail("calibration_analysis_function_proofs.json schema_version mismatch")
    proof_items = proofs.get("proofs")
    if not isinstance(proof_items, list) or len(proof_items) != 1:
        _fail("calibration analysis proofs must contain exactly one proof")
    proof = proof_items[0]
    if proof.get("id") != "load_cell_linear_analysis_function":
        _fail("calibration analysis proof id mismatch")
    if proof.get("route") != "calibration_analysis_function":
        _fail("load_cell_linear_analysis_function route mismatch")
    if proof.get("status") != "admissible_analysis_function_record":
        _fail("load_cell_linear_analysis_function status mismatch")
    if proof.get("fit", {}).get("residual_gate", {}).get("passed_simple_gate") is not True:
        _fail("load_cell_linear_analysis_function residual gate did not pass")
    readout = proof.get("analysis_readout", {})
    if float(readout.get("total_local_variance", 0.0)) <= 0.0:
        _fail("load_cell_linear_analysis_function total_local_variance must be positive")
    if float(readout.get("local_precision", 0.0)) <= 0.0:
        _fail("load_cell_linear_analysis_function local_precision must be positive")
    validity = proof.get("validity_range", {})
    new_indication = float(readout.get("new_indication", 0.0))
    if not (float(validity.get("indication_min", 1.0)) <= new_indication <= float(validity.get("indication_max", -1.0))):
        _fail("load_cell_linear_analysis_function new indication must be inside validity range")
    if "hidden_load without a ceiling" not in proof.get("refused_routes", []):
        _fail("load_cell_linear_analysis_function must refuse hidden_load without a ceiling")

    controls = proofs.get("negative_controls")
    if not isinstance(controls, list) or len(controls) != 2:
        _fail("calibration analysis proofs must contain exactly two negative controls")
    control_by_id = {control.get("id"): control for control in controls}
    extrapolation = control_by_id.get("load_cell_analysis_extrapolation_refusal")
    if extrapolation is None or extrapolation.get("status") != "refused":
        _fail("load-cell analysis extrapolation control must refuse")
    span = [float(item) for item in extrapolation.get("calibrated_indication_span", [])]
    requested = float(extrapolation.get("requested_indication", 0.0))
    if len(span) != 2 or span[0] <= requested <= span[1]:
        _fail("load-cell analysis extrapolation control should request outside span")
    rank = control_by_id.get("load_cell_analysis_rank_deficient_refusal")
    if rank is None or rank.get("status") != "refused":
        _fail("load-cell analysis rank-deficient control must refuse")
    rank_eigs = [float(item) for item in rank.get("fisher_eigenvalues", [])]
    if len(rank_eigs) != 2 or min(abs(value) for value in rank_eigs) > 1.0e-8:
        _fail("load-cell analysis rank-deficient control should have a near-zero eigenvalue")


def validate_process_control_gate_proofs(proofs: dict[str, Any]) -> None:
    if proofs.get("schema_version") != "process_control_gate_proofs.v0":
        _fail("process_control_gate_proofs.json schema_version mismatch")
    proof_items = proofs.get("proofs")
    if not isinstance(proof_items, list) or len(proof_items) != 1:
        _fail("process-control gate proofs must contain exactly one proof")
    proof = proof_items[0]
    if proof.get("id") != "check_standard_process_gate":
        _fail("process-control gate proof id mismatch")
    if proof.get("route") != "process_control_gate":
        _fail("check_standard_process_gate route mismatch")
    if proof.get("status") != "gate_passed_no_module_call":
        _fail("check_standard_process_gate status mismatch")
    if proof.get("allowed_route") != "gate_or_refusal_only":
        _fail("check_standard_process_gate allowed route mismatch")
    if proof.get("theorem_local_calls_licensed") is not False:
        _fail("process-control gate must not license theorem-local calls")
    readout = proof.get("gate_readout", {})
    if readout.get("passed_action_limit") is not True:
        _fail("check_standard_process_gate should pass the declared action limit")
    if float(readout.get("max_abs_z", 99.0)) > float(readout.get("action_limit_z", 0.0)):
        _fail("check_standard_process_gate max_abs_z should be inside the action limit")
    for refused in ("visible_precision", "hidden_load", "Fisher precision", "theorem-local kernel call requested"):
        if refused not in proof.get("refused_routes", []):
            _fail(f"check_standard_process_gate missing refusal {refused}")

    controls = proofs.get("negative_controls")
    if not isinstance(controls, list) or len(controls) != 2:
        _fail("process-control gate proofs must contain exactly two negative controls")
    control_by_id = {control.get("id"): control for control in controls}
    drift = control_by_id.get("check_standard_drift_refusal")
    if drift is None or drift.get("status") != "refused":
        _fail("check_standard_drift_refusal missing or not refused")
    drift_readout = drift.get("gate_readout", {})
    if drift_readout.get("passed_action_limit") is not False:
        _fail("check_standard_drift_refusal should fail the declared action limit")
    if float(drift_readout.get("max_abs_z", 0.0)) <= float(drift_readout.get("action_limit_z", 99.0)):
        _fail("check_standard_drift_refusal max_abs_z should exceed the action limit")
    if drift.get("theorem_local_calls_licensed") is not False:
        _fail("check_standard_drift_refusal must not license theorem-local calls")

    missing = control_by_id.get("check_standard_missing_baseline_refusal")
    if missing is None or missing.get("status") != "refused":
        _fail("check_standard_missing_baseline_refusal missing or not refused")
    required_missing = {"historical_baseline_mean", "repeatability_sigma", "bias_or_drift_rule"}
    if set(missing.get("missing_declarations", [])) != required_missing:
        _fail("check_standard_missing_baseline_refusal missing declarations mismatch")
    if missing.get("theorem_local_calls_licensed") is not False:
        _fail("check_standard_missing_baseline_refusal must not license theorem-local calls")


def validate_admissibility_contract(contract: dict[str, Any], root: Path) -> None:
    if contract.get("schema_version") != "admissibility_contract.v0":
        _fail("admissibility_contract_v0.json schema_version mismatch")
    routes = contract.get("routes")
    if not isinstance(routes, list) or len(routes) != 8:
        _fail("admissibility contract must contain exactly eight routes")
    expected_routes = {
        "exact_quadratic_storage",
        "local_quadratic_storage",
        "declared_ceiling_hidden_load",
        "measurement_equation_local_linearisation",
        "observation_equation_fisher",
        "finite_modal_eigenbasis",
        "calibration_analysis_function",
        "process_control_gate",
    }
    actual_routes = {route.get("route_id") for route in routes}
    if actual_routes != expected_routes:
        _fail(f"admissibility contract route ids mismatch: {actual_routes}")
    for route in routes:
        route_id = route["route_id"]
        for field in ("maturity", "covers", "emitted_object", "allowed_surfaces", "required_declarations", "mandatory_refusals", "validated_by"):
            if field not in route:
                _fail(f"{route_id}: missing contract field {field}")
        for list_field in ("covers", "allowed_surfaces", "required_declarations", "mandatory_refusals", "validated_by"):
            if not isinstance(route[list_field], list) or not route[list_field]:
                _fail(f"{route_id}: contract field {list_field} must be non-empty")
        if route_id != "process_control_gate" and not any("hidden_load" in item for item in route["mandatory_refusals"]):
            _fail(f"{route_id}: must include a hidden_load refusal")
        for ref in route["validated_by"]:
            if "*" in ref:
                continue
            candidate = root / ref
            if not candidate.exists():
                _fail(f"{route_id}: validated_by path does not exist: {ref}")
    if "promote" not in str(contract.get("promotion_rule", "")).lower():
        _fail("admissibility contract promotion_rule must be present")


def validate_empirical_dataset_templates(root: Path) -> int:
    protocol = root / "empirical_dataset_protocol_v0.md"
    if not protocol.exists():
        _fail("empirical_dataset_protocol_v0.md missing")
    protocol_text = protocol.read_text(encoding="utf-8")
    for required in ("rc_decay_trace", "beer_lambert_calibration", "hidden-load request"):
        if required not in protocol_text:
            _fail(f"empirical_dataset_protocol_v0.md missing {required}")

    templates_dir = root / "dataset_templates"
    if not templates_dir.exists():
        _fail("dataset_templates directory missing")
    readme = templates_dir / "README.md"
    if not readme.exists():
        _fail("dataset_templates/README.md missing")
    incoming_readme = root / "incoming_data" / "README.md"
    if not incoming_readme.exists():
        _fail("incoming_data/README.md missing")

    expected = {
        "rc_decay_trace": {
            "metadata": templates_dir / "rc_decay_metadata_template.json",
            "csv": templates_dir / "rc_decay_trace_template.csv",
            "columns": ["time_s", "voltage_V"],
            "route": "observation_equation_fisher",
        },
        "beer_lambert_calibration": {
            "metadata": templates_dir / "beer_lambert_metadata_template.json",
            "csv": templates_dir / "beer_lambert_calibration_template.csv",
            "columns": ["concentration", "absorbance"],
            "route": "measurement_equation_local_linearisation_or_observation_equation_fisher",
        },
    }
    for dataset_kind, spec in expected.items():
        metadata = _load_json(spec["metadata"])
        if metadata.get("schema_version") != "empirical_dataset_metadata.v0":
            _fail(f"{dataset_kind}: metadata schema_version mismatch")
        if metadata.get("dataset_kind") != dataset_kind:
            _fail(f"{dataset_kind}: dataset_kind mismatch")
        if metadata.get("route") != spec["route"]:
            _fail(f"{dataset_kind}: route mismatch")
        required_columns = [item.get("name") for item in metadata.get("required_columns", [])]
        if required_columns != spec["columns"]:
            _fail(f"{dataset_kind}: metadata required_columns mismatch")
        csv_header = spec["csv"].read_text(encoding="utf-8").splitlines()[0].split(",")
        if csv_header != spec["columns"]:
            _fail(f"{dataset_kind}: CSV template header mismatch")
        for field in ("measurement_model", "controls", "admission_checks", "allowed_output", "mandatory_refusals"):
            if field not in metadata:
                _fail(f"{dataset_kind}: metadata missing {field}")
        if "hidden_load" not in metadata["mandatory_refusals"]:
            _fail(f"{dataset_kind}: hidden_load must be a mandatory refusal")
    return len(expected)


def validate_online_empirical_probes(root: Path) -> int:
    log_path = root / "online_dataset_acquisition_log.md"
    if not log_path.exists():
        _fail("online_dataset_acquisition_log.md missing")
    log_text = log_path.read_text(encoding="utf-8")
    for required in ("edinburgh_data_driven_chemistry_unit09_section3", "beer_lambert_uvvis_probe", "rc_decay_q7_probe", "Q7 capacitor voltage decay"):
        if required not in log_text:
            _fail(f"online_dataset_acquisition_log.md missing {required}")

    dataset_dir = root / "datasets" / "edinburgh_data_driven_chemistry_unit09_section3"
    if not dataset_dir.exists():
        _fail("Edinburgh UV-Vis dataset directory missing")
    metadata = _load_json(dataset_dir / "metadata.json")
    if metadata.get("schema_version") != "empirical_dataset_metadata.v0":
        _fail("Edinburgh UV-Vis metadata schema_version mismatch")
    if metadata.get("dataset_kind") != "beer_lambert_uvvis_rhodamine_6g":
        _fail("Edinburgh UV-Vis metadata dataset_kind mismatch")
    spectra = sorted(dataset_dir.glob("Rhodamine_6G_in_methanol_*.csv"))
    if len(spectra) != 15:
        _fail(f"Edinburgh UV-Vis dataset should contain 15 Rhodamine spectra, found {len(spectra)}")
    for required_file in ("concentrations.csv", "methanol_baseline.csv"):
        if not (dataset_dir / required_file).exists():
            _fail(f"Edinburgh UV-Vis dataset missing {required_file}")

    probe = _load_json(root / "outputs" / "empirical_datasets" / "beer_lambert_uvvis_probe.json")
    if probe.get("schema_version") != "beer_lambert_uvvis_probe.v0":
        _fail("Beer-Lambert UV-Vis probe schema_version mismatch")
    if probe.get("status") != "admissible_empirical_probe":
        _fail("Beer-Lambert UV-Vis probe should be admissible")
    if float(probe.get("lambda_star_nm", 0.0)) <= 0.0:
        _fail("Beer-Lambert UV-Vis probe lambda_star_nm must be positive")
    calibration = probe.get("calibration", {})
    if int(calibration.get("known_observation_count", 0)) != 12:
        _fail("Beer-Lambert UV-Vis probe should have 12 known calibration observations")
    if float(calibration.get("sigma_absorbance_residual", 0.0)) <= 0.0:
        _fail("Beer-Lambert UV-Vis probe residual sigma must be positive")
    if float(calibration.get("beta_absorbance_per_mM", 0.0)) <= 0.0:
        _fail("Beer-Lambert UV-Vis probe calibration slope must be positive")
    estimates = probe.get("unknown_sample_estimates")
    if not isinstance(estimates, list) or len(estimates) != 3:
        _fail("Beer-Lambert UV-Vis probe should have three unknown-sample estimates")
    if probe.get("unknown_sample_summary", {}).get("inside_known_calibration_span") is not True:
        _fail("Beer-Lambert UV-Vis unknown estimate should be inside calibration span")
    refusals = set(probe.get("mandatory_refusals", []))
    for required_refusal in ("static storage Hessian", "hidden_load", "thermodynamic concentration model"):
        if required_refusal not in refusals:
            _fail(f"Beer-Lambert UV-Vis probe missing refusal {required_refusal}")

    q7_metadata = _load_json(root / "datasets" / "Q7_CapacitorVoltageDecay_metadata.json")
    if q7_metadata.get("schema_version") != "empirical_dataset_metadata.v0":
        _fail("Q7 RC metadata schema_version mismatch")
    if q7_metadata.get("dataset_kind") != "rc_decay_trace":
        _fail("Q7 RC metadata dataset_kind mismatch")
    q7_probe = _load_json(root / "outputs" / "empirical_datasets" / "rc_decay_q7_probe.json")
    if q7_probe.get("schema_version") != "rc_decay_q7_probe.v0":
        _fail("Q7 RC probe schema_version mismatch")
    if q7_probe.get("status") != "admissible_empirical_probe":
        _fail("Q7 RC probe should be admissible")
    if int(q7_probe.get("observation_count", 0)) != 9:
        _fail("Q7 RC probe should contain nine observations")
    fixed_fit = q7_probe.get("fixed_V0_voltage_fit", {})
    if float(fixed_fit.get("tau_hat_s", 0.0)) <= 0.0:
        _fail("Q7 RC tau_hat_s must be positive")
    if float(fixed_fit.get("sigma_voltage_residual_V", 0.0)) <= 0.0:
        _fail("Q7 RC residual sigma must be positive")
    if float(fixed_fit.get("fisher_tau", 0.0)) <= 0.0:
        _fail("Q7 RC fisher_tau must be positive")
    if len(q7_probe.get("data_warnings", [])) != 1:
        _fail("Q7 RC probe should record one nonmonotone data warning")
    q7_refusals = set(q7_probe.get("mandatory_refusals", []))
    for required_refusal in ("static storage Hessian", "hidden_load", "capacitance inference without separately declared resistance"):
        if required_refusal not in q7_refusals:
            _fail(f"Q7 RC probe missing refusal {required_refusal}")
    return 2


def validate_declared_local_positive_contact(root: Path) -> dict[str, int]:
    md_path = root / "declared_local_positive_contact_v0.md"
    json_path = root / "declared_local_positive_contact_v0.json"
    if not md_path.exists():
        _fail("declared_local_positive_contact_v0.md missing")
    if not json_path.exists():
        _fail("declared_local_positive_contact_v0.json missing")

    text = md_path.read_text(encoding="utf-8")
    for required in (
        "What finite declared positive object can the module touch here?",
        "The Gold Quartet",
        "Refusal Pairs",
        "The same local-positive interface can be reached through different physical routes",
    ):
        if required not in text:
            _fail(f"declared_local_positive_contact_v0.md missing {required}")

    synthesis = _load_json(json_path)
    if synthesis.get("schema_version") != "declared_local_positive_contact.v0":
        _fail("declared_local_positive_contact_v0.json schema_version mismatch")
    if synthesis.get("status") != "scientific_synthesis":
        _fail("declared_local_positive_contact_v0.json status mismatch")

    pattern = synthesis.get("core_pattern")
    if not isinstance(pattern, list) or len(pattern) != 5:
        _fail("declared local positive contact core_pattern must have five stages")

    expected_routes = {
        "storage_curvature",
        "local_regime_curvature",
        "metric_curvature",
        "finite_modal_curvature",
        "modal_observer_curvature",
        "reference_ceiling_hidden_load",
        "information_curvature",
        "calibration_analysis_curvature",
        "process_control_gate",
        "thermodynamic_fluctuation_candidate",
    }
    routes = synthesis.get("routes")
    if not isinstance(routes, list):
        _fail("declared local positive contact routes must be a list")
    actual_routes = {route.get("route_id") for route in routes}
    if actual_routes != expected_routes:
        _fail(f"declared local positive contact route ids mismatch: {actual_routes}")
    for route in routes:
        route_id = route["route_id"]
        for field in (
            "curvature_provenance",
            "declared_object",
            "licensed_surface",
            "representative_successes",
            "neighbouring_refusals",
        ):
            if field not in route:
                _fail(f"{route_id}: missing declared contact route field {field}")
        if not isinstance(route["licensed_surface"], list) or not route["licensed_surface"]:
            _fail(f"{route_id}: licensed_surface must be non-empty")
        if route_id != "thermodynamic_fluctuation_candidate" and not isinstance(route["representative_successes"], list):
            _fail(f"{route_id}: representative_successes must be a list")
        if not isinstance(route["neighbouring_refusals"], list) or not route["neighbouring_refusals"]:
            _fail(f"{route_id}: neighbouring_refusals must be non-empty")

    expected_gold = {
        "cross_substrate_hidden_renormalisation",
        "finite_modal_sensor_reference",
        "rc_decay_information_curvature",
        "beer_lambert_sensor_information_curvature",
    }
    quartet = synthesis.get("gold_quartet")
    if not isinstance(quartet, list):
        _fail("declared local positive contact gold_quartet must be a list")
    actual_gold = {item.get("id") for item in quartet}
    if actual_gold != expected_gold:
        _fail(f"declared local positive contact gold quartet ids mismatch: {actual_gold}")
    for item in quartet:
        item_id = item["id"]
        if not isinstance(item.get("artifacts"), list) or not item["artifacts"]:
            _fail(f"{item_id}: gold quartet artifacts must be non-empty")
        if not isinstance(item.get("claim"), str) or not item["claim"]:
            _fail(f"{item_id}: gold quartet claim must be present")

    expected_uncertainty = {
        "instrument_declared",
        "replicate_estimated",
        "residual_estimated",
        "assumed_toy",
    }
    uncertainty = synthesis.get("uncertainty_provenance")
    if not isinstance(uncertainty, list):
        _fail("declared local positive contact uncertainty_provenance must be a list")
    actual_uncertainty = {item.get("id") for item in uncertainty}
    if actual_uncertainty != expected_uncertainty:
        _fail(f"declared local positive contact uncertainty ids mismatch: {actual_uncertainty}")

    synthesis_text = json.dumps(synthesis, sort_keys=True)
    for artifact in (
        "mass_spring_coupled_observe_one",
        "coupled_lc_resonators_observe_one_node",
        "string_modal_hidden_load_with_declared_ceiling",
        "rc_decay_q7_probe",
        "beer_lambert_uvvis_probe",
        "ideal_gas_law_boundary",
    ):
        if artifact not in synthesis_text:
            _fail(f"declared_local_positive_contact_v0.json missing artifact {artifact}")

    return {"routes": len(routes), "gold_quartet": len(quartet), "uncertainty_provenance": len(uncertainty)}


def validate_gold_quartet_demonstration(root: Path) -> int:
    md_path = root / "gold_quartet_demonstration_v0.md"
    json_path = root / "gold_quartet_demonstration_v0.json"
    if not md_path.exists():
        _fail("gold_quartet_demonstration_v0.md missing")
    if not json_path.exists():
        _fail("gold_quartet_demonstration_v0.json missing")

    text = md_path.read_text(encoding="utf-8")
    for required in (
        "Cross-Substrate Hidden Renormalisation",
        "Finite Modal Sensor Reference",
        "RC Decay Information Curvature",
        "Beer-Lambert Sensor Information Curvature",
        "The refusals are part of the result.",
    ):
        if required not in text:
            _fail(f"gold_quartet_demonstration_v0.md missing {required}")

    demo = _load_json(json_path)
    if demo.get("schema_version") != "gold_quartet_demonstration.v0":
        _fail("gold_quartet_demonstration_v0.json schema_version mismatch")
    if demo.get("status") != "external_facing_demonstration":
        _fail("gold_quartet_demonstration_v0.json status mismatch")
    demonstrations = demo.get("demonstrations")
    if not isinstance(demonstrations, list) or len(demonstrations) != 4:
        _fail("gold quartet demonstration must contain exactly four demonstrations")

    expected_ids = {
        "cross_substrate_hidden_renormalisation",
        "finite_modal_sensor_reference",
        "rc_decay_information_curvature",
        "beer_lambert_sensor_information_curvature",
    }
    by_id = {item.get("id"): item for item in demonstrations}
    if set(by_id) != expected_ids:
        _fail(f"gold quartet demonstration ids mismatch: {set(by_id)}")
    for item in demonstrations:
        item_id = item["id"]
        for field in ("route", "source_outputs", "key_readout", "lesson", "mandatory_refusal"):
            if field not in item:
                _fail(f"{item_id}: missing gold quartet field {field}")
        for ref in item["source_outputs"]:
            if not (root / ref).exists():
                _fail(f"{item_id}: source output does not exist: {ref}")

    cross = by_id["cross_substrate_hidden_renormalisation"]["key_readout"]
    spring = _load_json(root / "outputs" / "mass_spring_coupled_observe_one.json")
    lc = _load_json(root / "outputs" / "coupled_lc_resonators_observe_one_node.json")
    _assert_close(spring["readout"]["visible_precision"][0][0], cross["spring_visible_precision"], "spring visible precision")
    _assert_close(spring["readout"]["ceiling_minus_visible_precision"][0][0], cross["spring_ceiling_minus_visible_precision"], "spring ceiling minus visible")
    _assert_close(spring["readout"]["hidden_load"][0][0], cross["spring_hidden_load"], "spring hidden load")
    _assert_close(lc["readout"]["visible_precision"][0][0], cross["lc_visible_precision"], "LC visible precision")
    _assert_close(lc["readout"]["ceiling_minus_visible_precision"][0][0], cross["lc_ceiling_minus_visible_precision"], "LC ceiling minus visible")
    _assert_close(lc["readout"]["hidden_load"][0][0], cross["lc_hidden_load"], "LC hidden load")

    modal = by_id["finite_modal_sensor_reference"]["key_readout"]
    string_fixed = _load_json(root / "outputs" / "string_fixed_end_three_mode.json")
    point_sensor = _load_json(root / "outputs" / "string_point_sensor_modal_observer.json")
    modal_hidden = _load_json(root / "outputs" / "string_modal_hidden_load_with_declared_ceiling.json")
    for index, expected_value in enumerate(modal["three_mode_stiffness_diagonal"]):
        _assert_close(string_fixed["readout"]["visible_precision"][index][index], expected_value, f"string modal stiffness diagonal {index + 1}")
    _assert_close(point_sensor["readout"]["visible_precision"][0][0], modal["point_sensor_visible_precision"], "point sensor visible precision")
    _assert_close(modal_hidden["readout"]["ceiling_minus_visible_precision"][0][0], modal["reference_ceiling_minus_visible_precision"], "modal reference ceiling minus visible")
    _assert_close(modal_hidden["readout"]["hidden_load"][0][0], modal["reference_hidden_load"], "modal reference hidden load")

    rc = by_id["rc_decay_information_curvature"]["key_readout"]
    rc_probe = _load_json(root / "outputs" / "empirical_datasets" / "rc_decay_q7_probe.json")
    fixed_fit = rc_probe["fixed_V0_voltage_fit"]
    _assert_close(fixed_fit["tau_hat_s"], rc["tau_hat_s"], "RC tau_hat_s")
    _assert_close(fixed_fit["fisher_tau"], rc["fisher_tau"], "RC fisher_tau")
    _assert_close(fixed_fit["local_variance_tau_s2"], rc["local_variance_tau_s2"], "RC local_variance_tau_s2")
    if len(rc_probe["data_warnings"]) != int(rc["data_warnings"]):
        _fail("RC data warning count mismatch")

    beer = by_id["beer_lambert_sensor_information_curvature"]["key_readout"]
    beer_probe = _load_json(root / "outputs" / "empirical_datasets" / "beer_lambert_uvvis_probe.json")
    _assert_close(beer_probe["lambda_star_nm"], beer["lambda_star_nm"], "Beer-Lambert lambda_star_nm")
    if int(beer_probe["calibration"]["known_observation_count"]) != int(beer["known_observation_count"]):
        _fail("Beer-Lambert known observation count mismatch")
    _assert_close(beer_probe["calibration"]["beta_absorbance_per_mM"], beer["beta_absorbance_per_mM"], "Beer-Lambert beta")
    _assert_close(beer_probe["calibration"]["sigma_absorbance_residual"], beer["sigma_absorbance_residual"], "Beer-Lambert sigma residual")
    _assert_close(beer_probe["unknown_sample_summary"]["mean_estimated_concentration_mM"], beer["mean_estimated_concentration_mM"], "Beer-Lambert mean concentration")
    if bool(beer_probe["unknown_sample_summary"]["inside_known_calibration_span"]) is not bool(beer["inside_known_calibration_span"]):
        _fail("Beer-Lambert calibration span flag mismatch")

    return len(demonstrations)


def validate_main_module_contact_bridge(root: Path) -> int:
    md_path = root / "main_module_contact_bridge_v0.md"
    json_path = root / "main_module_contact_bridge_v0.json"
    if not md_path.exists():
        _fail("main_module_contact_bridge_v0.md missing")
    if not json_path.exists():
        _fail("main_module_contact_bridge_v0.json missing")

    text = md_path.read_text(encoding="utf-8")
    for required in (
        "Main Module Contact Bridge v0",
        "declared local positive contact",
        "Do not broaden by adding more formula cards during wrap.",
        "record-only Fisher/covariance objects",
    ):
        if required not in text:
            _fail(f"main_module_contact_bridge_v0.md missing {required}")

    bridge = _load_json(json_path)
    if bridge.get("schema_version") != "main_module_contact_bridge.v0":
        _fail("main_module_contact_bridge_v0.json schema_version mismatch")
    if bridge.get("status") != "wrap_ready_synthesis":
        _fail("main_module_contact_bridge_v0.json status mismatch")
    handoffs = bridge.get("module_handoffs")
    if not isinstance(handoffs, list) or len(handoffs) != 5:
        _fail("main module contact bridge must contain exactly five handoffs")
    expected = {
        "visible_precision_direct",
        "hidden_load_with_declared_ceiling",
        "local_positive_record",
        "process_gate",
        "principled_refusal",
    }
    actual = {item.get("handoff_id") for item in handoffs}
    if actual != expected:
        _fail(f"main module contact bridge handoff ids mismatch: {actual}")
    for item in handoffs:
        handoff_id = item["handoff_id"]
        for field in ("module_surface", "input_contract", "current_examples", "licensed_claim", "mandatory_boundary"):
            if field not in item:
                _fail(f"{handoff_id}: missing main-module bridge field {field}")
        if not isinstance(item["current_examples"], list) or not item["current_examples"]:
            _fail(f"{handoff_id}: current_examples must be non-empty")

    bridge_text = json.dumps(bridge, sort_keys=True)
    for required in (
        "visible_precision",
        "hidden_load",
        "T >= Phi",
        "rc_decay_q7_probe",
        "beer_lambert_uvvis_probe",
        "check_standard_process_gate",
        "ideal_gas_law_boundary",
        "instrument-declared uncertainty",
        "thermodynamic positive contact",
    ):
        if required not in bridge_text:
            _fail(f"main_module_contact_bridge_v0.json missing {required}")
    return len(handoffs)


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate operationalise metadata artifacts.")
    parser.add_argument("--root", default=str(ROOT))
    args = parser.parse_args()
    root = Path(args.root)

    backlog = _load_json(root / "candidate_backlog.json")
    diagnostics = _load_json(root / "non_integration_diagnostics.json")
    sweep = _load_json(root / "outputs" / "sweeps" / "coupled_hidden_renormalisation_sweep.json")
    probes = _load_json(root / "outputs" / "failure_probes" / "failure_probe_suite.json")
    proofs = _load_json(root / "outputs" / "local_curvature" / "local_positive_curvature_proofs.json")
    stress = _load_json(root / "outputs" / "local_curvature" / "local_curvature_stress_tests.json")
    measurement_routes = _load_json(root / "measurement_route_backlog.json")
    measurement_equation_proofs = _load_json(root / "outputs" / "measurement_equation" / "measurement_equation_linearisation_proofs.json")
    observation_equation_proofs = _load_json(root / "outputs" / "observation_equation" / "observation_equation_fisher_proofs.json")
    finite_modal_proofs = _load_json(root / "outputs" / "finite_modal" / "finite_modal_eigenbasis_proofs.json")
    calibration_analysis_proofs = _load_json(root / "outputs" / "calibration_analysis" / "calibration_analysis_function_proofs.json")
    process_control_gate_proofs = _load_json(root / "outputs" / "process_control" / "process_control_gate_proofs.json")
    admissibility_contract = _load_json(root / "admissibility_contract_v0.json")

    validate_backlog(backlog)
    validate_diagnostics(backlog, diagnostics)
    validate_sweep(sweep)
    validate_failure_probes(backlog, diagnostics, probes)
    validate_local_curvature_proofs(proofs)
    validate_local_curvature_stress_tests(stress)
    validate_signoff(root / "card_signoff_v0.md")
    validate_measurement_route_backlog(measurement_routes)
    validate_measurement_equation_proofs(measurement_equation_proofs)
    validate_observation_equation_proofs(observation_equation_proofs)
    validate_finite_modal_proofs(finite_modal_proofs)
    validate_calibration_analysis_proofs(calibration_analysis_proofs)
    validate_process_control_gate_proofs(process_control_gate_proofs)
    validate_admissibility_contract(admissibility_contract, root)
    dataset_template_count = validate_empirical_dataset_templates(root)
    online_probe_count = validate_online_empirical_probes(root)
    declared_contact_counts = validate_declared_local_positive_contact(root)
    gold_quartet_demonstrations = validate_gold_quartet_demonstration(root)
    main_module_handoffs = validate_main_module_contact_bridge(root)

    print(
        json.dumps(
            {
                "candidate_backlog_items": len(backlog["items"]),
                "diagnosed_non_carded_items": len(diagnostics["items"]),
                "sweep_rows": len(sweep["rows"]),
                "failure_probes": len(probes["probes"]),
                "local_curvature_proofs": len(proofs["proofs"]),
                "local_curvature_stress_tests": len(stress["tests"]),
                "measurement_route_items": len(measurement_routes["items"]),
                "measurement_equation_proofs": len(measurement_equation_proofs["proofs"]),
                "observation_equation_proofs": len(observation_equation_proofs["proofs"]),
                "observation_equation_negative_controls": len(observation_equation_proofs["negative_controls"]),
                "finite_modal_proofs": len(finite_modal_proofs["proofs"]),
                "finite_modal_negative_controls": len(finite_modal_proofs["negative_controls"]),
                "calibration_analysis_proofs": len(calibration_analysis_proofs["proofs"]),
                "calibration_analysis_negative_controls": len(calibration_analysis_proofs["negative_controls"]),
                "process_control_gate_proofs": len(process_control_gate_proofs["proofs"]),
                "process_control_gate_negative_controls": len(process_control_gate_proofs["negative_controls"]),
                "admissibility_contract_routes": len(admissibility_contract["routes"]),
                "empirical_dataset_templates": dataset_template_count,
                "online_empirical_probes": online_probe_count,
                "declared_contact_routes": declared_contact_counts["routes"],
                "declared_contact_gold_quartet": declared_contact_counts["gold_quartet"],
                "declared_contact_uncertainty_provenance": declared_contact_counts["uncertainty_provenance"],
                "gold_quartet_demonstrations": gold_quartet_demonstrations,
                "main_module_handoffs": main_module_handoffs,
                "status": "passed",
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
