from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import venv
from dataclasses import asdict, dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tools" / "outputs"


@dataclass(frozen=True)
class CommandResult:
    label: str
    cwd: str
    command: list[str]
    returncode: int
    stdout_tail: str
    stderr_tail: str


@dataclass(frozen=True)
class InstallSoakSummary:
    editable_install: tuple[CommandResult, ...]
    test_runs: tuple[CommandResult, ...]
    example_runs: tuple[CommandResult, ...]
    wheel_builds: tuple[CommandResult, ...]
    wheel_install: tuple[CommandResult, ...]
    installed_surface_smoke: tuple[CommandResult, ...]


def run_install_soak() -> InstallSoakSummary:
    OUT.mkdir(parents=True, exist_ok=True)
    wheel_dir = OUT / "wheel_dist"
    if wheel_dir.exists():
        shutil.rmtree(wheel_dir)
    wheel_dir.mkdir(parents=True, exist_ok=True)

    editable_install = (
        _run("editable_root", [sys.executable, "-m", "pip", "install", "-e", "."], ROOT),
        _run("editable_nomodescent", [sys.executable, "-m", "pip", "install", "-e", "."], ROOT / "nomodescent"),
        _run("editable_evidence", [sys.executable, "-m", "pip", "install", "-e", "."], ROOT / "evidence"),
    )
    _require_success(editable_install)

    test_runs = (
        _run("pytest_root", [sys.executable, "-m", "pytest", "-q"], ROOT),
        _run("pytest_nomodescent", [sys.executable, "-m", "pytest"], ROOT / "nomodescent"),
        _run("pytest_evidence", [sys.executable, "-m", "pytest"], ROOT / "evidence"),
    )
    _require_success(test_runs)

    example_runs = (
        _run("example_entanglement", [sys.executable, "-m", "examples.entanglement_hidden_load.run_all"], ROOT),
        _run("example_bell", [sys.executable, "-m", "examples.bell_common_gluing.run_all"], ROOT),
        _run("example_arrow", [sys.executable, "-m", "examples.arrow_rank_deficiency.run_all"], ROOT),
        _run("descent_bell", [sys.executable, "-m", "worked_examples.bell_descent.run_main"], ROOT / "nomodescent"),
        _run("descent_rg", [sys.executable, "-m", "worked_examples.free_gaussian_rg.run_main"], ROOT / "nomodescent"),
        _run("descent_replication", [sys.executable, "-m", "worked_examples.replication_fragility.run_main"], ROOT / "nomodescent"),
        _run("evidence_bell", [sys.executable, "-m", "worked_examples.bell_evidence_encoding.run_main"], ROOT / "evidence"),
        _run("evidence_replication", [sys.executable, "-m", "worked_examples.replication_protocol_encoding.run_main"], ROOT / "evidence"),
        _run("evidence_benchmark", [sys.executable, "-m", "worked_examples.benchmark_blindness_encoding.run_main"], ROOT / "evidence"),
        _run("micro_real_bell", [sys.executable, "-m", "micro_real_bundles.bell_counts_bundle.run_main"], ROOT / "evidence"),
        _run("micro_real_iris", [sys.executable, "-m", "micro_real_bundles.iris_protocol_mismatch.run_main"], ROOT / "evidence"),
        _run("micro_real_benchmark", [sys.executable, "-m", "micro_real_bundles.leaderboard_benchmark_slice.run_main"], ROOT / "evidence"),
    )
    _require_success(example_runs)

    wheel_builds = (
        _run("wheel_root", [sys.executable, "-m", "pip", "wheel", ".", "-w", str(wheel_dir), "--no-deps"], ROOT),
        _run("wheel_nomodescent", [sys.executable, "-m", "pip", "wheel", ".", "-w", str(wheel_dir), "--no-deps"], ROOT / "nomodescent"),
        _run("wheel_evidence", [sys.executable, "-m", "pip", "wheel", ".", "-w", str(wheel_dir), "--no-deps"], ROOT / "evidence"),
    )
    _require_success(wheel_builds)

    with tempfile.TemporaryDirectory(dir=OUT) as temp_dir:
        env_dir = Path(temp_dir) / "stack_wheel_env"
        builder = venv.EnvBuilder(with_pip=True, clear=True, system_site_packages=True)
        builder.create(env_dir)
        py = env_dir / ("Scripts" if os.name == "nt" else "bin") / ("python.exe" if os.name == "nt" else "python")

        wheel_install = (
            _run("wheel_install_root", [str(py), "-m", "pip", "install", "--no-deps", str(_find_wheel(wheel_dir, "nomogeo"))], ROOT),
            _run("wheel_install_nomodescent", [str(py), "-m", "pip", "install", "--no-deps", str(_find_wheel(wheel_dir, "nomodescent"))], ROOT),
            _run("wheel_install_evidence", [str(py), "-m", "pip", "install", "--no-deps", str(_find_wheel(wheel_dir, "evidence"))], ROOT),
        )
        _require_success(wheel_install)

        installed_surface_smoke = (
            _run(
                "wheel_imports",
                [
                    str(py),
                    "-c",
                    (
                        "import numpy as np; "
                        "from nomogeo import visible_precision; "
                        "from nomodescent import ObserverSpec, ProblemSpec, VisibleEvidenceSpec, GoalSpec, common_descent_test; "
                        "from evidence import EvidenceBundle, ProtocolObservation, ObserverHypothesis, SourceRef, ExtractionRecord, encode_matrix_observation, assemble_problem_spec; "
                        "H=np.array([[2.0,0.3],[0.3,1.5]]); "
                        "C=np.array([[1.0,0.0]]); "
                        "phi=visible_precision(H,C); "
                        "problem=ProblemSpec(name='wheel', latent_dim=2, observers=(ObserverSpec(name='obs', matrix=[[1.0,0.0]]),), evidence=(VisibleEvidenceSpec(name='cov', observer='obs', kind='covariance', matrix=[[1.0]]),), goals=(GoalSpec(kind='common_completion'),)); "
                        "result=common_descent_test(problem); "
                        "bundle=EvidenceBundle(name='bundle', latent_dim=2, protocol_observations=(ProtocolObservation(name='p', facts=('x',), candidate_family=('h',), source_ref=SourceRef(source='wheel'), extraction=ExtractionRecord(extraction_mode='exact_extraction', epistemic_status='exact', authoritative=True)),), observer_hypotheses=(ObserverHypothesis(name='h', matrix=[[1.0,0.0]], protocol_name='p', features=('x',), source_ref=SourceRef(source='wheel'), extraction=ExtractionRecord(extraction_mode='exact_extraction', epistemic_status='exact', authoritative=True)),), matrix_observations=(encode_matrix_observation(name='cov_p', matrix_role='visible_object', matrix_kind='covariance', matrix=((1.0,),), observer_name='p', source='wheel'),)); "
                        "assembled=assemble_problem_spec(bundle, observer_selection={'p':'h'}); "
                        "print(round(float(phi[0,0]), 6), result.classification, assembled.classification)"
                    ),
                ],
                ROOT,
            ),
        )
        _require_success(installed_surface_smoke)

        summary = InstallSoakSummary(
            editable_install=editable_install,
            test_runs=test_runs,
            example_runs=example_runs,
            wheel_builds=wheel_builds,
            wheel_install=wheel_install,
            installed_surface_smoke=installed_surface_smoke,
        )

    return summary


def write_install_soak_outputs(summary: InstallSoakSummary) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "install_soak_summary.json").write_text(json.dumps(asdict(summary), indent=2), encoding="utf-8")


def main() -> None:
    summary = run_install_soak()
    write_install_soak_outputs(summary)
    print(json.dumps(asdict(summary), indent=2))


def _run(label: str, command: list[str], cwd: Path) -> CommandResult:
    temp_root = Path(tempfile.mkdtemp(dir=OUT, prefix="subprocess_tmp_"))
    env = os.environ.copy()
    env["TMP"] = str(temp_root)
    env["TEMP"] = str(temp_root)
    env.pop("PIP_BUILD_TRACKER", None)
    try:
        completed = subprocess.run(command, cwd=cwd, capture_output=True, text=True, check=False, env=env)
        return CommandResult(
            label=label,
            cwd=str(cwd),
            command=command,
            returncode=completed.returncode,
            stdout_tail=_tail(completed.stdout),
            stderr_tail=_tail(completed.stderr),
        )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def _tail(text: str, limit: int = 2000) -> str:
    if len(text) <= limit:
        return text
    return text[-limit:]


def _require_success(results: tuple[CommandResult, ...]) -> None:
    failures = [result for result in results if result.returncode != 0]
    if failures:
        labels = ", ".join(result.label for result in failures)
        raise RuntimeError(f"install soak commands failed: {labels}")


def _find_wheel(directory: Path, prefix: str) -> Path:
    matches = sorted(directory.glob(f"{prefix}-*.whl"))
    if not matches:
        raise FileNotFoundError(f"no wheel found for prefix '{prefix}' in {directory}")
    return matches[-1]


if __name__ == "__main__":
    main()

