from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PYTHON = Path(r"C:\Users\14566\miniconda3\envs\trkb_openmm\python.exe")

CANDIDATES = [
    ("2.3", "L006"),
    ("6.2", "L008"),
    ("17.2", "L010"),
]
SYSTEMS = ["20chol", "40chol"]


def run_stage(script_name: str, system: str, ligand_id: str, resname: str) -> int:
    env = os.environ.copy()
    env["TRKB_LIGAND_ID"] = ligand_id
    env["TRKB_RESNAME"] = resname
    script = SCRIPT_DIR / script_name
    completed = subprocess.run([str(PYTHON), "-B", str(script), system], env=env, cwd=str(Path(r"C:\TRKB_WP2")))
    return completed.returncode


def load_json(path: Path) -> dict:
    import json

    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def stage_paths(ligand_id: str, system: str) -> dict[str, Path]:
    root = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build") / f"{ligand_id}_{system}"
    return {
        "gate": root / "reports" / "pilot_gate_summary.json",
        "md50": root / "md50ps_R2" / "reports" / f"md50ps_{system}_R2_result.json",
        "md100": root / "md100ps_R2" / "reports" / f"md100ps_{system}_R2_result.json",
    }


def main() -> int:
    for ligand_id, resname in CANDIDATES:
        for system in SYSTEMS:
            print(f"=== {ligand_id} / {resname} / {system}: gate pipeline start ===")
            run_stage("run_nonproduction_pilot_generic.py", system, ligand_id, resname)
            paths = stage_paths(ligand_id, system)
            gate = load_json(paths["gate"])
            if not gate.get("recommend_50ps", False):
                print(f"=== {ligand_id} / {resname} / {system}: stop after 10 ps gate ===")
                subprocess.run([str(PYTHON), "-B", str(SCRIPT_DIR / "generate_batch_status_report.py")], cwd=str(Path(r"C:\TRKB_WP2")))
                continue

            run_stage("run_50ps_generic.py", system, ligand_id, resname)
            md50 = load_json(paths["md50"])
            if not md50.get("recommend_100ps", False):
                print(f"=== {ligand_id} / {resname} / {system}: stop after 50 ps gate ===")
                subprocess.run([str(PYTHON), "-B", str(SCRIPT_DIR / "generate_batch_status_report.py")], cwd=str(Path(r"C:\TRKB_WP2")))
                continue

            run_stage("run_100ps_generic.py", system, ligand_id, resname)
            md100 = load_json(paths["md100"])
            if not md100.get("completed_100ps", False):
                print(f"=== {ligand_id} / {resname} / {system}: stopped during 100 ps ===")
            else:
                print(f"=== {ligand_id} / {resname} / {system}: completed 100 ps ===")
            subprocess.run([str(PYTHON), "-B", str(SCRIPT_DIR / "generate_batch_status_report.py")], cwd=str(Path(r"C:\TRKB_WP2")))

    subprocess.run([str(PYTHON), "-B", str(SCRIPT_DIR / "generate_batch_status_report.py")], cwd=str(Path(r"C:\TRKB_WP2")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
