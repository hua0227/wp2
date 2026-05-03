from __future__ import annotations

import csv
from pathlib import Path

from common_40chol import (
    BOUND_PDB,
    BOUND_PSF,
    CLEANED_PDB,
    CLEANED_PSF,
    LOGS_DIR,
    MINIMIZED_PDB,
    REPORTS_DIR,
)


SUMMARY_MD = REPORTS_DIR / "pilot_40chol_summary.md"
CLASH_CSV = REPORTS_DIR / "clash_diagnosis_40chol.csv"
CLASH_MD = REPORTS_DIR / "clash_diagnosis_40chol.md"
MIN_CSV = REPORTS_DIR / "R2_minimization_40chol_summary.csv"
PROBE_CSV = REPORTS_DIR / "very_short_probe_40chol.csv"
READTEST_LOG = LOGS_DIR / "openmm_readtest_40chol.log"
PSFGEN_LOG = LOGS_DIR / "psfgen_build_40chol.log"
MIN_LOG = LOGS_DIR / "R2_pbcaware_minimization_40chol.log"
PROBE_LOG = LOGS_DIR / "very_short_probe_40chol.log"


def read_first_csv(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    return rows[0] if rows else {}


def read_probe_status(path: Path) -> tuple[bool, str]:
    if not path.exists():
        return False, "not executed"
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return True, "no rows"
    last = rows[-1]
    passed = last.get("step") == "4000" and last.get("status") == "ok"
    return passed, f"last step {last.get('step')} status {last.get('status')}"


def find_initial_energy(log_path: Path) -> str:
    if not log_path.exists():
        return "not available"
    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "Initial potential energy (kJ/mol):" in line:
            return line.split(":", 1)[1].strip()
    return "not available"


def count_clashes(path: Path) -> tuple[int, int]:
    if not path.exists():
        return 0, 0
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    severe = sum(1 for row in rows if row.get("contact_type") == "severe_clash")
    close = sum(1 for row in rows if row.get("contact_type") in {"severe_clash", "close_contact"})
    return severe, close


def cleanup_summary() -> str:
    if not CLASH_MD.exists():
        return "not available"
    lines = CLASH_MD.read_text(encoding="utf-8", errors="replace").splitlines()
    decision = [line for line in lines if line.startswith("- Decision:")]
    selected = [line.strip() for line in lines if line.startswith("  - ")]
    if selected:
        return f"{decision[0].replace('- Decision: ', '')}; deleted {', '.join(item[2:] for item in selected)}"
    if decision:
        return decision[0].replace("- Decision: ", "")
    return "not available"


def main() -> int:
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    initial_energy = find_initial_energy(READTEST_LOG)
    severe, close = count_clashes(CLASH_CSV)
    min_row = read_first_csv(MIN_CSV)
    probe_passed, probe_status = read_probe_status(PROBE_CSV)
    min_success = bool(min_row) and min_row.get("energy_finite") == "1" and min_row.get("coordinate_nan") == "0" and min_row.get("force_nan") == "0"
    gate_passed = min_row.get("probe_gate_passed") == "1"
    recommend_10ps = min_success and gate_passed and (not PROBE_CSV.exists() or probe_passed)
    lines = [
        "# 14.2 / L002 / 40chol Non-Production Pilot Summary",
        "",
        "## Purpose",
        "",
        "This stage tests whether the ligand-bound 40% CHOL system can be built, read by OpenMM, diagnosed for clashes, and minimized with R2 PBC-aware restraints before any longer restrained MD.",
        "",
        "## Relationship To 20chol Pilot",
        "",
        "The 20chol pilot is used only as a workflow reference. No 20chol coordinates, cleaned systems, or extended 20chol MD products were copied into this 40chol build.",
        "",
        "## Build And Read Test",
        "",
        f"- Ligand-bound PSF generated: {BOUND_PSF.exists()}",
        f"- Ligand-bound PDB generated: {BOUND_PDB.exists()}",
        f"- psfgen log: `{PSFGEN_LOG}`",
        f"- OpenMM read test log: `{READTEST_LOG}`",
        f"- OpenMM initial energy: {initial_energy} kJ/mol",
        "",
        "## Clash Diagnosis And Cleanup",
        "",
        f"- Severe ligand/non-ligand heavy-atom clashes < 1.0 A: {severe}",
        f"- Close ligand/non-ligand heavy-atom contacts < 1.5 A: {close}",
        f"- Targeted cleanup: {cleanup_summary()}",
        f"- Cleaned PSF: `{CLEANED_PSF}`",
        f"- Cleaned PDB: `{CLEANED_PDB}`",
        "",
        "## R2 PBC-Aware Minimization",
        "",
        f"- R2 minimization success: {min_success}",
        f"- Minimized PDB: `{MINIMIZED_PDB}`",
        f"- Final energy: {min_row.get('final_energy_kJ_mol', 'not available')} kJ/mol",
        f"- Max force: {min_row.get('max_force_kJ_mol_nm', 'not available')} kJ/mol/nm",
        f"- Top force atom: {min_row.get('top_force_atom', 'not available')}",
        f"- Ligand min distance to TYR13: {min_row.get('ligand_min_distance_TYR13_A', 'not available')} A",
        f"- Ligand min distance to VAL17: {min_row.get('ligand_min_distance_VAL17_A', 'not available')} A",
        f"- Ligand min distance to SER20: {min_row.get('ligand_min_distance_SER20_A', 'not available')} A",
        f"- Ligand near chains count: {min_row.get('ligand_near_chains_count', 'not available')}",
        f"- Severe/close after minimization: {min_row.get('severe_lt_1A', 'not available')}/{min_row.get('close_lt_1p5A', 'not available')}",
        "",
        "## Very-Short Probe",
        "",
        f"- Very-short 0.2 ps 50 K probe executed: {PROBE_CSV.exists()}",
        f"- Very-short probe passed: {probe_passed}",
        f"- Very-short probe status: {probe_status}",
        f"- Very-short probe log: `{PROBE_LOG}`",
        "",
        "## Recommendation",
        "",
        f"- Recommend entering 10 ps 40chol restrained MD: {recommend_10ps}",
        "- This is not production MD.",
    ]
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"pilot summary: {SUMMARY_MD}")
    print(f"recommend_10ps_40chol_restrained_MD: {recommend_10ps}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
