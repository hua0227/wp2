from __future__ import annotations

import csv
import json
import math
from pathlib import Path


BASE_DIR = Path(r"C:\TRKB_WP2\ligand_bound_MD")
PILOT_BUILD_DIR = BASE_DIR / "pilot_build"
LIGAND_ID = "9.2"
RESNAME = "L004"
SYSTEMS = ["20chol", "40chol"]
SUMMARY_MD = PILOT_BUILD_DIR / "9.2_gate_summary.md"
SUMMARY_CSV = PILOT_BUILD_DIR / "9.2_gate_summary.csv"

FIELDS = [
    "candidate_id",
    "resname",
    "system",
    "PSF_PDB_generated",
    "read_test_succeeded",
    "initial_energy_kJ_mol",
    "severe_before_cleanup",
    "close_before_cleanup",
    "cleanup_executed",
    "deleted_molecules",
    "minimization_final_energy_kJ_mol",
    "very_short_probe_passed",
    "completed_10ps",
    "completed_50ps",
    "completed_100ps",
    "final_temperature_K",
    "final_potential_energy_kJ_mol",
    "final_TYR13_distance_A",
    "final_VAL17_distance_A",
    "final_SER20_distance_A",
    "near_chains_count",
    "ligand_rmsd_final_A",
    "ligand_rmsd_max_A",
    "final_severe_contacts",
    "final_close_contacts",
    "stage_reached",
    "stop_reason",
    "recommend_comparison_analysis",
]


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def to_float(value: object) -> float:
    if value is None:
        return math.nan
    text = str(value).strip()
    if text == "":
        return math.nan
    try:
        return float(text)
    except ValueError:
        return math.nan


def fmt(value: object, digits: int = 3) -> str:
    number = to_float(value)
    if math.isfinite(number):
        return f"{number:.{digits}f}"
    if isinstance(value, bool):
        return "YES" if value else "NO"
    return str(value)


def stage_reached(row: dict[str, object]) -> str:
    if row.get("completed_100ps"):
        return "100ps"
    if row.get("completed_50ps"):
        return "50ps"
    if row.get("completed_10ps"):
        return "10ps"
    if row.get("very_short_probe_passed"):
        return "very_short_probe"
    if row.get("minimization_succeeded"):
        return "minimization"
    if row.get("read_test_succeeded"):
        return "openmm_read_test"
    if row.get("PSF_PDB_generated"):
        return "psfgen_build"
    if row.get("preflight_ok"):
        return "preflight"
    return "not_started_or_failed_preflight"


def row_for_system(system: str) -> dict[str, object]:
    root = PILOT_BUILD_DIR / f"{LIGAND_ID}_{system}"
    gate = read_json(root / "reports" / "pilot_gate_summary.json")
    md50 = read_json(root / "md50ps_R2" / "reports" / f"md50ps_{system}_R2_result.json")
    md100 = read_json(root / "md100ps_R2" / "reports" / f"md100ps_{system}_R2_result.json")

    final_row_100 = md100.get("final_row", {}) if isinstance(md100.get("final_row", {}), dict) else {}
    final_row_50 = md50.get("final_row", {}) if isinstance(md50.get("final_row", {}), dict) else {}
    final_source = final_row_100 or final_row_50 or {}

    completed_100 = bool(md100.get("completed_100ps"))
    completed_50 = bool(md50.get("completed_50ps"))
    completed_10 = bool(gate.get("completed_10ps"))
    recommend_comparison = bool(md100.get("recommend_comparison_analysis")) if md100 else False
    row = {
        "candidate_id": LIGAND_ID,
        "resname": RESNAME,
        "system": system,
        "preflight_ok": bool(gate.get("preflight_ok")),
        "PSF_PDB_generated": bool(gate.get("psf_pdb_generated")),
        "read_test_succeeded": bool(gate.get("read_test_succeeded")),
        "initial_energy_kJ_mol": gate.get("initial_energy_kJ_mol", float("nan")),
        "severe_before_cleanup": gate.get("severe_before_cleanup", ""),
        "close_before_cleanup": gate.get("close_before_cleanup", ""),
        "cleanup_executed": bool(gate.get("cleanup_executed")),
        "deleted_molecules": gate.get("deleted_molecules", []),
        "minimization_final_energy_kJ_mol": gate.get("minimization_final_energy_kJ_mol", float("nan")),
        "very_short_probe_passed": bool(gate.get("very_short_probe_passed")),
        "completed_10ps": completed_10,
        "completed_50ps": completed_50,
        "completed_100ps": completed_100,
        "final_temperature_K": final_source.get("temperature_K", gate.get("final_temperature_K", float("nan"))),
        "final_potential_energy_kJ_mol": final_source.get("potential_energy_kJ_mol", gate.get("final_potential_energy_kJ_mol", float("nan"))),
        "final_TYR13_distance_A": final_source.get("ligand_distance_TYR13_A", gate.get("ligand_distance_TYR13_A", float("nan"))),
        "final_VAL17_distance_A": final_source.get("ligand_distance_VAL17_A", gate.get("ligand_distance_VAL17_A", float("nan"))),
        "final_SER20_distance_A": final_source.get("ligand_distance_SER20_A", gate.get("ligand_distance_SER20_A", float("nan"))),
        "near_chains_count": final_source.get("ligand_near_chains_count", gate.get("ligand_near_chains_count", "")),
        "ligand_rmsd_final_A": md100.get("ligand_rmsd_final_A", md50.get("ligand_rmsd_final_A", gate.get("ligand_rmsd_A", float("nan")))),
        "ligand_rmsd_max_A": md100.get("ligand_rmsd_max_A", md50.get("ligand_rmsd_max_A", gate.get("ligand_rmsd_A", float("nan")))),
        "final_severe_contacts": final_source.get("severe_lt_1A", gate.get("severe_lt_1A", "")),
        "final_close_contacts": final_source.get("close_lt_1p5A", gate.get("close_lt_1p5A", "")),
        "stop_reason": md100.get("stop_reason") or md50.get("stop_reason") or gate.get("stop_reason") or gate.get("error") or "missing",
        "recommend_comparison_analysis": recommend_comparison,
        "has_nan": bool(md100.get("has_nan")) or bool(md50.get("has_nan")) or bool(gate.get("has_nan")),
    }
    row["stage_reached"] = stage_reached(row)
    return row


def collect_stage_rows() -> list[dict[str, object]]:
    return [row_for_system(system) for system in SYSTEMS]


def write_csv(rows: list[dict[str, object]]) -> None:
    with SUMMARY_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_report(rows: list[dict[str, object]]) -> None:
    lines = [
        f"# {LIGAND_ID} / {RESNAME} Gate Summary",
        "",
        f"This file summarizes the gate-based non-production restrained MD pilot for {LIGAND_ID} / {RESNAME} in 20chol and 40chol.",
        "",
        "These results are from non-production R2 PBC-aware restrained MD pilot runs only. They are not production MD conclusions.",
        "",
    ]
    for row in rows:
        lines.extend(
            [
                f"## {row['system']}",
                "",
                f"- stage reached: {row['stage_reached']}",
                f"- PSF/PDB generated: {row['PSF_PDB_generated']}",
                f"- read test succeeded: {row['read_test_succeeded']}",
                f"- initial energy: {row['initial_energy_kJ_mol']}",
                f"- severe/close before cleanup: {row['severe_before_cleanup']}/{row['close_before_cleanup']}",
                f"- cleanup executed: {row['cleanup_executed']}",
                f"- deleted molecules: {row['deleted_molecules']}",
                f"- minimization final energy: {row['minimization_final_energy_kJ_mol']}",
                f"- very-short probe passed: {row['very_short_probe_passed']}",
                f"- 10 ps completed: {row['completed_10ps']}",
                f"- 50 ps completed: {row['completed_50ps']}",
                f"- 100 ps completed: {row['completed_100ps']}",
                f"- final temperature: {row['final_temperature_K']}",
                f"- final potential energy: {row['final_potential_energy_kJ_mol']}",
                f"- final TYR13 / VAL17 / SER20 distances: {row['final_TYR13_distance_A']} / {row['final_VAL17_distance_A']} / {row['final_SER20_distance_A']}",
                f"- near chains count: {row['near_chains_count']}",
                f"- ligand RMSD final/max: {row['ligand_rmsd_final_A']} / {row['ligand_rmsd_max_A']}",
                f"- final severe/close contacts: {row['final_severe_contacts']}/{row['final_close_contacts']}",
                f"- NaN: {row['has_nan']}",
                f"- stop reason: {row['stop_reason']}",
                f"- recommend comparison analysis: {row['recommend_comparison_analysis']}",
                "",
            ]
        )
    SUMMARY_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows = collect_stage_rows()
    write_csv(rows)
    write_report(rows)
    for row in rows:
        print(f"{row['system']}: stage={row['stage_reached']} completed_100ps={row['completed_100ps']} stop_reason={row['stop_reason']}")


if __name__ == "__main__":
    main()
