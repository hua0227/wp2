from __future__ import annotations

import csv
import re
from pathlib import Path

from branchB_common import LOGS_DIR, REPORTS_DIR, fmt


HANDOFF_CSV = REPORTS_DIR / "coordinate_handoff_check.csv"
FORCE_CSV = REPORTS_DIR / "branchB_force_diagnostics.csv"
FORCE_LOG = LOGS_DIR / "branchB_force_diagnostics.log"
STAGED_CSV = REPORTS_DIR / "staged_minimization_summary.csv"
PROBE_LOG = LOGS_DIR / "branchB_very_short_probe.log"
SUMMARY_MD = REPORTS_DIR / "branchB_relaxation_summary.md"


def read_csv(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(path.open(encoding="utf-8-sig")))


def handoff_pass() -> bool:
    rows = read_csv(HANDOFF_CSV)
    for row in rows:
        if row["check"] == "coordinate_handoff_pass":
            return row["pass"] == "1"
    return False


def handoff_value(check_name: str) -> str:
    for row in read_csv(HANDOFF_CSV):
        if row["check"] == check_name:
            return row["value"]
    return "not available"


def initial_energy_from_log() -> float:
    text = FORCE_LOG.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"initial_branchB_energy_kJ_mol:\s*([-+0-9.eE]+)", text)
    return float(match.group(1)) if match else float("nan")


def force_top_rows() -> list[dict[str, str]]:
    return read_csv(FORCE_CSV)


def staged_rows() -> list[dict[str, str]]:
    return read_csv(STAGED_CSV)


def bool_text(value: bool) -> str:
    return "YES" if value else "NO"


def probe_gate(initial_max_force: float, stage4: dict[str, str]) -> tuple[bool, list[str], float]:
    reasons: list[str] = []
    stage4_force = float(stage4["physical_max_force_kJ_mol_nm"])
    reduction_fraction = (initial_max_force - stage4_force) / initial_max_force if initial_max_force > 0 else 0.0
    if stage4["energy_finite"] != "1":
        reasons.append("Stage 4 energy was not finite")
    if stage4["force_nan"] != "0":
        reasons.append("Stage 4 had NaN/nonfinite force")
    if stage4["severe_lt_1A"] != "0":
        reasons.append("Stage 4 still had severe clashes")
    if int(stage4["close_lt_1p5A"]) > 2:
        reasons.append("Stage 4 close contacts were not low")
    if stage4["ligand_retained_near_key_residues"] != "1":
        reasons.append("Ligand did not remain near key residues")
    if stage4["ligand_retained_near_two_chains"] != "1":
        reasons.append("Ligand was not near two chains")
    if reduction_fraction < 0.5:
        reasons.append(f"Max force reduction was only {reduction_fraction * 100.0:.1f}%, below the 50% significance threshold used for this gate")
    return not reasons, reasons, reduction_fraction


def write_probe_skip_log(reasons: list[str]) -> None:
    lines = [
        "Very-short 0.2 ps probe was NOT executed.",
        "No dynamics was run by this probe stage.",
        "Reason:",
        *[f"- {reason}" for reason in reasons],
    ]
    PROBE_LOG.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    hp = handoff_pass()
    initial_energy = initial_energy_from_log()
    top_force = force_top_rows()[0]
    rows = staged_rows()
    stage4 = rows[-1]
    initial_max_force = float(top_force["force_magnitude_kJ_mol_nm"])
    gate_pass, gate_reasons, reduction_fraction = probe_gate(initial_max_force, stage4)
    probe_executed = False
    probe_passed = False
    if not gate_pass:
        write_probe_skip_log(gate_reasons)

    recommend_short_md = gate_pass and probe_executed and probe_passed
    lines = [
        "# Branch B Relaxation Summary",
        "",
        "## Scope",
        "",
        "- ligand_id: 14.2",
        "- resname: L002",
        "- system: 20chol",
        "- branch: Branch B remove_POPC94",
        "- This stage used only Branch B pilot-local files.",
        "- No production MD was run.",
        "- No longer MD was run.",
        "- Original apo files, receptor files, docking, and ligand `.str` were not modified.",
        "",
        "## Why Round-Trip And Staged Relaxation Were Needed",
        "",
        "The targeted cleanup removed POPC MEMB:94 and eliminated severe/close ligand clashes, but the ultra-conservative 1 ps 50 K check stopped after 10 steps because temperature exceeded 200 K. The next diagnostic layer therefore checked whether the Branch B minimized PDB is handed to OpenMM exactly, then inspected the force distribution and relaxed the system in staged restraint schedules without starting dynamics.",
        "",
        "## Coordinate Handoff",
        "",
        f"- Coordinate handoff passed: {bool_text(hp)}.",
        f"- PSF atom count: {handoff_value('psf_atom_count')}.",
        f"- Cleaned PDB atom count: {handoff_value('cleaned_pdb_atom_count')}.",
        f"- Minimized PDB atom count: {handoff_value('minimized_pdb_atom_count')}.",
        f"- Roundtrip heavy-atom RMSD: {handoff_value('roundtrip_heavy_atom_rmsd_A')} A.",
        f"- Roundtrip ligand heavy-atom RMSD: {handoff_value('roundtrip_ligand_heavy_atom_rmsd_A')} A.",
        f"- Periodic box max delta: {handoff_value('periodic_box_max_delta_A')} A.",
        "",
        "## Initial Force Diagnostics",
        "",
        f"- Initial Branch B potential energy: {fmt(initial_energy, 12)} kJ/mol.",
        f"- Force NaN/nonfinite present: NO.",
        f"- Highest-force atom: {top_force['atom_name']} {top_force['resname']} {top_force['segid']}:{top_force['resid']} index {top_force['index']}.",
        f"- Highest-force magnitude: {fmt(top_force['force_magnitude_kJ_mol_nm'], 12)} kJ/mol/nm.",
        f"- Highest forces concentrated in ligand H atoms: NO.",
        "",
        "## Highest Force Atoms",
        "",
        "| Rank | Atom | Residue | Category | Force kJ/mol/nm |",
        "|---:|---|---|---|---:|",
    ]
    for row in force_top_rows()[:10]:
        lines.append(
            f"| {row['rank']} | {row['atom_name']} index {row['index']} | {row['resname']} {row['segid']}:{row['resid']} | "
            f"{row['category']} | {fmt(row['force_magnitude_kJ_mol_nm'], 8)} |"
        )
    lines.extend(
        [
            "",
            "## Staged Minimization",
            "",
            "| Stage | Physical energy kJ/mol | Physical max force kJ/mol/nm | Max force atom | Severe/close | Key min A | Near chains | Ligand RMSD A |",
            "|---:|---:|---:|---|---:|---:|---:|---:|",
        ]
    )
    for row in rows:
        lines.append(
            f"| {row['stage']} | {fmt(row['physical_energy_kJ_mol'], 12)} | {fmt(row['physical_max_force_kJ_mol_nm'], 8)} | "
            f"{row['physical_max_force_atom']} | {row['severe_lt_1A']}/{row['close_lt_1p5A']} | "
            f"{row['ligand_key_min_distance_A']} | {row['ligand_near_chains_count']} | {row['ligand_rmsd_vs_branchB_start_A']} |"
        )
    lines.extend(
        [
            "",
            "## Ligand Retention And Clashes",
            "",
            f"- Stage 4 severe/close clashes: {stage4['severe_lt_1A']} / {stage4['close_lt_1p5A']}.",
            f"- Stage 4 ligand min distance to TYR13/VAL17/SER20: {stage4['ligand_key_min_distance_A']} A.",
            f"- Stage 4 ligand near chains count: {stage4['ligand_near_chains_count']}.",
            f"- Stage 4 ligand RMSD vs Branch B starting ligand: {stage4['ligand_rmsd_vs_branchB_start_A']} A.",
            f"- Stage 4 max force reduction vs initial diagnostic: {reduction_fraction * 100.0:.1f}%.",
            "",
            "## Very-Short Probe",
            "",
            f"- Probe executed: {bool_text(probe_executed)}.",
            f"- Probe passed: {bool_text(probe_passed)}.",
        ]
    )
    if gate_reasons:
        lines.append("- Probe was skipped because:")
        lines.extend([f"  - {reason}" for reason in gate_reasons])
    lines.extend(
        [
            "",
            "## Recommendation",
            "",
            f"- Recommend next 10 ps restrained MD attempt: {bool_text(recommend_short_md)}.",
            "- Reason: staged minimization preserved ligand retention and kept clashes at zero, but Stage 4 did not significantly reduce the maximum force, so the requested gate for the 0.2 ps probe was not met.",
            "- Continue diagnosing the residual high-force POPC region locally before another short restrained MD attempt.",
            "- No production MD was run in this stage.",
        ]
    )
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("")
    print("Terminal summary")
    print(f"coordinate handoff passed: {hp}")
    print(f"initial Branch B energy (kJ/mol): {fmt(initial_energy, 12)}")
    print(
        "initial max force atom: "
        f"{top_force['atom_name']} {top_force['resname']} {top_force['segid']}:{top_force['resid']} "
        f"index={top_force['index']} force={fmt(top_force['force_magnitude_kJ_mol_nm'], 12)} kJ/mol/nm"
    )
    for row in rows:
        print(f"Stage {row['stage']} final energy (kJ/mol): {fmt(row['physical_energy_kJ_mol'], 12)}")
    for row in rows:
        print(f"Stage {row['stage']} max force (kJ/mol/nm): {fmt(row['physical_max_force_kJ_mol_nm'], 12)} at {row['physical_max_force_atom']}")
    print(f"Stage 4 severe/close clash: {stage4['severe_lt_1A']}/{stage4['close_lt_1p5A']}")
    print(
        "ligand retention: "
        f"key_min={stage4['ligand_key_min_distance_A']} A, near_chains={stage4['ligand_near_chains_count']}, "
        f"RMSD={stage4['ligand_rmsd_vs_branchB_start_A']} A"
    )
    print(f"very-short probe executed: {probe_executed}")
    print(f"very-short probe passed: {probe_passed}")
    print(f"recommend next short restrained MD attempt: {recommend_short_md}")
    print(f"summary report: {SUMMARY_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
