from __future__ import annotations

import csv
import json
from pathlib import Path

from residual_common import OUTPUTS_DIR, REPORTS_DIR, fmt


FORCE_DECOMP_SUMMARY = REPORTS_DIR / "force_decomposition_summary.csv"
FORCE_DECOMP_TOP50 = REPORTS_DIR / "force_decomposition_top50.csv"
GEOM_SUMMARY_JSON = REPORTS_DIR / "MEMB49_geometry_summary.json"
NO_RESTRAINT_TOP50 = REPORTS_DIR / "no_restraint_force_top50.csv"
NO_RESTRAINT_MD = REPORTS_DIR / "no_restraint_force_summary.md"
SUMMARY_MD = REPORTS_DIR / "residual_force_diagnostics_summary.md"
PML_OUT = OUTPUTS_DIR / "view_L002_MEMB49.pml"


def rows(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(path.open(encoding="utf-8-sig")))


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_no_restraint_summary_value(prefix: str) -> str:
    for line in NO_RESTRAINT_MD.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(prefix):
            return line[len(prefix) :].strip().lstrip(":").strip().strip(".")
    return "not available"


def main() -> int:
    decomp = rows(FORCE_DECOMP_SUMMARY)
    top50 = rows(FORCE_DECOMP_TOP50)
    geom = load_json(GEOM_SUMMARY_JSON)
    no_top50 = rows(NO_RESTRAINT_TOP50)
    dominant = max(decomp, key=lambda row: float(row["MEMB49_max_force_kJ_mol_nm"]))
    non_restraint_decomp = [row for row in decomp if row["force_label"] != "CustomExternalForce/restraints"]
    dominant_non_restraint = max(non_restraint_decomp, key=lambda row: float(row["MEMB49_max_force_kJ_mol_nm"]))
    no_restraint_top = no_top50[0]
    no_restraint_memb49 = [row for row in no_top50 if row["is_MEMB49"] == "1"]
    memb49_in_no_restraint_top50 = bool(no_restraint_memb49)
    no_restraint_memb49_max = read_no_restraint_summary_value("- MEMB:49 max force magnitude")
    no_restraint_memb49_atom = read_no_restraint_summary_value("- MEMB:49 max force atom")
    no_restraint_energy = read_no_restraint_summary_value("- Potential energy")
    recommendations = {
        "A_continue_minimization": "NO - MEMB:49 has no abnormal bonds and is not the no-restraint top force source; blind minimization is unlikely to target the root cause.",
        "B_adjust_restraints": "YES - the with-restraint MEMB:49 force is dominated by CustomExternalForce/restraints, consistent with absolute position restraints interacting with periodic/wrapped coordinates.",
        "C_delete_or_process_MEMB49": "NO - MEMB:49 is far from ligand, has no abnormal bonds, and drops out of the no-restraint top 50.",
        "D_abandon_14_2_pilot": "NO - ligand pose remains clean; diagnose restraint/reference imaging and the no-restraint top protein force before abandoning.",
    }
    lines = [
        "# Residual Force Diagnostics Summary",
        "",
        "## Scope",
        "",
        "- ligand_id: 14.2",
        "- resname: L002",
        "- system: 20chol",
        "- branch: Branch B remove_POPC94",
        "- This stage is diagnostic only; no structure cleanup or deletion was performed.",
        "- No MD of any kind was run.",
        "- No original apo files, receptor files, docking outputs, or ligand `.str` files were modified.",
        "",
        "## MEMB:49 High Force Status",
        "",
        f"- MEMB:49 high force with Stage 4 restraints: YES, max {dominant['MEMB49_max_force_kJ_mol_nm']} kJ/mol/nm.",
        f"- Dominant with-restraint force group: {dominant['force_label']} ({dominant['force_class']}, group {dominant['force_group']}).",
        f"- Dominant with-restraint MEMB:49 atom: {dominant['MEMB49_max_force_atom']}.",
        f"- Dominant non-restraint MEMB:49 force group: {dominant_non_restraint['force_label']} ({dominant_non_restraint['force_class']}), max {dominant_non_restraint['MEMB49_max_force_kJ_mol_nm']} kJ/mol/nm.",
        f"- No-restraint MEMB:49 max force: {no_restraint_memb49_max}.",
        f"- No-restraint MEMB:49 in top 50: {'YES' if memb49_in_no_restraint_top50 else 'NO'}.",
        "",
        "## Force Source Interpretation",
        "",
        "- Restraints: YES, the MEMB:49 high force is dominated by the CustomExternalForce restraint group when Stage 4-style restraints are present.",
        f"- Bond/angle/internal geometry: no abnormal MEMB:49 bonds were found; the largest non-restraint MEMB:49 group is {dominant_non_restraint['force_label']}.",
        f"- Nonbonded contacts: nonbonded contribution exists but is not the dominant with-restraint source; no-restraint top force atom is {no_restraint_top['atom_name']} {no_restraint_top['resname']} {no_restraint_top['segid']}:{no_restraint_top['resid']} index {no_restraint_top['index']}.",
        "- Ligand-lipid interaction: not supported by the distance/contact diagnostics.",
        "",
        "## MEMB:49 Geometry And Contacts",
        "",
        f"- Abnormal MEMB:49 bond count: {geom['abnormal_bond_count']}.",
        f"- MEMB:49 to ligand L002 minimum distance: {fmt(geom['ligand_min_distance_A'], 6)} A.",
        f"- MEMB:49 still clashes with ligand: {'YES' if float(geom['ligand_min_distance_A']) < 1.5 else 'NO'}.",
        f"- MEMB:49 has NaN coordinates: {'YES' if geom['has_nan_coords'] else 'NO'}.",
        f"- MEMB:49 near box boundary: {'YES' if geom['near_boundary'] else 'NO'}; min boundary margin {fmt(geom['box_margin_A'], 6)} A.",
        f"- Original apo MEMB:49 heavy-atom RMSD: {fmt(geom['apo_heavy_rmsd_A'], 6)} A; mapped atoms {geom['apo_mapped_heavy_atom_count']}.",
        "",
        "## Top With-Restraint Force Atoms",
        "",
        "| Rank | Atom | Residue | Category | Force kJ/mol/nm | MEMB49 |",
        "|---:|---|---|---|---:|---:|",
    ]
    for row in top50[:10]:
        lines.append(
            f"| {row['rank']} | {row['atom_name']} index {row['index']} | {row['resname']} {row['segid']}:{row['resid']} | "
            f"{row['category']} | {fmt(row['force_magnitude_kJ_mol_nm'], 8)} | {row['is_MEMB49']} |"
        )
    lines.extend(
        [
            "",
            "## Top No-Restraint Force Atoms",
            "",
            "| Rank | Atom | Residue | Category | Force kJ/mol/nm | MEMB49 |",
            "|---:|---|---|---|---:|---:|",
        ]
    )
    for row in no_top50[:10]:
        lines.append(
            f"| {row['rank']} | {row['atom_name']} index {row['index']} | {row['resname']} {row['segid']}:{row['resid']} | "
            f"{row['category']} | {fmt(row['force_magnitude_kJ_mol_nm'], 8)} | {row['is_MEMB49']} |"
        )
    lines.extend(
        [
            "",
            "## Recommendations",
            "",
            f"- A. Continue minimization: {recommendations['A_continue_minimization']}",
            f"- B. Adjust restraints: {recommendations['B_adjust_restraints']}",
            f"- C. Locally delete/process MEMB:49: {recommendations['C_delete_or_process_MEMB49']}",
            f"- D. Abandon current 14.2 pilot: {recommendations['D_abandon_14_2_pilot']}",
            "",
            "Recommended next step: fix the restraint/reference-coordinate handling before any new short restrained MD attempt. In particular, avoid absolute lipid position restraints on coordinates that may be imaged across the periodic box, or regenerate restraint reference coordinates in the same periodic image as the OpenMM Context. Also inspect the no-restraint top protein force around PROA:437 before restarting a stability probe.",
            "",
            f"PyMOL view script: `{PML_OUT}`",
            "",
            "No MD was run in this diagnostic stage.",
        ]
    )
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("Terminal summary")
    print(f"MEMB49 max force with restraints: {dominant['MEMB49_max_force_kJ_mol_nm']} kJ/mol/nm")
    print(f"MEMB49 max force without restraints: {no_restraint_memb49_max}")
    print(f"dominant force group: {dominant['force_label']} ({dominant['force_class']})")
    print(f"MEMB49 abnormal bonds: {geom['abnormal_bond_count']}")
    print(f"MEMB49 ligand min distance: {fmt(geom['ligand_min_distance_A'], 6)} A")
    print(f"no-restraint top force atom: {no_restraint_top['atom_name']} {no_restraint_top['resname']} {no_restraint_top['segid']}:{no_restraint_top['resid']} index={no_restraint_top['index']} force={fmt(no_restraint_top['force_magnitude_kJ_mol_nm'], 12)} kJ/mol/nm")
    print("recommendation: adjust restraint/reference-coordinate handling; do not delete MEMB49 and do not restart short MD yet.")
    print(f"summary report: {SUMMARY_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
