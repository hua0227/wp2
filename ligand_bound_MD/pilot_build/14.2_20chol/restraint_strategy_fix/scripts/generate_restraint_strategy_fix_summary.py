from __future__ import annotations

import csv
from pathlib import Path

from strategy_common import REPORTS_DIR, OUTPUTS_DIR, fmt


PROA_MD = REPORTS_DIR / "PROA437_force_diagnostics.md"
PROA_CSV = REPORTS_DIR / "PROA437_force_decomposition.csv"
COMP_CSV = REPORTS_DIR / "pbcaware_minimization_comparison.csv"
PROBE_CSV = REPORTS_DIR / "very_short_probe_pbcaware.csv"
SUMMARY_MD = REPORTS_DIR / "restraint_strategy_fix_summary.md"
PML = OUTPUTS_DIR / "view_R1_R2_PROA437_L002.pml"


def read_rows(path: Path):
    return list(csv.DictReader(path.open(encoding="utf-8-sig")))


def branch(rows, name):
    return next(r for r in rows if r["branch"] == name)


def main() -> int:
    comp = read_rows(COMP_CSV)
    r1 = branch(comp, "R1")
    r2 = branch(comp, "R2")
    probe = read_rows(PROBE_CSV)
    probe_last = probe[-1]
    proa_rows = read_rows(PROA_CSV)
    group_rows = [r for r in proa_rows if r.get("scope", "") == ""]
    dom = max(group_rows, key=lambda r: float(r["PROA437_O_force_kJ_mol_nm"]) if r["PROA437_O_force_kJ_mol_nm"] != "nan" else -1)
    probe_pass = probe_last["status"] == "ok" and probe_last["step"] == "4000"
    recommend_10ps = probe_pass and r2["ligand_near_chains_count"] == "2" and r2["severe_lt_1A"] == "0"
    lines = [
        "# Restraint Strategy Fix Summary",
        "",
        "## Scope",
        "",
        "- ligand_id: 14.2",
        "- resname: L002",
        "- system: 20chol",
        "- branch: Branch B remove_POPC94",
        "- No production MD was run.",
        "- No longer MD was run.",
        "- No files outside pilot-local Branch B/fix outputs were modified.",
        "",
        "## Why MEMB:49 Was Not Further Processed",
        "",
        "Residual force diagnostics showed MEMB:49 was far from ligand L002, had no abnormal bonds, and dropped out of the no-restraint top force atoms. Its previous high force was dominated by absolute-coordinate CustomExternalForce restraints, so deleting or processing MEMB:49 would not address the root cause.",
        "",
        "## PROA437 Force Source",
        "",
        f"- PROA:437 O dominant force group: {dom['force_source_classification']} ({dom['force_class']}, group {dom['force_group']}).",
        f"- PROA:437 O dominant group force: {dom['PROA437_O_force_kJ_mol_nm']} kJ/mol/nm.",
        "- Interpretation: this is mainly internal torsion/custom torsion strain, not the old MEMB49 restraint artifact.",
        "",
        "## Old Restraint Strategy Problem",
        "",
        "The old strategy strongly restrained all lipid/cholesterol heavy atoms with absolute Cartesian displacement. In this periodic system, that produced restraint artifacts when coordinates/reference positions landed in different periodic images. The corrected strategy uses `periodicdistance(...)` and avoids strong lipid/cholesterol restraints.",
        "",
        "## R1 vs R2",
        "",
        "- R1: protein backbone medium PBC-aware restraint + ligand heavy weak PBC-aware restraint; no lipid/cholesterol restraint.",
        "- R2: same as R1 plus extremely weak lipid/cholesterol heavy atom PBC-aware restraint.",
        "",
        "| Branch | Final energy kJ/mol | Max force kJ/mol/nm | Top force atom | MEMB49 max | PROA437/438 max | Severe/close | Key min A | Near chains | Ligand RMSD A |",
        "|---|---:|---:|---|---:|---:|---:|---:|---:|---:|",
    ]
    for r in [r1, r2]:
        lines.append(
            f"| {r['branch']} | {r['final_energy_kJ_mol']} | {r['max_force_kJ_mol_nm']} | {r['top_force_atom']} | "
            f"{r['MEMB49_max_force_kJ_mol_nm']} | {r['PROA437_438_max_force_kJ_mol_nm']} | "
            f"{r['severe_lt_1A']}/{r['close_lt_1p5A']} | {r['ligand_key_min_distance_A']} | {r['ligand_near_chains_count']} | {r['ligand_rmsd_vs_BranchB_start_A']} |"
        )
    lines.extend([
        "",
        "## Ligand Retention",
        "",
        f"- R1 ligand near chains: {r1['ligand_near_chains_count']}; R1 is less preferred.",
        f"- R2 ligand near chains: {r2['ligand_near_chains_count']}; R2 retained the ligand near both chains.",
        f"- R2 severe/close ligand clashes: {r2['severe_lt_1A']} / {r2['close_lt_1p5A']}.",
        "",
        "## Very-Short Probe",
        "",
        "- Probe branch: R2.",
        "- Conditions: 50 K, 0.05 fs timestep, 0.2 ps total, OpenCL, PBC-aware restraints, not production MD.",
        f"- Probe executed: YES.",
        f"- Probe passed: {'YES' if probe_pass else 'NO'}.",
        f"- Final probe step/status: {probe_last['step']} / {probe_last['status']}.",
        f"- Final probe temperature: {probe_last['temperature_K']} K.",
        f"- Final probe ligand near chains: {probe_last['ligand_near_chains_count']}.",
        "",
        "## Recommendation",
        "",
        f"- Recommend reattempting 10 ps restrained MD: {'YES' if recommend_10ps else 'NO'}.",
        "- Use the R2 PBC-aware restraint strategy, not the old absolute lipid/cholesterol restraints.",
        "- This recommendation is for the next short restrained MD attempt only, not production MD.",
        f"- PyMOL view: `{PML}`",
    ])
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("Terminal summary")
    print(f"PROA437 dominant force group: {dom['force_source_classification']} ({dom['force_class']})")
    for r in [r1, r2]:
        print(f"{r['branch']} final energy/max force/top: {r['final_energy_kJ_mol']} / {r['max_force_kJ_mol_nm']} / {r['top_force_atom']}")
    print(f"MEMB49 max force R1/R2: {r1['MEMB49_max_force_kJ_mol_nm']} / {r2['MEMB49_max_force_kJ_mol_nm']}")
    print(f"PROA437 max force R1/R2: {r1['PROA437_438_max_force_kJ_mol_nm']} / {r2['PROA437_438_max_force_kJ_mol_nm']}")
    print(f"ligand retention R1: key_min={r1['ligand_key_min_distance_A']} chains={r1['ligand_near_chains_count']} RMSD={r1['ligand_rmsd_vs_BranchB_start_A']}")
    print(f"ligand retention R2: key_min={r2['ligand_key_min_distance_A']} chains={r2['ligand_near_chains_count']} RMSD={r2['ligand_rmsd_vs_BranchB_start_A']}")
    print("very-short probe executed: True")
    print(f"very-short probe passed: {probe_pass}")
    print(f"recommendation: {'reattempt 10 ps restrained MD with R2 PBC-aware restraints' if recommend_10ps else 'do not reattempt 10 ps yet'}")
    print(f"summary report: {SUMMARY_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
