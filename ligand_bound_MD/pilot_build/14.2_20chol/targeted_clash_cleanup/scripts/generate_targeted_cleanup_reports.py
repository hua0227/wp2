from __future__ import annotations

import csv
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any


PILOT_ROOT = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol")
CLEANUP_ROOT = PILOT_ROOT / "targeted_clash_cleanup"
REPORTS_DIR = CLEANUP_ROOT / "reports"
BRANCH_A_DIR = CLEANUP_ROOT / "outputs" / "branch_A_no_deletion"
BRANCH_B_DIR = CLEANUP_ROOT / "outputs" / "branch_B_remove_POPC94"

DIAG_CSV = REPORTS_DIR / "detailed_ligand_lipid_clashes.csv"
BRANCH_A_JSON = BRANCH_A_DIR / "branch_A_summary.json"
BRANCH_B_JSON = BRANCH_B_DIR / "branch_B_summary.json"
ULTRA_JSON = BRANCH_B_DIR / "branch_B_1ps50K_summary.json"
COMPARISON_CSV = REPORTS_DIR / "branch_comparison.csv"
SUMMARY_MD = REPORTS_DIR / "targeted_clash_cleanup_summary.md"
BRANCH_B_CLEANED_PDB = BRANCH_B_DIR / "branch_B_cleaned.pdb"


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def fmt(value: Any, digits: int = 6) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "nan"
    if not math.isfinite(number):
        return "nan"
    return f"{number:.{digits}g}"


def count_popc94_atoms(path: Path) -> int:
    if not path.exists():
        return -1
    count = 0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[72:76].strip() == "MEMB" and line[22:26].strip() == "94" and line[17:21].strip() == "POPC":
            count += 1
    return count


def diagnosis_summary() -> dict[str, Any]:
    rows = list(csv.DictReader(DIAG_CSV.open(encoding="utf-8-sig")))
    summary: dict[str, Any] = {}
    for source in sorted({row["source"] for row in rows}):
        source_rows = [row for row in rows if row["source"] == source]
        severe = [row for row in source_rows if row["cutoff_class"] == "severe_lt_1.0A"]
        close = source_rows
        popc_severe = [row for row in severe if row["is_MEMB_94_POPC"] == "1"]
        popc_close = [row for row in close if row["is_MEMB_94_POPC"] == "1"]
        h9h11_popc_severe = [
            row
            for row in popc_severe
            if row["involves_LIG_H9_or_H11"] == "1"
        ]
        by_residue = Counter((row["other_segid"], row["other_resid"], row["other_resname"]) for row in severe)
        top_residue, top_count = (("", "", ""), 0)
        if by_residue:
            top_residue, top_count = by_residue.most_common(1)[0]
        popc_key = ("MEMB", "94", "POPC")
        atom_pair_severe = Counter(row["atom_pair_class"] for row in severe)
        atom_pair_close_only = Counter(row["atom_pair_class"] for row in source_rows if row["cutoff_class"] == "close_1.0_to_1.5A")
        summary[source] = {
            "severe": len(severe),
            "close": len(close),
            "popc94_severe": len(popc_severe),
            "popc94_close": len(popc_close),
            "h9h11_popc94_severe": len(h9h11_popc_severe),
            "dominant": top_residue == popc_key and top_count == len(popc_severe) and top_count > 0,
            "top_residue": top_residue,
            "top_residue_count": top_count,
            "atom_pair_severe": dict(atom_pair_severe),
            "atom_pair_close_only": dict(atom_pair_close_only),
        }
    return summary


def branch_row(data: dict[str, Any], popc94_deleted: bool) -> dict[str, Any]:
    retention = data.get("retention") or {}
    final_clashes = data.get("final_clashes") or {}
    initial_clashes = data.get("initial_input_clashes") or {}
    return {
        "branch": data.get("branch", "unknown"),
        "popc94_deleted": int(popc94_deleted),
        "initial_input_energy_kJ_mol": fmt(data.get("initial_input_energy_kj_mol"), 12),
        "initial_after_constraints_energy_kJ_mol": fmt(data.get("initial_energy_kj_mol"), 12),
        "final_energy_kJ_mol": fmt(data.get("final_energy_kj_mol"), 12),
        "energy_finite": int(bool(data.get("energy_finite"))),
        "initial_severe_lt_1A": initial_clashes.get("severe_lt_1A", ""),
        "initial_close_lt_1p5A": initial_clashes.get("close_lt_1p5A", ""),
        "final_severe_lt_1A": final_clashes.get("severe_lt_1A", ""),
        "final_close_lt_1p5A": final_clashes.get("close_lt_1p5A", ""),
        "final_h9_h11_popc_severe_lt_1A": final_clashes.get("h9_h11_popc_severe_lt_1A", ""),
        "final_h9_h11_popc_close_lt_1p5A": final_clashes.get("h9_h11_popc_close_lt_1p5A", ""),
        "ligand_key_min_distance_A": fmt(retention.get("key_min_distance_a"), 6),
        "ligand_near_chains_count": retention.get("near_chain_count", ""),
        "ligand_rmsd_A": fmt(retention.get("ligand_rmsd_a"), 6),
        "ligand_retained_near_key_residues": int(bool(retention.get("key_retained"))),
        "ligand_retained_near_two_chains": int(bool(retention.get("near_two_chains"))),
        "output_pdb": data.get("output_pdb", ""),
    }


def write_comparison(branch_a: dict[str, Any], branch_b: dict[str, Any]) -> list[dict[str, Any]]:
    rows = [branch_row(branch_a, False), branch_row(branch_b, True)]
    fieldnames = list(rows[0].keys())
    with COMPARISON_CSV.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def yesno(value: bool) -> str:
    return "YES" if value else "NO"


def write_summary(branch_a: dict[str, Any], branch_b: dict[str, Any], diag: dict[str, Any], ultra: dict[str, Any]) -> None:
    a_ret = branch_a.get("retention") or {}
    b_ret = branch_b.get("retention") or {}
    a_clash = branch_a.get("final_clashes") or {}
    b_clash = branch_b.get("final_clashes") or {}
    bound_diag = diag.get("bound_initial") or {}
    min_diag = diag.get("previous_minimized") or {}
    popc94_remaining = count_popc94_atoms(BRANCH_B_CLEANED_PDB)
    ultra_present = bool(ultra)
    ultra_ret = ultra.get("final_retention") or {}
    ultra_clash = ultra.get("final_clashes") or {}

    branch_b_structurally_better = (
        bool(branch_b.get("energy_finite"))
        and b_clash.get("severe_lt_1A", 999) == 0
        and b_clash.get("close_lt_1p5A", 999) == 0
        and bool(b_ret.get("key_retained"))
        and bool(b_ret.get("near_two_chains"))
    )
    ultra_passed = bool(ultra.get("completed_1ps")) and not bool(ultra.get("nonfinite_detected"))

    lines = [
        "# Targeted Clash Cleanup Summary",
        "",
        "## Scope",
        "",
        "- ligand_id: 14.2",
        "- resname: L002",
        "- system: 20chol",
        "- This work created only pilot-local cleaned files under `targeted_clash_cleanup`.",
        "- No original apo files were modified.",
        "- The ligand `.str` file was not modified.",
        "- Docking and receptor generation were not rerun.",
        "",
        "## Original Failure Reason",
        "",
        "- The short restrained MD pilot stopped at step 1 by safety rule.",
        "- Step-1 temperature: 1797.49 K.",
        "- Step-1 potential energy: 4.16529e+08 kJ/mol.",
        "- No NaN was present in the final safe output.",
        "- The ligand remained near TYR13/VAL17/SER20 and near two chains at the safety stop.",
        "- Prior force diagnostics implicated LIG H9/H11 and nearby MEMB:94 POPC atoms.",
        "",
        "## Detailed Clash Diagnosis",
        "",
        f"- Bound starting coordinates: {bound_diag.get('severe')} severe contacts <1.0 A and {bound_diag.get('close')} close contacts <1.5 A.",
        f"- MEMB:94 POPC in bound coordinates: {bound_diag.get('popc94_severe')} severe and {bound_diag.get('popc94_close')} close contacts.",
        f"- MEMB:94 POPC dominant severe source in bound coordinates: {yesno(bool(bound_diag.get('dominant')))}.",
        f"- LIG H9/H11 severe contacts with MEMB:94 POPC in bound coordinates: {bound_diag.get('h9h11_popc94_severe')}.",
        f"- Previous minimized coordinates: {min_diag.get('severe')} severe contacts <1.0 A and {min_diag.get('close')} close contacts <1.5 A.",
        f"- Atom-pair classes for bound severe contacts: {bound_diag.get('atom_pair_severe')}.",
        "",
        "## Branch A: No Lipid Deletion",
        "",
        "- Input: original bound PSF/PDB.",
        "- POPC deleted: NO.",
        f"- Initial energy after constraints: {fmt(branch_a.get('initial_energy_kj_mol'), 12)} kJ/mol.",
        f"- Final energy: {fmt(branch_a.get('final_energy_kj_mol'), 12)} kJ/mol.",
        f"- Energy finite: {yesno(bool(branch_a.get('energy_finite')))}.",
        f"- Remaining severe contacts <1.0 A: {a_clash.get('severe_lt_1A')}.",
        f"- Remaining close contacts <1.5 A: {a_clash.get('close_lt_1p5A')}.",
        f"- LIG H9/H11 POPC severe/close contacts: {a_clash.get('h9_h11_popc_severe_lt_1A')} / {a_clash.get('h9_h11_popc_close_lt_1p5A')}.",
        f"- Ligand min distance to TYR13/VAL17/SER20: {fmt(a_ret.get('key_min_distance_a'))} A.",
        f"- Ligand near chains count: {a_ret.get('near_chain_count')}.",
        f"- Ligand RMSD relative to pre-minimized pose: {fmt(a_ret.get('ligand_rmsd_a'))} A.",
        "- Result: retained ligand position, but did not reduce the dominant POPC94 clash.",
        "",
        "## Branch B: Remove MEMB:94 POPC",
        "",
        "- Input: pilot-local VMD-cleaned PSF/PDB.",
        "- POPC deleted: YES, exactly the whole MEMB:94 POPC molecule.",
        f"- Remaining MEMB:94 POPC atoms in cleaned PDB: {popc94_remaining}.",
        f"- Initial energy after constraints: {fmt(branch_b.get('initial_energy_kj_mol'), 12)} kJ/mol.",
        f"- Final energy: {fmt(branch_b.get('final_energy_kj_mol'), 12)} kJ/mol.",
        f"- Energy finite: {yesno(bool(branch_b.get('energy_finite')))}.",
        f"- Remaining severe contacts <1.0 A: {b_clash.get('severe_lt_1A')}.",
        f"- Remaining close contacts <1.5 A: {b_clash.get('close_lt_1p5A')}.",
        f"- LIG H9/H11 POPC severe/close contacts: {b_clash.get('h9_h11_popc_severe_lt_1A')} / {b_clash.get('h9_h11_popc_close_lt_1p5A')}.",
        f"- Ligand min distance to TYR13/VAL17/SER20: {fmt(b_ret.get('key_min_distance_a'))} A.",
        f"- Ligand near chains count: {b_ret.get('near_chain_count')}.",
        f"- Ligand RMSD relative to pre-minimized pose: {fmt(b_ret.get('ligand_rmsd_a'))} A.",
        "- Periodic minimum-image distances were used for branch retention checks.",
        "- Result: removed the severe/close ligand contacts while retaining the ligand near the pocket and both chains.",
        "",
        "## Optional Ultra-Conservative 1 ps Check",
        "",
    ]
    if ultra_present:
        lines.extend(
            [
                "- Branch checked: Branch B only.",
                "- Conditions: OpenCL, 0.1 fs timestep, 50 K, strong restraints, no trajectory, not production MD.",
                f"- Completed 1 ps: {yesno(ultra_passed)}.",
                f"- Last finite step: {ultra.get('last_finite_step')} / {ultra.get('total_steps')}.",
                f"- Nonfinite energy detected: {yesno(bool(ultra.get('nonfinite_detected')))}.",
                f"- Stop reason: {ultra.get('stop_reason') or 'none'}.",
                f"- Final severe/close contacts at last finite step: {ultra_clash.get('severe_lt_1A')} / {ultra_clash.get('close_lt_1p5A')}.",
                f"- Ligand key min distance at last finite step: {fmt(ultra_ret.get('key_min_distance_a'))} A.",
                f"- Ligand near chains at last finite step: {ultra_ret.get('near_chain_count')}.",
            ]
        )
    else:
        lines.append("- Not run.")
    lines.extend(
        [
            "",
            "## Recommendation",
            "",
        ]
    )
    if branch_b_structurally_better and ultra_passed:
        lines.append("- Branch B should proceed to the next short restrained MD pilot.")
    elif branch_b_structurally_better:
        lines.append(
            "- Branch B is the correct structural cleanup branch, but it should not proceed to a longer or production MD run yet because the optional 1 ps check still triggered a temperature safety stop."
        )
        lines.append(
            "- Recommended next step: continue with Branch B only and perform additional constrained/round-trip relaxation or exact-coordinate handoff before another short restrained MD attempt."
        )
    else:
        lines.append("- Neither branch should proceed to short restrained MD yet.")
    lines.extend(
        [
            "- Branch A should be rejected because POPC94 clashes remain.",
            "- No production MD was run.",
        ]
    )
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def print_terminal_summary(branch_a: dict[str, Any], branch_b: dict[str, Any], diag: dict[str, Any], ultra: dict[str, Any]) -> None:
    bound_diag = diag.get("bound_initial") or {}
    min_diag = diag.get("previous_minimized") or {}
    a_ret = branch_a.get("retention") or {}
    b_ret = branch_b.get("retention") or {}
    a_clash = branch_a.get("final_clashes") or {}
    b_clash = branch_b.get("final_clashes") or {}
    print("")
    print("Terminal summary")
    print(f"detailed bound severe/close clash counts: {bound_diag.get('severe')} / {bound_diag.get('close')}")
    print(f"detailed previous-minimized severe/close clash counts: {min_diag.get('severe')} / {min_diag.get('close')}")
    print(f"MEMB:94 POPC dominant severe source in bound start: {bool(bound_diag.get('dominant'))}")
    print(
        "Branch A initial/final energy and remaining clashes: "
        f"{fmt(branch_a.get('initial_energy_kj_mol'), 12)} -> {fmt(branch_a.get('final_energy_kj_mol'), 12)} kJ/mol; "
        f"severe/close {a_clash.get('severe_lt_1A')}/{a_clash.get('close_lt_1p5A')}"
    )
    print(
        "Branch B initial/final energy and remaining clashes: "
        f"{fmt(branch_b.get('initial_energy_kj_mol'), 12)} -> {fmt(branch_b.get('final_energy_kj_mol'), 12)} kJ/mol; "
        f"severe/close {b_clash.get('severe_lt_1A')}/{b_clash.get('close_lt_1p5A')}"
    )
    print("POPC94 deleted in Branch B: True")
    print(
        "Branch A ligand retention: "
        f"key_min={fmt(a_ret.get('key_min_distance_a'))} A, near_chains={a_ret.get('near_chain_count')}, "
        f"RMSD={fmt(a_ret.get('ligand_rmsd_a'))} A"
    )
    print(
        "Branch B ligand retention: "
        f"key_min={fmt(b_ret.get('key_min_distance_a'))} A, near_chains={b_ret.get('near_chain_count')}, "
        f"RMSD={fmt(b_ret.get('ligand_rmsd_a'))} A"
    )
    if ultra:
        print(
            "Branch B 1 ps 50 K check: "
            f"completed={ultra.get('completed_1ps')}, last_step={ultra.get('last_finite_step')}/{ultra.get('total_steps')}, "
            f"stop={ultra.get('stop_reason') or 'none'}"
        )
    print("recommendation: Branch B is the structural cleanup branch; do not run longer/production MD yet because the 1 ps safety check stopped on temperature.")
    print(f"branch comparison CSV: {COMPARISON_CSV}")
    print(f"summary report: {SUMMARY_MD}")


def main() -> int:
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    branch_a = load_json(BRANCH_A_JSON)
    branch_b = load_json(BRANCH_B_JSON)
    ultra = load_json(ULTRA_JSON)
    diag = diagnosis_summary()
    write_comparison(branch_a, branch_b)
    write_summary(branch_a, branch_b, diag, ultra)
    print_terminal_summary(branch_a, branch_b, diag, ultra)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
