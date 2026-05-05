from __future__ import annotations

import importlib.util
from pathlib import Path


BASE_SCRIPT = Path(r"C:\TRKB_WP2\ligand_bound_MD\WP2_Gate1_summary\scripts\generate_WP2_Gate1_summary.py")
SPEC = importlib.util.spec_from_file_location("generate_WP2_Gate1_summary_base", BASE_SCRIPT)
BASE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(BASE)

BASE.OUT_DIR = BASE.BASE_DIR / "WP2_Gate1_summary_with_9.2"
BASE.SCRIPTS_DIR = BASE.OUT_DIR / "scripts"
BASE.TABLES_DIR = BASE.OUT_DIR / "tables"
BASE.FIGURES_DIR = BASE.OUT_DIR / "figures"
BASE.REPORTS_DIR = BASE.OUT_DIR / "reports"
BASE.PYMOL_DIR = BASE.OUT_DIR / "pymol"
BASE.CANDIDATES = {
    "14.2": "L002",
    "19.1": "L003",
    "2.3": "L006",
    "6.2": "L008",
    "17.2": "L010",
    "9.2": "L004",
}
BASE.CANDIDATE_SUMMARY_CSV = BASE.TABLES_DIR / "WP2_Gate1_candidate_summary.csv"
BASE.LIGAND_LEVEL_CSV = BASE.TABLES_DIR / "WP2_Gate1_ligand_level_summary.csv"
BASE.REPORT_READY_CSV = BASE.TABLES_DIR / "WP2_Gate1_report_ready_summary.csv"
BASE.REPORT_PATH = BASE.REPORTS_DIR / "WP2_Gate1_summary_report.md"
BASE.PYMOL_PATH = BASE.PYMOL_DIR / "view_gate1_final_structures.pml"

ORIGINAL_WRITE_REPORT = BASE.write_report


OUT_DIR = BASE.OUT_DIR
CANDIDATES = BASE.CANDIDATES
classify_contact_mode = BASE.classify_contact_mode
classify_cholesterol_response = BASE.classify_cholesterol_response
classify_priority = BASE.classify_priority


def write_pymol_script() -> None:
    candidate_colors = {
        ("14.2", "20chol"): "marine",
        ("14.2", "40chol"): "bluewhite",
        ("19.1", "20chol"): "forest",
        ("19.1", "40chol"): "limegreen",
        ("2.3", "20chol"): "violet",
        ("2.3", "40chol"): "magenta",
        ("6.2", "20chol"): "orange",
        ("6.2", "40chol"): "tv_orange",
        ("17.2", "20chol"): "teal",
        ("17.2", "40chol"): "cyan",
        ("9.2", "20chol"): "salmon",
        ("9.2", "40chol"): "red",
    }
    lines = ["reinitialize"]
    for candidate_id, resname in BASE.CANDIDATES.items():
        for system in BASE.SYSTEMS:
            pdb = BASE.locate_final_pdb(candidate_id, resname, system)
            obj = f"{candidate_id.replace('.', '_')}_{system}_{resname}"
            lines.append(f"load {pdb.as_posix()}, {obj}")
            lines.append(f"color {candidate_colors[(candidate_id, system)]}, {obj}")
            lines.append(f"show sticks, {obj} and resn {resname}")
    lines.extend(
        [
            "hide everything",
            "show cartoon, polymer.protein",
            "show sticks, resn L002+L003+L004+L006+L008+L010",
            "show sticks, (resn TYR and resi 434) or (resn VAL and resi 438) or (resn SER and resi 441)",
            "color yellow, resn TYR+VAL+SER and resi 434+438+441",
            "set cartoon_transparency, 0.45",
            "set stick_radius, 0.18",
            "zoom (resn L002+L003+L004+L006+L008+L010) or (resi 434+438+441), 10",
            "orient (resn L002+L003+L004+L006+L008+L010) or (resi 434+438+441)",
            "# Visualization only; not used for simulation.",
            "",
        ]
    )
    BASE.PYMOL_PATH.write_text("\n".join(lines), encoding="utf-8")


def write_report(candidate_rows, ligand_rows) -> None:
    ORIGINAL_WRITE_REPORT(candidate_rows, ligand_rows)
    text = BASE.REPORT_PATH.read_text(encoding="utf-8")
    text = text.replace("总结五个 moderate-penalty 候选物", "总结六个 moderate-penalty 候选物")
    text = text.replace("- 五个候选物的 20chol / 40chol 体系都完成了 100 ps。", "- 六个候选物的 20chol / 40chol 体系都完成了 100 ps。")
    text = text.replace("- near chains 在 10 个 candidate-system 中都保持为 2。", "- near chains 在 12 个 candidate-system 中都保持为 2。")
    text = text.replace("- severe / close clash 在五个候选物的 comparison 表中都没有再次出现。", "- severe / close clash 在六个候选物的 comparison 表中都没有再次出现。")
    BASE.REPORT_PATH.write_text(text, encoding="utf-8")
    priority_a = [f"{row['candidate_id']} / {row['resname']}" for row in ligand_rows if row["WP2_Gate1_priority_class"] == "Priority A"]
    nine_two = next(row for row in ligand_rows if row["candidate_id"] == "9.2")
    note = [
        "",
        "## Update Note: Added 9.2 / L004",
        "",
        "- Original Gate 1 summary covered 5 candidates.",
        "- This updated version adds 9.2 / L004.",
        f"- 9.2 / L004 priority class: {nine_two['WP2_Gate1_priority_class']}",
        f"- 9.2 / L004 cholesterol response label: {nine_two['cholesterol_response_label']}",
        f"- 9.2 / L004 rationale: {nine_two['rationale']}",
        f"- Priority A after adding 9.2: {', '.join(priority_a) if priority_a else 'none'}",
        "- Existing Priority A ranking for 19.1 / L003 and 17.2 / L010 is unchanged if 9.2 is not also assigned Priority A.",
    ]
    with BASE.REPORT_PATH.open("a", encoding="utf-8") as handle:
        handle.write("\n".join(note) + "\n")


BASE.write_pymol_script = write_pymol_script
BASE.write_report = write_report


def main() -> None:
    BASE.main()


if __name__ == "__main__":
    main()
