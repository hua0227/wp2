from __future__ import annotations

import json
from pathlib import Path


SUMMARY_PATH = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_100ps_summary.md")
JSON_20 = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_20chol\md100ps_R2\reports\md100ps_20chol_R2_result.json")
JSON_40 = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_40chol\md100ps_R2\reports\md100ps_40chol_R2_result.json")


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    data20 = load(JSON_20)
    data40 = load(JSON_40)
    final20 = data20.get("final_row", {})
    final40 = data40.get("final_row", {})
    lines = [
        "# 19.1 / L003 100 ps Restrained MD Pilot Summary",
        "",
        "## Completion",
        "",
        f"- 20chol completed 100 ps: {data20['completed_100ps']}",
        f"- 40chol completed 100 ps: {data40['completed_100ps']}",
        f"- Both continued from checkpoint: {data20['continued_from_checkpoint'] and data40['continued_from_checkpoint']}",
        "",
        "## 20chol",
        "",
        f"- Final temperature: {final20.get('temperature_K', 'nan')} K",
        f"- Final potential energy: {final20.get('potential_energy_kJ_mol', 'nan')} kJ/mol",
        f"- NaN detected: {data20['has_nan']}",
        f"- Final TYR13 / VAL17 / SER20: {final20.get('ligand_distance_TYR13_A', 'nan')} / {final20.get('ligand_distance_VAL17_A', 'nan')} / {final20.get('ligand_distance_SER20_A', 'nan')} A",
        f"- Near chains count: {final20.get('ligand_near_chains_count', 'nan')}",
        f"- Ligand RMSD final/max: {data20.get('ligand_rmsd_final_A', 'nan')} / {data20.get('ligand_rmsd_max_A', 'nan')} A",
        f"- Final severe/close: {final20.get('severe_lt_1A', 'nan')}/{final20.get('close_lt_1p5A', 'nan')}",
        f"- Recommend 20chol vs 40chol comparison analysis: {data20['recommend_comparison_analysis']}",
        "",
        "## 40chol",
        "",
        f"- Final temperature: {final40.get('temperature_K', 'nan')} K",
        f"- Final potential energy: {final40.get('potential_energy_kJ_mol', 'nan')} kJ/mol",
        f"- NaN detected: {data40['has_nan']}",
        f"- Final TYR13 / VAL17 / SER20: {final40.get('ligand_distance_TYR13_A', 'nan')} / {final40.get('ligand_distance_VAL17_A', 'nan')} / {final40.get('ligand_distance_SER20_A', 'nan')} A",
        f"- Near chains count: {final40.get('ligand_near_chains_count', 'nan')}",
        f"- Ligand RMSD final/max: {data40.get('ligand_rmsd_final_A', 'nan')} / {data40.get('ligand_rmsd_max_A', 'nan')} A",
        f"- Final severe/close: {final40.get('severe_lt_1A', 'nan')}/{final40.get('close_lt_1p5A', 'nan')}",
        f"- Recommend 20chol vs 40chol comparison analysis: {data40['recommend_comparison_analysis']}",
        "",
        "## Interpretation",
        "",
        f"- Both systems recommend 19.1 20chol vs 40chol comparison analysis: {data20['recommend_comparison_analysis'] and data40['recommend_comparison_analysis']}",
        "- Initial comparison with 14.2 suggests 19.1 continues to keep TYR13, VAL17, and SER20 in a more balanced 3-5 A contact band across both cholesterol conditions, rather than showing the clearer VAL17/SER20 weakening seen for 14.2 in 40chol.",
        "- This remains a non-production restrained MD pilot assessment, not a production MD conclusion.",
    ]
    SUMMARY_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"summary_path: {SUMMARY_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
