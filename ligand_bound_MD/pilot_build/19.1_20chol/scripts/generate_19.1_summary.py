from __future__ import annotations

import json
from pathlib import Path


SUMMARY_PATH = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_summary.md")
JSON_20 = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_20chol\reports\pilot_gate_summary.json")
JSON_40 = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_40chol\reports\pilot_gate_summary.json")


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    data20 = load(JSON_20)
    data40 = load(JSON_40)
    lines = [
        "# 19.1 / L003 Non-production Restrained MD Pilot Summary",
        "",
        "## Why 19.1 / L003",
        "",
        "19.1 / L003 was selected as the second ligand-bound non-production pilot after the validated 14.2 workflow, so the same gate structure could be reused while forcing a fresh topology build, clash diagnosis, cleanup decision, R2 minimization, very-short probe, and 10 ps restrained MD check.",
        "",
        "## 20chol Pilot",
        "",
        f"- PSF/PDB generated: {data20['psf_pdb_generated']}",
        f"- Read test succeeded: {data20['read_test_succeeded']}",
        f"- Cleanup executed: {data20['cleanup_executed']}",
        f"- Deleted molecules: {data20['deleted_molecules']}",
        f"- Minimization succeeded: {data20['minimization_succeeded']}",
        f"- Very-short probe passed: {data20['very_short_probe_passed']}",
        f"- 10 ps completed: {data20['completed_10ps']}",
        f"- Final TYR13/VAL17/SER20 distances: {data20.get('ligand_distance_TYR13_A', 'nan')} / {data20.get('ligand_distance_VAL17_A', 'nan')} / {data20.get('ligand_distance_SER20_A', 'nan')}",
        f"- Near chains count: {data20.get('ligand_near_chains_count', 'nan')}",
        f"- Ligand RMSD: {data20.get('ligand_rmsd_A', 'nan')}",
        f"- Recommend 50 ps: {data20['recommend_50ps']}",
        "",
        "## 40chol Pilot",
        "",
        f"- PSF/PDB generated: {data40['psf_pdb_generated']}",
        f"- Read test succeeded: {data40['read_test_succeeded']}",
        f"- Cleanup executed: {data40['cleanup_executed']}",
        f"- Deleted molecules: {data40['deleted_molecules']}",
        f"- Minimization succeeded: {data40['minimization_succeeded']}",
        f"- Very-short probe passed: {data40['very_short_probe_passed']}",
        f"- 10 ps completed: {data40['completed_10ps']}",
        f"- Final TYR13/VAL17/SER20 distances: {data40.get('ligand_distance_TYR13_A', 'nan')} / {data40.get('ligand_distance_VAL17_A', 'nan')} / {data40.get('ligand_distance_SER20_A', 'nan')}",
        f"- Near chains count: {data40.get('ligand_near_chains_count', 'nan')}",
        f"- Ligand RMSD: {data40.get('ligand_rmsd_A', 'nan')}",
        f"- Recommend 50 ps: {data40['recommend_50ps']}",
        "",
        "## Comparison",
        "",
        f"- Both systems completed 10 ps: {data20['completed_10ps'] and data40['completed_10ps']}",
        f"- Preliminary contact pattern suggests 20chol distances: {data20.get('ligand_distance_TYR13_A', 'nan')} / {data20.get('ligand_distance_VAL17_A', 'nan')} / {data20.get('ligand_distance_SER20_A', 'nan')} A and 40chol distances: {data40.get('ligand_distance_TYR13_A', 'nan')} / {data40.get('ligand_distance_VAL17_A', 'nan')} / {data40.get('ligand_distance_SER20_A', 'nan')} A.",
        f"- Recommend continuing to 50 ps: 20chol={data20['recommend_50ps']}, 40chol={data40['recommend_50ps']}.",
        "",
        "This is a non-production restrained MD pilot summary, not a production MD conclusion.",
    ]
    SUMMARY_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"summary_path: {SUMMARY_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
