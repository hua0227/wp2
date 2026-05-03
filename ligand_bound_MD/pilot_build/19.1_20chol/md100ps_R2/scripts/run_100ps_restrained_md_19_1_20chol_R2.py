from __future__ import annotations

import sys
from pathlib import Path


COMMON_DIR = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\19.1_20chol\md100ps_R2\scripts")
if str(COMMON_DIR) not in sys.path:
    sys.path.insert(0, str(COMMON_DIR))

from md100ps_common_19_1 import make_config, run_100ps


def main() -> int:
    config = make_config("20chol")
    result = run_100ps(config)
    final = result.get("final_row", {})
    print("For 20chol:")
    print(f"- completed 100 ps: {'YES' if result['completed_100ps'] else 'NO'}")
    print(f"- continued from checkpoint: {'YES' if result['continued_from_checkpoint'] else 'NO'}")
    print(f"- final temperature: {final.get('temperature_K', 'nan')}")
    print(f"- final potential energy: {final.get('potential_energy_kJ_mol', 'nan')}")
    print(f"- NaN: {'YES' if result['has_nan'] else 'NO'}")
    print(
        "- final TYR13 / VAL17 / SER20 distances: "
        f"{final.get('ligand_distance_TYR13_A', 'nan')} / "
        f"{final.get('ligand_distance_VAL17_A', 'nan')} / "
        f"{final.get('ligand_distance_SER20_A', 'nan')}"
    )
    print(f"- near chains count: {final.get('ligand_near_chains_count', 'nan')}")
    print(f"- ligand RMSD final/max: {result.get('ligand_rmsd_final_A', float('nan'))} / {result.get('ligand_rmsd_max_A', float('nan'))}")
    print(f"- final severe/close contacts: {final.get('severe_lt_1A', 'nan')}/{final.get('close_lt_1p5A', 'nan')}")
    print(f"- recommend comparison analysis: {'YES' if result['recommend_comparison_analysis'] else 'NO'}")
    return 0 if result["completed_100ps"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
