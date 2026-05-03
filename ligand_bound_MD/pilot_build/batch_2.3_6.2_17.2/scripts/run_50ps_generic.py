from __future__ import annotations

import sys

from generic_md50ps_common import LIGAND_ID, RESNAME, make_config, run_50ps


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_50ps_generic.py <20chol|40chol>")
    system_name = sys.argv[1]
    config = make_config(system_name)
    result = run_50ps(config)
    final = result.get("final_row", {})
    print(f"For {LIGAND_ID} / {RESNAME} / {system_name}:")
    print(f"- completed 50 ps: {'YES' if result['completed_50ps'] else 'NO'}")
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
    print(f"- recommend 100 ps: {'YES' if result['recommend_100ps'] else 'NO'}")
    return 0 if result["completed_50ps"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
