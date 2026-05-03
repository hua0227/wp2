from __future__ import annotations

from pilot_common_19_1 import LIGAND_ID, RESNAME, make_config, run_pipeline


def main() -> int:
    config = make_config("20chol")
    result = run_pipeline(config)
    print("For 20chol:")
    print(f"- PSF/PDB generated: {'YES' if result['psf_pdb_generated'] else 'NO'}")
    print(f"- read test succeeded: {'YES' if result['read_test_succeeded'] else 'NO'}")
    print(f"- initial energy: {result['initial_energy_kJ_mol']}")
    print(f"- severe/close before cleanup: {result['severe_before_cleanup']}/{result['close_before_cleanup']}")
    print(f"- cleanup executed and deleted molecules: {result['cleanup_executed']} {result['deleted_molecules']}")
    print(f"- minimization final energy: {result.get('minimization_final_energy_kJ_mol', float('nan'))}")
    print(f"- very-short probe passed: {'YES' if result['very_short_probe_passed'] else 'NO'}")
    print(f"- 10 ps completed: {'YES' if result['completed_10ps'] else 'NO'}")
    print(f"- final temperature: {result.get('final_temperature_K', float('nan'))}")
    print(
        "- final distances to TYR13 / VAL17 / SER20: "
        f"{result.get('ligand_distance_TYR13_A', float('nan'))} / "
        f"{result.get('ligand_distance_VAL17_A', float('nan'))} / "
        f"{result.get('ligand_distance_SER20_A', float('nan'))}"
    )
    print(f"- near chains count: {result.get('ligand_near_chains_count', 'nan')}")
    print(f"- ligand RMSD: {result.get('ligand_rmsd_A', float('nan'))}")
    print(f"- severe/close final: {result.get('severe_lt_1A', 'nan')}/{result.get('close_lt_1p5A', 'nan')}")
    print(f"- recommend 50 ps: {'YES' if result['recommend_50ps'] else 'NO'}")
    if result.get("error"):
        print(f"- error: {result['error']}")
    return 0 if result["completed_10ps"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
