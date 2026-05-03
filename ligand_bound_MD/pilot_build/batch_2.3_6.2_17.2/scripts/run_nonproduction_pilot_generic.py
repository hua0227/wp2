from __future__ import annotations

import os
import sys

from generic_pilot_common import LIGAND_ID, RESNAME, make_config, run_pipeline


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_nonproduction_pilot_generic.py <20chol|40chol>")
    system_name = sys.argv[1]
    config = make_config(system_name)
    result = run_pipeline(config)
    print(f"For {LIGAND_ID} / {RESNAME} / {system_name}:")
    print(f"- preflight ok: {'YES' if result.get('preflight_ok') else 'NO'}")
    print(f"- PSF/PDB generated: {'YES' if result.get('psf_pdb_generated') else 'NO'}")
    print(f"- read test succeeded: {'YES' if result.get('read_test_succeeded') else 'NO'}")
    print(f"- cleanup executed: {'YES' if result.get('cleanup_executed') else 'NO'}")
    print(f"- minimization succeeded: {'YES' if result.get('minimization_succeeded') else 'NO'}")
    print(f"- very-short probe passed: {'YES' if result.get('very_short_probe_passed') else 'NO'}")
    print(f"- 10 ps completed: {'YES' if result.get('completed_10ps') else 'NO'}")
    print(f"- recommend 50 ps: {'YES' if result.get('recommend_50ps') else 'NO'}")
    if result.get("error"):
        print(f"- error: {result['error']}")
    return 0 if result.get("completed_10ps") else 1


if __name__ == "__main__":
    raise SystemExit(main())
