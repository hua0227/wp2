from __future__ import annotations

import importlib.util
from pathlib import Path

MODULE_PATH = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\batch_2.3_6.2_17.2\scripts\batch_compare_candidates.py")
SPEC = importlib.util.spec_from_file_location('batch_compare_candidates', MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)

if __name__ == '__main__':
    MODULE.run_single_ligand('9.2', 'L004')
