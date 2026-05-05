from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT_PATH = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\generate_9.2_gate_summary.py")


def load_module():
    spec = importlib.util.spec_from_file_location("generate_9_2_gate_summary", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_collect_stage_rows_reads_20chol_and_40chol_layout():
    mod = load_module()
    rows = mod.collect_stage_rows()
    assert {row["system"] for row in rows} == {"20chol", "40chol"}
    for row in rows:
        assert row["candidate_id"] == "9.2"
        assert row["resname"] == "L004"


def test_stage_label_prefers_highest_completed_stage():
    mod = load_module()
    assert mod.stage_reached({"completed_10ps": True, "completed_50ps": False, "completed_100ps": False}) == "10ps"
    assert mod.stage_reached({"completed_10ps": True, "completed_50ps": True, "completed_100ps": False}) == "50ps"
    assert mod.stage_reached({"completed_10ps": True, "completed_50ps": True, "completed_100ps": True}) == "100ps"
