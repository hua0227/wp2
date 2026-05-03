from pathlib import Path
import importlib.util
import os
import sys


COMMON_PATH = Path(__file__).with_name("generic_pilot_common.py")
SUMMARY_PATH = Path(__file__).with_name("generate_batch_status_report.py")


def load_common(ligand_id: str, resname: str):
    os.environ["TRKB_LIGAND_ID"] = ligand_id
    os.environ["TRKB_RESNAME"] = resname
    spec = importlib.util.spec_from_file_location(f"generic_pilot_common_{ligand_id}_{resname}", COMMON_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_summary():
    spec = importlib.util.spec_from_file_location("generate_batch_status_report", SUMMARY_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_make_config_uses_env_selected_ligand_paths():
    mod = load_common("2.3", "L006")
    cfg20 = mod.make_config("20chol")
    assert str(cfg20.pilot_root).endswith(r"pilot_build\2.3_20chol")
    assert str(cfg20.ligand_str).endswith(r"returned_str\L006_2.3.str")
    assert str(cfg20.input_ligand_pdb).endswith(r"ligand_pose_allatom_named\2.3\20chol\L006_2.3_20chol_allatom_named.pdb")

    mod = load_common("17.2", "L010")
    cfg40 = mod.make_config("40chol")
    assert str(cfg40.pilot_root).endswith(r"pilot_build\17.2_40chol")
    assert str(cfg40.ligand_str).endswith(r"returned_str\L010_17.2.str")
    assert str(cfg40.input_ligand_pdb).endswith(r"ligand_pose_allatom_named\17.2\40chol\L010_17.2_40chol_allatom_named.pdb")


def test_batch_status_fields_are_stable():
    mod = load_summary()
    assert mod.STATUS_FIELDS == [
        "ligand_id",
        "resname",
        "system",
        "preflight_ok",
        "psf_pdb_generated",
        "read_test_succeeded",
        "cleanup_executed",
        "minimization_succeeded",
        "very_short_probe_passed",
        "completed_10ps",
        "recommend_50ps",
        "completed_50ps",
        "recommend_100ps",
        "completed_100ps",
        "recommend_comparison_analysis",
        "final_stage",
        "stop_reason",
        "error",
    ]


def test_final_stage_prefers_furthest_successful_gate():
    mod = load_summary()
    row = mod.summarize_system_status(
        "6.2",
        "L008",
        "20chol",
        {
            "gate": {"completed_10ps": True, "recommend_50ps": True, "error": ""},
            "md50": {"completed_50ps": True, "recommend_100ps": True, "stop_reason": "none"},
            "md100": {"completed_100ps": False, "stop_reason": "persistent_severe_clash_gt_0", "recommend_comparison_analysis": False},
        },
    )
    assert row["final_stage"] == "md50ps_failed_during_md100ps"
    assert row["stop_reason"] == "persistent_severe_clash_gt_0"
