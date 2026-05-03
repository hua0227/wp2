from pathlib import Path
import importlib.util
import sys


SCRIPT_PATH = Path(__file__).with_name("md100ps_common_19_1.py")


def load_module():
    spec = importlib.util.spec_from_file_location("md100ps_common_19_1", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_100ps_budget_and_output_cadence_are_bounded():
    mod = load_module()
    assert mod.TOTAL_STEPS_100PS == 200000
    assert mod.TIMESTEP_100PS == 0.0005
    assert mod.REPORT_INTERVAL_STEPS == 2000
    assert mod.DCD_INTERVAL_STEPS == 4000


def test_checkpoint_preferred_and_fallback_is_explicit():
    mod = load_module()
    text = mod.build_continuation_note(True, False, "")
    assert "continued from checkpoint" in text.lower()
    fallback = mod.build_continuation_note(False, True, "checkpoint load failed")
    assert "fallback used" in fallback.lower()
    assert "checkpoint load failed" in fallback


def test_monitor_fields_cover_requested_metrics():
    mod = load_module()
    assert mod.MONITOR_FIELDS == [
        "step",
        "time_ps",
        "temperature_K",
        "potential_energy_kJ_mol",
        "kinetic_energy_kJ_mol",
        "total_energy_kJ_mol",
        "ligand_distance_TYR13_A",
        "ligand_distance_VAL17_A",
        "ligand_distance_SER20_A",
        "ligand_min_distance_to_key_residues_A",
        "ligand_near_chains_count",
        "ligand_heavy_rmsd_vs_100ps_input_A",
        "ligand_heavy_rmsd_vs_original_R2_minimized_input_A",
        "severe_lt_1A",
        "close_lt_1p5A",
        "status",
    ]
