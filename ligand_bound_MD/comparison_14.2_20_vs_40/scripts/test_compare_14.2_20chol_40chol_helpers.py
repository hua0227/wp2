from pathlib import Path
import importlib.util
import math


SCRIPT_PATH = Path(__file__).with_name("compare_14.2_20chol_40chol.py")


def load_module():
    spec = importlib.util.spec_from_file_location("compare_14_2_20_40", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_normalize_monitor_row_supports_20chol_and_40chol_headers():
    mod = load_module()
    row20 = {
        "time_ps": "100",
        "actual_temperature_K": "310.952",
        "potential_energy_kJ_mol": "-453165",
        "total_energy_kJ_mol": "-312000",
        "ligand_key_min_distance_A": "3.6668",
        "ligand_near_chains_count": "2",
        "ligand_heavy_rmsd_vs_100ps_input_A": "1.19815",
        "severe_lt_1A": "0",
        "close_lt_1p5A": "0",
        "status": "ok",
    }
    row40 = {
        "time_ps": "100",
        "temperature_K": "311.505",
        "potential_energy_kJ_mol": "-447460.465845",
        "total_energy_kJ_mol": "-305924.840845",
        "ligand_min_distance_to_key_residues_A": "4.47124",
        "ligand_distance_TYR13_A": "4.47124",
        "ligand_distance_VAL17_A": "8.74193",
        "ligand_distance_SER20_A": "7.98056",
        "ligand_near_chains_count": "2",
        "ligand_heavy_rmsd_vs_100ps_input_A": "1.55302",
        "severe_lt_1A": "0",
        "close_lt_1p5A": "0",
        "status": "ok",
    }

    n20 = mod.normalize_monitor_row("20chol", row20)
    n40 = mod.normalize_monitor_row("40chol", row40)

    assert n20["temperature_K"] == 310.952
    assert n20["key_min_distance_A"] == 3.6668
    assert math.isnan(n20["TYR13_distance_A"])
    assert n40["temperature_K"] == 311.505
    assert n40["VAL17_distance_A"] == 8.74193


def test_required_table_headers_are_stable():
    mod = load_module()
    assert mod.SUMMARY_FIELDS == [
        "system",
        "completed_ps",
        "temperature_min",
        "temperature_max",
        "temperature_mean",
        "temperature_final",
        "potential_energy_initial",
        "potential_energy_final",
        "potential_energy_mean",
        "ligand_rmsd_initial",
        "ligand_rmsd_final",
        "ligand_rmsd_max",
        "key_min_distance_initial",
        "key_min_distance_final",
        "key_min_distance_min",
        "key_min_distance_max",
        "key_min_distance_mean",
        "TYR13_distance_initial",
        "TYR13_distance_final",
        "TYR13_distance_mean",
        "VAL17_distance_initial",
        "VAL17_distance_final",
        "VAL17_distance_mean",
        "SER20_distance_initial",
        "SER20_distance_final",
        "SER20_distance_mean",
        "near_chains_all_2",
        "severe_clash_ever",
        "close_clash_ever",
        "NaN_detected",
    ]
    assert mod.COMBINED_FIELDS == [
        "time_ps",
        "system",
        "temperature_K",
        "potential_energy_kJ_mol",
        "ligand_rmsd_A",
        "key_min_distance_A",
        "TYR13_distance_A",
        "VAL17_distance_A",
        "SER20_distance_A",
        "near_chains",
        "severe_clashes",
        "close_contacts",
        "status",
    ]


def test_merge_recomputed_distances_does_not_fabricate_missing_times():
    mod = load_module()
    rows = [
        {"time_ps": 0.0, "system": "20chol", "TYR13_distance_A": math.nan},
        {"time_ps": 1.0, "system": "20chol", "TYR13_distance_A": math.nan},
        {"time_ps": 2.0, "system": "20chol", "TYR13_distance_A": math.nan},
    ]
    recomputed = {2.0: {"TYR13_distance_A": 3.1, "VAL17_distance_A": 3.8, "SER20_distance_A": 4.0}}

    mod.merge_recomputed_per_residue(rows, recomputed)

    assert math.isnan(rows[0]["TYR13_distance_A"])
    assert math.isnan(rows[1]["TYR13_distance_A"])
    assert rows[2]["TYR13_distance_A"] == 3.1
