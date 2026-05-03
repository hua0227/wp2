from pathlib import Path
import importlib.util


SCRIPT_PATH = Path(__file__).with_name("compare_19.1_20chol_40chol.py")


def load_module():
    spec = importlib.util.spec_from_file_location("compare_19_1_20_40", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_normalize_monitor_row_supports_19_1_100ps_headers():
    mod = load_module()
    row = {
        "time_ps": "100",
        "temperature_K": "312.407",
        "potential_energy_kJ_mol": "-456382.600128",
        "total_energy_kJ_mol": "-310473.537628",
        "ligand_distance_TYR13_A": "3.7754",
        "ligand_distance_VAL17_A": "3.19104",
        "ligand_distance_SER20_A": "3.65355",
        "ligand_min_distance_to_key_residues_A": "3.19104",
        "ligand_near_chains_count": "2",
        "ligand_heavy_rmsd_vs_100ps_input_A": "1.88172",
        "severe_lt_1A": "0",
        "close_lt_1p5A": "0",
        "status": "ok",
    }
    normalized = mod.normalize_monitor_row("20chol", row)
    assert normalized["temperature_K"] == 312.407
    assert normalized["TYR13_distance_A"] == 3.7754
    assert normalized["VAL17_distance_A"] == 3.19104
    assert normalized["SER20_distance_A"] == 3.65355
    assert normalized["key_min_distance_A"] == 3.19104


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
        "TYR13_distance_min",
        "TYR13_distance_max",
        "TYR13_distance_mean",
        "VAL17_distance_initial",
        "VAL17_distance_final",
        "VAL17_distance_min",
        "VAL17_distance_max",
        "VAL17_distance_mean",
        "SER20_distance_initial",
        "SER20_distance_final",
        "SER20_distance_min",
        "SER20_distance_max",
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


def test_contact_occupancy_summary_headers_are_stable():
    mod = load_module()
    assert mod.OCCUPANCY_FIELDS == [
        "system",
        "per_residue_frame_count",
        "TYR13_occ_lt4A",
        "TYR13_occ_lt5A",
        "TYR13_occ_lt6A",
        "TYR13_occ_lt8A",
        "VAL17_occ_lt4A",
        "VAL17_occ_lt5A",
        "VAL17_occ_lt6A",
        "VAL17_occ_lt8A",
        "SER20_occ_lt4A",
        "SER20_occ_lt5A",
        "SER20_occ_lt6A",
        "SER20_occ_lt8A",
        "key_min_occ_lt4A",
        "key_min_occ_lt5A",
        "key_min_occ_lt6A",
        "near_chains_eq2_fraction",
        "ligand_rmsd_mean",
        "ligand_rmsd_median",
        "ligand_rmsd_max",
        "ligand_rmsd_p95",
    ]
