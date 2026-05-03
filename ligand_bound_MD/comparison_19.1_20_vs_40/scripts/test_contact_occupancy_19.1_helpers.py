from pathlib import Path
import importlib.util


SCRIPT_PATH = Path(__file__).with_name("compare_19.1_20chol_40chol.py")


def load_module():
    spec = importlib.util.spec_from_file_location("compare_19_1_20_40_occ", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_occupancy_uses_only_finite_values_in_denominator():
    mod = load_module()
    values = [3.5, 4.5, float("nan"), 7.0]
    result = mod.occupancy(values, 5.0)
    assert result["hits"] == 2
    assert result["total"] == 3
    assert abs(result["fraction"] - (2.0 / 3.0)) < 1e-12


def test_system_summary_uses_real_per_residue_data_without_fabrication():
    mod = load_module()
    rows = [
        {
            "system": "20chol",
            "TYR13_distance_A": 3.9,
            "VAL17_distance_A": 4.4,
            "SER20_distance_A": 5.1,
            "key_min_distance_A": 3.9,
            "near_chains": 2,
            "ligand_rmsd_A": 1.0,
        },
        {
            "system": "20chol",
            "TYR13_distance_A": 4.8,
            "VAL17_distance_A": 3.7,
            "SER20_distance_A": 4.9,
            "key_min_distance_A": 3.7,
            "near_chains": 2,
            "ligand_rmsd_A": 1.2,
        },
        {
            "system": "20chol",
            "TYR13_distance_A": 5.3,
            "VAL17_distance_A": 4.2,
            "SER20_distance_A": 3.8,
            "key_min_distance_A": 3.8,
            "near_chains": 1,
            "ligand_rmsd_A": 1.4,
        },
    ]
    summary = mod.summarize_occupancy_system("20chol", rows)
    assert summary["per_residue_frame_count"] == 3
    assert abs(summary["TYR13_occ_lt5A"] - (2.0 / 3.0)) < 1e-12
    assert abs(summary["VAL17_occ_lt4A"] - (1.0 / 3.0)) < 1e-12
    assert abs(summary["SER20_occ_lt6A"] - 1.0) < 1e-12
    assert abs(summary["key_min_occ_lt4A"] - 1.0) < 1e-12
    assert abs(summary["near_chains_eq2_fraction"] - (2.0 / 3.0)) < 1e-12
    assert abs(summary["ligand_rmsd_median"] - 1.2) < 1e-12
