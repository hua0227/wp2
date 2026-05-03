from pathlib import Path
import importlib.util


SCRIPT_PATH = Path(__file__).with_name("contact_occupancy_analysis.py")


def load_module():
    spec = importlib.util.spec_from_file_location("contact_occupancy_analysis", SCRIPT_PATH)
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


def test_summary_fields_cover_requested_outputs():
    mod = load_module()
    assert mod.SUMMARY_FIELDS == [
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


def test_system_summary_respects_sparse_per_residue_data():
    mod = load_module()
    rows = [
        {
            "system": "20chol",
            "TYR13_distance_A": float("nan"),
            "VAL17_distance_A": float("nan"),
            "SER20_distance_A": float("nan"),
            "key_min_distance_A": 3.8,
            "near_chains": 2,
            "ligand_rmsd_A": 1.0,
        },
        {
            "system": "20chol",
            "TYR13_distance_A": 6.2,
            "VAL17_distance_A": 3.6,
            "SER20_distance_A": 4.4,
            "key_min_distance_A": 3.6,
            "near_chains": 2,
            "ligand_rmsd_A": 1.2,
        },
        {
            "system": "20chol",
            "TYR13_distance_A": 7.4,
            "VAL17_distance_A": 4.1,
            "SER20_distance_A": 3.9,
            "key_min_distance_A": 4.2,
            "near_chains": 1,
            "ligand_rmsd_A": 1.4,
        },
    ]
    summary = mod.summarize_system("20chol", rows)

    assert summary["per_residue_frame_count"] == 2
    assert abs(summary["VAL17_occ_lt4A"] - 0.5) < 1e-12
    assert abs(summary["SER20_occ_lt5A"] - 1.0) < 1e-12
    assert abs(summary["key_min_occ_lt4A"] - (2.0 / 3.0)) < 1e-12
    assert abs(summary["near_chains_eq2_fraction"] - (2.0 / 3.0)) < 1e-12
    assert abs(summary["ligand_rmsd_median"] - 1.2) < 1e-12
