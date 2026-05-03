from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT_PATH = Path(
    r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\batch_2.3_6.2_17.2\scripts\batch_compare_candidates.py"
)


def load_module():
    spec = importlib.util.spec_from_file_location("batch_compare_candidates", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_discover_inputs_prefers_real_md100ps_outputs():
    mod = load_module()
    inputs = mod.discover_inputs_for_ligand("2.3", "L006")
    assert set(inputs) == {"20chol", "40chol"}
    assert inputs["20chol"]["monitor_csv"].name == "md100ps_20chol_R2_monitor.csv"
    assert inputs["20chol"]["result_json"].name == "md100ps_20chol_R2_result.json"
    assert inputs["40chol"]["final_pdb"].name == "TRKB_40chol_L006_2.3_100ps_R2_final.pdb"
    assert inputs["40chol"]["summary_md"].name == "md100ps_40chol_R2_summary.md"


def test_contact_mode_labels_cover_balanced_shifted_and_weak_cases():
    mod = load_module()
    assert (
        mod.classify_contact_mode(
            {"TYR13_occ_lt5A": 0.91, "VAL17_occ_lt5A": 0.88, "SER20_occ_lt5A": 0.84}
        )
        == "balanced_three_residue_contact"
    )
    assert (
        mod.classify_contact_mode(
            {"TYR13_occ_lt5A": 0.18, "VAL17_occ_lt5A": 0.92, "SER20_occ_lt5A": 0.89}
        )
        == "VAL17_SER20_dominant"
    )
    assert (
        mod.classify_contact_mode(
            {"TYR13_occ_lt5A": 0.95, "VAL17_occ_lt5A": 0.12, "SER20_occ_lt5A": 0.55}
        )
        == "TYR13_shifted"
    )
    assert (
        mod.classify_contact_mode(
            {"TYR13_occ_lt5A": 0.61, "VAL17_occ_lt5A": 0.42, "SER20_occ_lt5A": 0.78}
        )
        == "weak_or_mixed_contact"
    )


def test_cholesterol_response_labels_cover_defined_gate_rules():
    mod = load_module()
    balanced = mod.classify_cholesterol_response(
        {
            "20chol": {"contact_mode_label": "balanced_three_residue_contact", "near_chains_all_2": True},
            "40chol": {"contact_mode_label": "balanced_three_residue_contact", "near_chains_all_2": True},
        }
    )
    assert balanced == "cholesterol_robust_balanced"

    shifted = mod.classify_cholesterol_response(
        {
            "20chol": {"contact_mode_label": "VAL17_SER20_dominant", "near_chains_all_2": True},
            "40chol": {"contact_mode_label": "TYR13_shifted", "near_chains_all_2": True},
        }
    )
    assert shifted == "cholesterol_sensitive_contact_shift"

    weaker = mod.classify_cholesterol_response(
        {
            "20chol": {
                "contact_mode_label": "balanced_three_residue_contact",
                "near_chains_all_2": True,
                "TYR13_occ_lt5A": 0.95,
                "VAL17_occ_lt5A": 0.96,
                "SER20_occ_lt5A": 0.93,
            },
            "40chol": {
                "contact_mode_label": "weak_or_mixed_contact",
                "near_chains_all_2": True,
                "TYR13_occ_lt5A": 0.94,
                "VAL17_occ_lt5A": 0.37,
                "SER20_occ_lt5A": 0.31,
            },
        }
    )
    assert weaker == "interface_retained_but_weaker_contact"

    inconclusive = mod.classify_cholesterol_response(
        {
            "20chol": {"contact_mode_label": "weak_or_mixed_contact", "near_chains_all_2": False},
            "40chol": {"contact_mode_label": "weak_or_mixed_contact", "near_chains_all_2": True},
        }
    )
    assert inconclusive == "weak_or_inconclusive"
