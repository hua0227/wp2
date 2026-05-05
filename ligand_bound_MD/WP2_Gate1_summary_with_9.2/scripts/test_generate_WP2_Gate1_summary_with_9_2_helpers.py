from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT_PATH = Path(
    r"C:\TRKB_WP2\ligand_bound_MD\WP2_Gate1_summary_with_9.2\scripts\generate_WP2_Gate1_summary_with_9_2.py"
)


def load_module():
    spec = importlib.util.spec_from_file_location("generate_WP2_Gate1_summary_with_9_2", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_updated_gate1_summary_uses_distinct_output_dir_and_includes_9_2():
    mod = load_module()
    assert str(mod.OUT_DIR).endswith("WP2_Gate1_summary_with_9.2")
    assert mod.CANDIDATES["9.2"] == "L004"
    assert len(mod.CANDIDATES) == 6


def test_updated_gate1_classification_keeps_expected_priority_examples():
    mod = load_module()
    ligand_19_1 = {
        "candidate_id": "19.1",
        "20chol": {
            "completed_100ps": True,
            "NaN_detected": False,
            "near_chains_all_2": True,
            "severe_clash_ever": False,
            "close_clash_ever": False,
            "contact_mode_label": "balanced_three_residue_contact",
            "TYR13_occ_lt5A": 1.0,
            "VAL17_occ_lt5A": 1.0,
            "SER20_occ_lt5A": 1.0,
        },
        "40chol": {
            "completed_100ps": True,
            "NaN_detected": False,
            "near_chains_all_2": True,
            "severe_clash_ever": False,
            "close_clash_ever": False,
            "contact_mode_label": "balanced_three_residue_contact",
            "TYR13_occ_lt5A": 0.97,
            "VAL17_occ_lt5A": 0.99,
            "SER20_occ_lt5A": 1.0,
        },
        "cholesterol_response_label": "cholesterol_robust_balanced",
    }
    assert mod.classify_priority(ligand_19_1) == "Priority A"

    ligand_14_2 = {
        "candidate_id": "14.2",
        "20chol": {
            "completed_100ps": True,
            "NaN_detected": False,
            "near_chains_all_2": True,
            "severe_clash_ever": False,
            "close_clash_ever": False,
            "contact_mode_label": "VAL17_SER20_dominant",
        },
        "40chol": {
            "completed_100ps": True,
            "NaN_detected": False,
            "near_chains_all_2": True,
            "severe_clash_ever": False,
            "close_clash_ever": False,
            "contact_mode_label": "TYR13_shifted",
        },
        "cholesterol_response_label": "cholesterol_sensitive_contact_shift",
    }
    assert mod.classify_priority(ligand_14_2) == "Priority B"
