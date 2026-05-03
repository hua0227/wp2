from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT_PATH = Path(
    r"C:\TRKB_WP2\ligand_bound_MD\WP2_Gate1_summary\scripts\generate_WP2_Gate1_summary.py"
)


def load_module():
    spec = importlib.util.spec_from_file_location("generate_WP2_Gate1_summary", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_contact_mode_rules_match_gate1_definitions():
    mod = load_module()
    assert (
        mod.classify_contact_mode({"TYR13_occ_lt5A": 0.95, "VAL17_occ_lt5A": 0.92, "SER20_occ_lt5A": 0.88})
        == "balanced_three_residue_contact"
    )
    assert (
        mod.classify_contact_mode({"TYR13_occ_lt5A": 0.20, "VAL17_occ_lt5A": 0.91, "SER20_occ_lt5A": 0.93})
        == "VAL17_SER20_dominant"
    )
    assert (
        mod.classify_contact_mode({"TYR13_occ_lt5A": 0.97, "VAL17_occ_lt5A": 0.69, "SER20_occ_lt5A": 0.01})
        == "TYR13_shifted"
    )
    assert (
        mod.classify_contact_mode({"TYR13_occ_lt5A": 0.74, "VAL17_occ_lt5A": 0.58, "SER20_occ_lt5A": 0.61})
        == "weak_or_mixed_contact"
    )


def test_cholesterol_response_rules_cover_expected_ligand_patterns():
    mod = load_module()
    assert (
        mod.classify_cholesterol_response(
            {
                "20chol": {
                    "contact_mode_label": "balanced_three_residue_contact",
                    "near_chains_all_2": True,
                    "TYR13_occ_lt5A": 1.0,
                    "VAL17_occ_lt5A": 1.0,
                    "SER20_occ_lt5A": 1.0,
                },
                "40chol": {
                    "contact_mode_label": "balanced_three_residue_contact",
                    "near_chains_all_2": True,
                    "TYR13_occ_lt5A": 0.97,
                    "VAL17_occ_lt5A": 0.99,
                    "SER20_occ_lt5A": 1.0,
                },
            }
        )
        == "cholesterol_robust_balanced"
    )
    assert (
        mod.classify_cholesterol_response(
            {
                "20chol": {
                    "contact_mode_label": "VAL17_SER20_dominant",
                    "near_chains_all_2": True,
                    "TYR13_occ_lt5A": 0.0,
                    "VAL17_occ_lt5A": 1.0,
                    "SER20_occ_lt5A": 0.98,
                },
                "40chol": {
                    "contact_mode_label": "TYR13_shifted",
                    "near_chains_all_2": True,
                    "TYR13_occ_lt5A": 0.96,
                    "VAL17_occ_lt5A": 0.0,
                    "SER20_occ_lt5A": 0.0,
                },
            }
        )
        == "cholesterol_sensitive_contact_shift"
    )
    assert (
        mod.classify_cholesterol_response(
            {
                "20chol": {
                    "contact_mode_label": "TYR13_shifted",
                    "near_chains_all_2": True,
                    "TYR13_occ_lt5A": 1.0,
                    "VAL17_occ_lt5A": 1.0,
                    "SER20_occ_lt5A": 0.16,
                },
                "40chol": {
                    "contact_mode_label": "TYR13_shifted",
                    "near_chains_all_2": True,
                    "TYR13_occ_lt5A": 1.0,
                    "VAL17_occ_lt5A": 0.69,
                    "SER20_occ_lt5A": 0.01,
                },
            }
        )
        == "TYR13_dominant_interface_retention"
    )
    assert (
        mod.classify_cholesterol_response(
            {
                "20chol": {
                    "contact_mode_label": "balanced_three_residue_contact",
                    "near_chains_all_2": True,
                    "TYR13_occ_lt5A": 0.99,
                    "VAL17_occ_lt5A": 1.0,
                    "SER20_occ_lt5A": 1.0,
                },
                "40chol": {
                    "contact_mode_label": "weak_or_mixed_contact",
                    "near_chains_all_2": True,
                    "TYR13_occ_lt5A": 1.0,
                    "VAL17_occ_lt5A": 0.98,
                    "SER20_occ_lt5A": 0.63,
                },
            }
        )
        == "interface_retained_but_weaker_contact"
    )


def test_priority_class_rules_promote_balanced_and_penalize_weaker_contacts():
    mod = load_module()
    priority_a = mod.classify_priority(
        {
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
    )
    assert priority_a == "Priority A"

    priority_b = mod.classify_priority(
        {
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
    )
    assert priority_b == "Priority B"

    priority_c = mod.classify_priority(
        {
            "candidate_id": "6.2",
            "20chol": {
                "completed_100ps": True,
                "NaN_detected": False,
                "near_chains_all_2": True,
                "severe_clash_ever": False,
                "close_clash_ever": False,
                "contact_mode_label": "weak_or_mixed_contact",
            },
            "40chol": {
                "completed_100ps": True,
                "NaN_detected": False,
                "near_chains_all_2": True,
                "severe_clash_ever": False,
                "close_clash_ever": False,
                "contact_mode_label": "TYR13_shifted",
            },
            "cholesterol_response_label": "interface_retained_but_weaker_contact",
        }
    )
    assert priority_c == "Priority C"
