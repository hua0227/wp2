import csv
import importlib.util
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(r"C:\TRKB_WP2\docking")
SCRIPT = ROOT / "scripts" / "make_integrated_scoring_table.py"


def load_module():
    spec = importlib.util.spec_from_file_location("make_integrated_scoring_table", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class IntegratedScoringTableTests(unittest.TestCase):
    def test_anchor_panel_matches_requested_references(self):
        module = load_module()

        self.assertEqual(module.ANCHOR_REFS, {"FLX_R", "FLX_S", "IMI", "HNK_RR"})

    def test_outputs_include_requested_summary_columns(self):
        subprocess.run(
            [sys.executable, str(SCRIPT)],
            cwd=str(ROOT),
            check=True,
            capture_output=True,
            text=True,
        )

        with open(ROOT / "integrated_scoring_by_system.csv", newline="", encoding="utf-8-sig") as f:
            system_fields = set(csv.DictReader(f).fieldnames or [])
        with open(ROOT / "integrated_candidate_summary.csv", newline="", encoding="utf-8-sig") as f:
            candidate_fields = set(csv.DictReader(f).fieldnames or [])

        self.assertTrue(
            {
                "group",
                "ligand_id",
                "system",
                "mode",
                "affinity",
                "min_distance_A",
                "contacts_4A",
                "contacts_6A",
                "chain_contacts",
                "system_score",
            }.issubset(system_fields)
        )
        self.assertTrue(
            {
                "final_score_100",
                "priority_class",
                "score_20chol",
                "score_40chol",
                "mode_20chol",
                "mode_40chol",
                "best_system",
                "best_mode",
            }.issubset(candidate_fields)
        )


if __name__ == "__main__":
    unittest.main()
