import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(r"C:\TRKB_WP2\ligand_bound_MD\scripts\sync_top10_best_poses.py")


def load_module():
    spec = importlib.util.spec_from_file_location("sync_top10_best_poses", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SyncTop10BestPosesTests(unittest.TestCase):
    def test_extract_pdbqt_mode_returns_requested_model_only(self):
        module = load_module()
        text = "\n".join(
            [
                "MODEL 1",
                "REMARK VINA RESULT: -1.0",
                "ATOM      1  C   UNL     1       1.000   1.000   1.000  0.00  0.00    +0.000 C",
                "ENDMDL",
                "MODEL 2",
                "REMARK VINA RESULT: -2.0",
                "ATOM      1  N   UNL     1       2.000   2.000   2.000  0.00  0.00    +0.000 N",
                "ENDMDL",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "multi.pdbqt"
            source.write_text(text, encoding="utf-8")

            extracted = module.extract_pdbqt_mode(source, 2)

        self.assertIn("MODEL 2", extracted)
        self.assertIn("REMARK VINA RESULT: -2.0", extracted)
        self.assertIn(" N   UNL ", extracted)
        self.assertNotIn("MODEL 1", extracted)
        self.assertTrue(extracted.endswith("ENDMDL\n"))


if __name__ == "__main__":
    unittest.main()
