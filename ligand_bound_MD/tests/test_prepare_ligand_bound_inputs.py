import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(r"C:\TRKB_WP2\ligand_bound_MD\scripts\prepare_ligand_bound_inputs.py")


def load_module():
    spec = importlib.util.spec_from_file_location("prepare_ligand_bound_inputs", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PrepareLigandBoundInputsTests(unittest.TestCase):
    def test_find_exact_ligand_file_does_not_match_neighbor_ids(self):
        module = load_module()

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "18.1.sdf").write_text("wrong\n", encoding="utf-8")
            (root / "8.1.sdf").write_text("right\n", encoding="utf-8")

            found = module.find_exact_ligand_file("8.1", [root], ".sdf")

        self.assertEqual(found.name, "8.1.sdf")


if __name__ == "__main__":
    unittest.main()
