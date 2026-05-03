import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(r"C:\TRKB_WP2\ligand_bound_MD\scripts\prepare_cgenff_upload_package.py")


def load_module():
    spec = importlib.util.spec_from_file_location("prepare_cgenff_upload_package", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PrepareCgenffUploadPackageTests(unittest.TestCase):
    def test_resname_for_rank_uses_four_character_charmm_names(self):
        module = load_module()

        self.assertEqual(module.resname_for_rank(1), "L001")
        self.assertEqual(module.resname_for_rank(10), "L010")
        self.assertLessEqual(len(module.resname_for_rank(10)), 4)

    def test_check_mol2_sections_detects_atom_and_bond_blocks(self):
        module = load_module()

        mol2_text = "\n".join(
            [
                "@<TRIPOS>MOLECULE",
                "L001",
                "@<TRIPOS>ATOM",
                "1 C1 0.0 0.0 0.0 C.3 1 L001 0.0",
                "@<TRIPOS>BOND",
                "1 1 1 1",
                "",
            ]
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.mol2"
            path.write_text(mol2_text, encoding="utf-8")

            found, has_atoms, has_bonds = module.check_mol2(path)

        self.assertTrue(found)
        self.assertTrue(has_atoms)
        self.assertTrue(has_bonds)


if __name__ == "__main__":
    unittest.main()
