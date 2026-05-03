import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(r"C:\TRKB_WP2\ligand_bound_MD\preflight\scripts\preflight_check_ligand_parameters.py")


def load_module():
    spec = importlib.util.spec_from_file_location("preflight_check_ligand_parameters", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PreflightCheckLigandParametersTests(unittest.TestCase):
    def test_parse_mol2_atoms_reads_names_elements_and_count(self):
        module = load_module()
        mol2 = "\n".join(
            [
                "@<TRIPOS>MOLECULE",
                "L001",
                "@<TRIPOS>ATOM",
                "1 C1 0.0 0.0 0.0 C.3 1 L001 0.0",
                "2 N2 1.0 0.0 0.0 N.am 1 L001 0.0",
                "@<TRIPOS>BOND",
                "1 1 2 1",
                "",
            ]
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "lig.mol2"
            path.write_text(mol2, encoding="utf-8")

            atoms = module.parse_mol2_atoms(path)

        self.assertEqual([atom["name"] for atom in atoms], ["C1", "N2"])
        self.assertEqual([atom["element"] for atom in atoms], ["C", "N"])

    def test_named_pose_uses_mol2_names_resname_chain_and_pose_coordinates(self):
        module = load_module()
        mol2_atoms = [
            {"name": "C1", "element": "C"},
            {"name": "N2", "element": "N"},
        ]
        coords = [
            {"x": 1.25, "y": -2.5, "z": 3.75},
            {"x": 4.0, "y": 5.5, "z": -6.25},
        ]

        lines = module.build_named_pose_pdb_lines(mol2_atoms, coords, "L001")

        self.assertEqual(len(lines), 2)
        self.assertTrue(lines[0].startswith("HETATM"))
        self.assertEqual(lines[0][17:21].strip(), "L001")
        self.assertEqual(lines[0][22], "Z")
        self.assertEqual(lines[0][24:28].strip(), "1")
        self.assertIn("C1", lines[0])
        self.assertIn("N2", lines[1])
        self.assertIn("   1.250", lines[0])
        self.assertIn("  -6.250", lines[1])


if __name__ == "__main__":
    unittest.main()
