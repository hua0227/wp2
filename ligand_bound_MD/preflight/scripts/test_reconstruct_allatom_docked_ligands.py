from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import reconstruct_allatom_docked_ligands as recon


class ReconstructionUtilityTests(unittest.TestCase):
    def test_parse_mol2_atoms_preserves_names_elements_and_heavy_count(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            mol2 = Path(tmp) / "mini.mol2"
            mol2.write_text(
                "\n".join(
                    [
                        "@<TRIPOS>MOLECULE",
                        "L001",
                        " 3 2 0 0 0",
                        "SMALL",
                        "GASTEIGER",
                        "",
                        "@<TRIPOS>ATOM",
                        "1 C1 0.000 0.000 0.000 C.3 1 L001 0.0",
                        "2 CL1 1.000 0.000 0.000 Cl 1 L001 0.0",
                        "3 H1 0.000 1.000 0.000 H 1 L001 0.0",
                        "@<TRIPOS>BOND",
                        "1 1 2 1",
                        "2 1 3 1",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            atoms = recon.parse_mol2_atoms(mol2)

        self.assertEqual([atom.name for atom in atoms], ["C1", "CL1", "H1"])
        self.assertEqual([atom.element for atom in atoms], ["C", "Cl", "H"])
        self.assertEqual(recon.count_heavy_atoms(atoms), 2)

    def test_build_named_pose_lines_uses_cgenff_names_residue_chain_and_end(self) -> None:
        atoms = [
            recon.Mol2Atom(1, "C1", "C", 1.0, 2.0, 3.0, "C.3"),
            recon.Mol2Atom(2, "H1", "H", 1.5, 2.5, 3.5, "H"),
        ]
        coords = [(11.0, 12.0, 13.0), (11.5, 12.5, 13.5)]

        lines = recon.build_named_pose_pdb_lines(atoms, coords, "L001")

        self.assertEqual(len(lines), 3)
        self.assertTrue(lines[0].startswith("HETATM"))
        self.assertIn(" C1 ", lines[0])
        self.assertIn("L001Z", lines[0])
        self.assertEqual(lines[0][17:21], "L001")
        self.assertEqual(lines[0][21], "Z")
        self.assertEqual(lines[-1], "END")

    def test_write_complex_preview_keeps_only_atom_records_and_marks_visual_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            apo = tmp_path / "apo.pdb"
            lig = tmp_path / "lig.pdb"
            out = tmp_path / "preview.pdb"
            apo.write_text(
                "REMARK ignore\nATOM      1  N   MET A   1       0.000   0.000   0.000  1.00  0.00           N  \nEND\n",
                encoding="utf-8",
            )
            lig.write_text(
                "HETATM    1  C1  L001 Z   1       1.000   1.000   1.000  1.00  0.00           C  \nEND\n",
                encoding="utf-8",
            )

            recon.write_complex_preview(apo, lig, out)
            text = out.read_text(encoding="utf-8")

        self.assertIn("visual-only complex preview", text)
        self.assertIn("ATOM      1  N", text)
        self.assertIn("HETATM    1  C1", text)
        self.assertTrue(text.rstrip().endswith("END"))


if __name__ == "__main__":
    unittest.main()
