from __future__ import annotations

import unittest

import classify_bound_clashes as classify


class BoundClashClassificationTests(unittest.TestCase):
    def test_parse_bound_pdb_atom_reads_resname_segid_and_element(self) -> None:
        line = "ATOM  52405  H27 L002Z   1       5.142   7.042 -10.726  1.00  0.00      LIG  H"

        atom = classify.parse_pdb_atom_line(line)

        self.assertIsNotNone(atom)
        assert atom is not None
        self.assertEqual(atom.resname, "L002")
        self.assertEqual(atom.segid, "LIG")
        self.assertEqual(atom.element, "H")

    def test_residue_classification_keeps_ligand_and_water_ion_distinct(self) -> None:
        self.assertEqual(classify.classify_atom(classify.PdbAtom(1, "CA", "TYR", "PROA", "434", "C", 0, 0, 0)), "protein")
        self.assertEqual(classify.classify_atom(classify.PdbAtom(2, "P", "POPC", "MEMB", "98", "P", 0, 0, 0)), "lipid")
        self.assertEqual(classify.classify_atom(classify.PdbAtom(3, "C3", "CHL1", "MEMB", "87", "C", 0, 0, 0)), "cholesterol")
        self.assertEqual(classify.classify_atom(classify.PdbAtom(4, "OH2", "TIP3", "WT1", "1", "O", 0, 0, 0)), "water")
        self.assertEqual(classify.classify_atom(classify.PdbAtom(5, "SOD", "SOD", "ION", "1", "Na", 0, 0, 0)), "ion")

    def test_cleaned_tcl_deletes_only_severe_water_or_ion_residues(self) -> None:
        severe_water = classify.Contact(
            "severe_clash",
            0.7,
            classify.PdbAtom(10, "C1", "L002", "LIG", "1", "C", 0, 0, 0),
            classify.PdbAtom(20, "OH2", "TIP3", "WT1", "55", "O", 0.7, 0, 0),
            "water",
        )
        close_water = classify.Contact(
            "close_contact",
            1.2,
            classify.PdbAtom(10, "C1", "L002", "LIG", "1", "C", 0, 0, 0),
            classify.PdbAtom(21, "OH2", "TIP3", "WT1", "56", "O", 1.2, 0, 0),
            "water",
        )

        deletions = classify.deletable_water_ion_residues([severe_water, close_water])

        self.assertEqual(deletions, [("WT1", "55", "TIP3", "water")])


if __name__ == "__main__":
    unittest.main()
