from __future__ import annotations

import unittest

import check_ligand_apo_clashes as clashes


class ClashPrecheckTests(unittest.TestCase):
    def test_parse_pdb_atom_handles_four_character_resname_and_chain(self) -> None:
        line = "HETATM    1  C1  L002Z   1       0.000   0.000   0.000  1.00  0.00           C"

        atom = clashes.parse_pdb_atom_line(line)

        self.assertIsNotNone(atom)
        assert atom is not None
        self.assertEqual(atom.resname, "L002")
        self.assertEqual(atom.chain, "Z")
        self.assertEqual(atom.element, "C")
        self.assertEqual(atom.resid, "1")

    def test_classify_residue_groups_expected_apo_components(self) -> None:
        self.assertEqual(clashes.classify_residue("MET"), "protein")
        self.assertEqual(clashes.classify_residue("POP"), "lipid")
        self.assertEqual(clashes.classify_residue("CHL"), "cholesterol")
        self.assertEqual(clashes.classify_residue("HOH"), "water")
        self.assertEqual(clashes.classify_residue("SOD"), "ion")
        self.assertEqual(clashes.classify_residue("XXX"), "unknown")

    def test_find_contacts_counts_severe_as_subset_of_close_contacts(self) -> None:
        ligand = [clashes.PdbAtom(1, "C1", "L002", "Z", "1", "C", 0.0, 0.0, 0.0, "")]
        apo = [
            clashes.PdbAtom(1, "CA", "MET", "A", "1", "C", 0.5, 0.0, 0.0, ""),
            clashes.PdbAtom(2, "C1", "POP", "C", "2", "C", 1.2, 0.0, 0.0, ""),
            clashes.PdbAtom(3, "O", "HOH", "D", "3", "O", 2.0, 0.0, 0.0, ""),
            clashes.PdbAtom(4, "H1", "HOH", "D", "3", "H", 0.1, 0.0, 0.0, ""),
        ]

        contacts = clashes.find_contacts(ligand, apo)
        severe = [contact for contact in contacts if contact.contact_type == "severe_clash"]
        close = [contact for contact in contacts if contact.distance < 1.5]

        self.assertEqual(len(severe), 1)
        self.assertEqual(len(close), 2)
        self.assertEqual([contact.apo_category for contact in contacts], ["protein", "lipid"])


if __name__ == "__main__":
    unittest.main()
