from __future__ import annotations

import unittest

import prepare_psfgen_apo_pdb as prep


class PreparePsfgenApoPdbTests(unittest.TestCase):
    def test_merge_template_identity_with_shortmd_coordinates(self) -> None:
        template = "ATOM      1  CA  MET   422     -16.402  -0.429 -14.329  1.00  0.00      PROA  "
        coords = "ATOM      1  CA  MET A   1     -15.616  -4.236 -14.419  1.00  0.00           C  "

        merged = prep.merge_atom_line(template, coords)

        self.assertEqual(merged[:30], template[:30])
        self.assertEqual(merged[30:54], coords[30:54])
        self.assertEqual(merged[54:], template[54:])


if __name__ == "__main__":
    unittest.main()
