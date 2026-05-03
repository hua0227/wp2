from __future__ import annotations

import unittest

import run_restrained_minimization as minim


class RestrainedMinimizationHelperTests(unittest.TestCase):
    def test_key_residue_selection_uses_protein_segment_ordinal_positions(self) -> None:
        atoms = [
            minim.PdbAtom(1, "CA", "MET", "PROA", "422", "C", 0, 0, 0),
            minim.PdbAtom(2, "CA", "TYR", "PROA", "434", "C", 0, 0, 0),
            minim.PdbAtom(3, "CA", "VAL", "PROA", "438", "C", 0, 0, 0),
            minim.PdbAtom(4, "CA", "SER", "PROA", "441", "C", 0, 0, 0),
        ]
        residue_order = {"PROA": [("422", "MET")] + [(str(i), "ALA") for i in range(423, 434)] + [("434", "TYR")] + [(str(i), "ALA") for i in range(435, 438)] + [("438", "VAL"), ("439", "ILE"), ("440", "ALA"), ("441", "SER")]}

        keys = minim.key_residue_keys(residue_order)

        self.assertIn(("PROA", "434", "TYR13"), keys)
        self.assertIn(("PROA", "438", "VAL17"), keys)
        self.assertIn(("PROA", "441", "SER20"), keys)

    def test_choose_platform_prefers_opencl(self) -> None:
        class FakePlatform:
            names = ["CPU", "OpenCL"]

            def __init__(self, name: str) -> None:
                self._name = name

            def getName(self) -> str:
                return self._name

            @classmethod
            def getNumPlatforms(cls) -> int:
                return len(cls.names)

            @classmethod
            def getPlatform(cls, index: int):
                return cls(cls.names[index])

            @classmethod
            def getPlatformByName(cls, name: str):
                if name not in cls.names:
                    raise Exception(name)
                return cls(name)

        self.assertEqual(minim.choose_platform(FakePlatform).getName(), "OpenCL")


if __name__ == "__main__":
    unittest.main()
