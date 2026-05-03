import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(r"C:\TRKB_WP2\ligand_bound_MD\scripts\check_cgenff_penalties.py")


def load_module():
    spec = importlib.util.spec_from_file_location("check_cgenff_penalties", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class CheckCgenffPenaltiesTests(unittest.TestCase):
    def test_parse_penalties_extracts_summary_and_section_maxima(self):
        module = load_module()
        text = "\n".join(
            [
                "RESI L001           0.000 ! param penalty=  54.800 ; charge penalty=  22.945",
                "ATOM C1     CG331  -0.264 !    0.000   6   C",
                "ATOM C2     CG311  -0.088 !   22.945   6   C",
                "read param card flex append",
                "BONDS",
                "CG A   100.0 1.0 ! L001 , penalty= 8.0",
                "ANGLES",
                "CG A B 40.0 110.0 ! L001 , penalty= 12.5",
                "DIHEDRALS",
                "A B C D 0.0 3 180.0 ! L001 , penalty= 54.8",
                "IMPROPERS",
                "A B C D 0.0 0 0.0 ! L001 , penalty= 101.0",
            ]
        )

        parsed = module.parse_str_penalties(text)

        self.assertEqual(parsed["param_penalty"], 54.8)
        self.assertEqual(parsed["charge_penalty"], 22.945)
        self.assertEqual(parsed["bond_penalty"], 8.0)
        self.assertEqual(parsed["angle_penalty"], 12.5)
        self.assertEqual(parsed["dihedral_penalty"], 54.8)
        self.assertEqual(parsed["improper_penalty"], 101.0)
        self.assertEqual(parsed["max_penalty"], 101.0)
        self.assertEqual(parsed["n_penalty_ge_10"], 6)
        self.assertEqual(parsed["n_penalty_ge_50"], 3)
        self.assertEqual(parsed["n_penalty_ge_100"], 1)

    def test_penalty_level_rules(self):
        module = load_module()

        self.assertEqual(module.penalty_level(9.9), "low")
        self.assertEqual(module.penalty_level(10), "moderate")
        self.assertEqual(module.penalty_level(50), "high")
        self.assertEqual(module.penalty_level(100), "very_high")
        self.assertEqual(module.penalty_level(None), "unknown")


if __name__ == "__main__":
    unittest.main()
