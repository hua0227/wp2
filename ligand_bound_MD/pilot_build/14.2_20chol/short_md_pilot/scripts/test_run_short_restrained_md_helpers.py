from __future__ import annotations

import importlib.util
import math
import unittest
from pathlib import Path


SCRIPT_PATH = Path(__file__).with_name("run_short_restrained_md.py")


def load_module():
    spec = importlib.util.spec_from_file_location("run_short_restrained_md", SCRIPT_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {SCRIPT_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ShortRestrainedMdHelperTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.md = load_module()

    def test_heating_schedule_is_exactly_10_ps_at_half_fs(self):
        total_steps = sum(stage.steps for stage in self.md.HEATING_STAGES)
        total_ps = total_steps * self.md.TIMESTEP_PS
        self.assertEqual(total_steps, 20000)
        self.assertTrue(math.isclose(total_ps, 10.0, rel_tol=0.0, abs_tol=1e-12))
        self.assertEqual([stage.temperature_k for stage in self.md.HEATING_STAGES], [50.0, 100.0, 200.0, 310.0])

    def test_platform_selection_requires_opencl_and_rejects_cuda_only(self):
        self.assertEqual(self.md.select_platform_name(["Reference", "CPU", "OpenCL", "CUDA"]), "OpenCL")
        with self.assertRaises(RuntimeError):
            self.md.select_platform_name(["Reference", "CPU", "CUDA"])

    def test_longer_restrained_md_recommendation_requires_stable_retention(self):
        self.assertTrue(
            self.md.recommend_longer_restrained_md(
                completed=True,
                has_nan=False,
                max_temperature_k=312.0,
                final_key_min_distance_a=3.2,
                final_near_chains=2,
                final_ligand_rmsd_a=2.0,
            )
        )

    def test_nonfinite_sample_is_not_safe_for_final_snapshot(self):
        self.assertTrue(self.md.sample_values_are_finite(1.0, 2.0, 3.0, 50.0))
        self.assertFalse(self.md.sample_values_are_finite(float("nan"), 2.0, 3.0, 50.0))
        self.assertFalse(self.md.sample_values_are_finite(1.0, 2.0, float("inf"), 50.0))
        self.assertFalse(
            self.md.recommend_longer_restrained_md(
                completed=True,
                has_nan=False,
                max_temperature_k=312.0,
                final_key_min_distance_a=3.2,
                final_near_chains=0,
                final_ligand_rmsd_a=2.0,
            )
        )


if __name__ == "__main__":
    unittest.main()
