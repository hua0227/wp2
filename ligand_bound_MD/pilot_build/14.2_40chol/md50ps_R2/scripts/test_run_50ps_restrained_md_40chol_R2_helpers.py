from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).with_name("run_50ps_restrained_md_40chol_R2.py")


def load_module():
    spec = importlib.util.spec_from_file_location("run_50ps_restrained_md_40chol_R2", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_50ps_step_budget_and_reporting_are_bounded():
    module = load_module()
    assert module.TIMESTEP_PS == 0.0005
    assert module.TOTAL_STEPS == 100000
    assert module.TOTAL_STEPS * module.TIMESTEP_PS == 50.0
    assert module.REPORT_INTERVAL_STEPS in {1000, 2000}
    assert module.DCD_INTERVAL_STEPS in {2000, 4000}


def test_monitor_fields_track_ser20_individually_and_two_rmsds():
    module = load_module()
    assert "ligand_distance_SER20_A" in module.MONITOR_FIELDS
    assert "ligand_min_distance_to_key_residues_A" in module.MONITOR_FIELDS
    assert "ligand_heavy_rmsd_vs_50ps_input_A" in module.MONITOR_FIELDS
    assert "ligand_heavy_rmsd_vs_R2_minimized_A" in module.MONITOR_FIELDS


def test_r2_restraints_exclude_water_and_ions_and_use_periodicdistance():
    module = load_module()
    atoms = [
        module.PdbAtom(0, 1, "CA", "VAL", "PROA", "438", "C", 0.0, 0.0, 0.0),
        module.PdbAtom(1, 2, "CB", "VAL", "PROA", "438", "C", 0.0, 0.0, 0.0),
        module.PdbAtom(2, 3, "C1", "L002", "LIG", "1", "C", 0.0, 0.0, 0.0),
        module.PdbAtom(3, 4, "H1", "L002", "LIG", "1", "H", 0.0, 0.0, 0.0),
        module.PdbAtom(4, 5, "C21", "POPC", "MEMB", "113", "C", 0.0, 0.0, 0.0),
        module.PdbAtom(5, 6, "C3", "CHL1", "MEMB", "9", "C", 0.0, 0.0, 0.0),
        module.PdbAtom(6, 7, "OH2", "TIP3", "TIP3", "1", "O", 0.0, 0.0, 0.0),
        module.PdbAtom(7, 8, "SOD", "SOD", "ION", "1", "Na", 0.0, 0.0, 0.0),
    ]
    assert module.restraint_map(atoms) == {0: 500.0, 2: 25.0, 4: 5.0, 5: 5.0}
    assert "periodicdistance" in module.PBC_RESTRAINT_EXPRESSION
    assert "x-x0" not in module.PBC_RESTRAINT_EXPRESSION
