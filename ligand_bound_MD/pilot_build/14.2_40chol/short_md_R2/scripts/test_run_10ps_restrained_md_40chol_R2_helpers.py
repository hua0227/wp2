from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).with_name("run_10ps_restrained_md_40chol_R2.py")


def load_module():
    spec = importlib.util.spec_from_file_location("run_10ps_restrained_md_40chol_R2", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_stages_are_exactly_10ps_at_half_fs():
    module = load_module()
    assert module.TIMESTEP_PS == 0.0005
    assert module.TOTAL_STEPS == 20000
    assert [(stage.name, stage.target_temperature_k, stage.steps) for stage in module.STAGES] == [
        ("50K_2ps", 50.0, 4000),
        ("100K_2ps", 100.0, 4000),
        ("200K_2ps", 200.0, 4000),
        ("310K_4ps", 310.0, 8000),
    ]
    assert sum(stage.steps for stage in module.STAGES) == module.TOTAL_STEPS


def test_r2_restraints_exclude_water_and_ions():
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
    restraints = module.restraint_map(atoms)
    assert restraints == {0: 500.0, 2: 25.0, 4: 5.0, 5: 5.0}


def test_restraint_expression_uses_periodicdistance():
    module = load_module()
    assert "periodicdistance" in module.PBC_RESTRAINT_EXPRESSION
    assert "x-x0" not in module.PBC_RESTRAINT_EXPRESSION
