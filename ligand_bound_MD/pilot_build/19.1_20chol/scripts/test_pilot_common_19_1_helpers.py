from pathlib import Path
import importlib.util
import sys


SCRIPT_PATH = Path(__file__).with_name("pilot_common_19_1.py")


def load_module():
    spec = importlib.util.spec_from_file_location("pilot_common_19_1", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_10ps_stage_schedule_is_bounded():
    mod = load_module()
    assert mod.TOTAL_STEPS_10PS == 20000
    assert sum(stage.steps for stage in mod.STAGES_10PS) == mod.TOTAL_STEPS_10PS
    assert mod.PROBE_TIMESTEP_PS * mod.PROBE_TOTAL_STEPS <= 0.2


def test_cleanup_targets_never_delete_protein_and_require_dominant_lipid_source():
    mod = load_module()
    contacts = [
        mod.ClashContact("severe_clash", 0.82, "lipid", "MEMB", "91", "POPC"),
        mod.ClashContact("severe_clash", 0.86, "lipid", "MEMB", "91", "POPC"),
        mod.ClashContact("severe_clash", 0.91, "lipid", "MEMB", "91", "POPC"),
        mod.ClashContact("severe_clash", 0.95, "protein", "PROA", "434", "TYR"),
    ]
    targets, reason, _counts = mod.choose_cleanup_targets(contacts)
    assert targets == [("lipid", "MEMB", "91", "POPC")]
    assert "dominant" in reason
    assert all(target[0] != "protein" for target in targets)


def test_restraint_assignment_excludes_water_and_ions():
    mod = load_module()
    atoms = [
        mod.MockAtom("PROA", "434", "TYR", "CA", "C"),
        mod.MockAtom("LIG", "1", "L003", "C1", "C"),
        mod.MockAtom("MEMB", "91", "POPC", "C2", "C"),
        mod.MockAtom("WT1", "10", "TIP3", "OH2", "O"),
        mod.MockAtom("IONA", "1", "SOD", "SOD", "Na"),
    ]
    restraints = mod.r2_restraint_assignment(atoms)
    assert restraints[0] == mod.PROTEIN_BACKBONE_K
    assert restraints[1] == mod.LIGAND_HEAVY_K
    assert restraints[2] == mod.LIPID_CHOL_HEAVY_K
    assert 3 not in restraints
    assert 4 not in restraints
