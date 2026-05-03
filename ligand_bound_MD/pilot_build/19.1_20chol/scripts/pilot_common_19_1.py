from __future__ import annotations

import csv
import json
import math
import shutil
import subprocess
import traceback
from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, NamedTuple

import openmm
from openmm import LangevinMiddleIntegrator, LocalEnergyMinimizer, Platform, VerletIntegrator, XmlSerializer
from openmm.app import CharmmParameterSet, CharmmPsfFile, DCDReporter, HBonds, PDBFile, PME, Simulation
from openmm.unit import angstrom, kelvin, kilojoules_per_mole, nanometer, picoseconds


VMD_EXE = Path(r"D:\software\VMD2\vmd.exe")
OPENMM_PYTHON = Path(r"C:\Users\14566\miniconda3\envs\trkb_openmm\python.exe")
LIGAND_ID = "19.1"
RESNAME = "L003"
SEGID = "LIG"
CHAIN_ID = "Z"
LIGAND_RESID = "1"
KEY_RESIDUES = {"434": "TYR13", "438": "VAL17", "441": "SER20"}

SEVERE_CUTOFF_A = 1.0
CLOSE_CUTOFF_A = 1.5
RETENTION_CUTOFF_A = 6.0
KEY_DISTANCE_GATE_A = 6.0
PROBE_ABORT_TEMPERATURE_K = 200.0
TENPS_ABORT_TEMPERATURE_K = 500.0

PROTEIN_BACKBONE_K = 500.0
LIGAND_HEAVY_K = 25.0
LIPID_CHOL_HEAVY_K = 5.0
PBC_RESTRAINT_EXPRESSION = "0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2"
MAX_MINIMIZATION_ITERATIONS = 20000

PROBE_TIMESTEP_PS = 0.00005
PROBE_TOTAL_STEPS = 4000
TOTAL_STEPS_10PS = 20000
TIMESTEP_10PS = 0.0005

PROTEIN_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYM", "CYX", "GLN", "GLU", "GLH", "GLY",
    "HIS", "HSD", "HSE", "HSP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
}
LIPID_RESNAMES = {"POP", "POPC", "POPE", "POPS", "POPG", "DOPC", "DPPC", "DMPC"}
CHOLESTEROL_RESNAMES = {"CHL", "CHL1", "CHOL", "CLR"}
WATER_RESNAMES = {"TIP3", "HOH", "WAT", "SOL"}
ION_RESNAMES = {"SOD", "CLA", "POT", "NA", "CL", "K", "CAL", "MG", "ZN"}
PROTEIN_BACKBONE_NAMES = {"N", "CA", "C", "O", "OT1", "OT2"}


class Stage(NamedTuple):
    name: str
    target_temperature_k: float
    steps: int


STAGES_10PS = [
    Stage("50K_2ps", 50.0, 4000),
    Stage("100K_2ps", 100.0, 4000),
    Stage("200K_2ps", 200.0, 4000),
    Stage("310K_4ps", 310.0, 8000),
]


@dataclass(frozen=True)
class MockAtom:
    segid: str
    resid: str
    resname: str
    name: str
    element: str


@dataclass(frozen=True)
class PdbAtom:
    index: int
    serial: int
    name: str
    resname: str
    segid: str
    resid: str
    element: str
    x: float
    y: float
    z: float


@dataclass(frozen=True)
class ClashContact:
    contact_type: str
    distance: float
    other_category: str
    other_segid: str
    other_resid: str
    other_resname: str
    ligand_atom_name: str = ""
    ligand_serial: int = 0
    other_atom_name: str = ""
    other_serial: int = 0


@dataclass(frozen=True)
class PilotConfig:
    system_name: str
    pilot_root: Path
    apo_psf: Path
    apo_template_pdb: Path
    apo_shortmd_pdb: Path
    gro_box: Path
    toppar_dir: Path
    ligand_str: Path
    input_ligand_pdb: Path
    seed: int


def make_config(system_name: str) -> PilotConfig:
    base = Path(r"C:\TRKB_WP2")
    trkb_root = base / ("TRKB_20CHOL" if system_name == "20chol" else "TRKB_40CHOL")
    pilot_root = base / "ligand_bound_MD" / "pilot_build" / f"{LIGAND_ID}_{system_name}"
    return PilotConfig(
        system_name=system_name,
        pilot_root=pilot_root,
        apo_psf=trkb_root / "step5_assembly.psf",
        apo_template_pdb=trkb_root / "step5_assembly.pdb",
        apo_shortmd_pdb=trkb_root / "openmm_short_md_output" / "short_md_final.pdb",
        gro_box=trkb_root / "gromacs" / "step5_input.gro",
        toppar_dir=trkb_root / "toppar",
        ligand_str=base / "ligand_bound_MD" / "cgenff_parameterization" / "returned_str" / f"{RESNAME}_{LIGAND_ID}.str",
        input_ligand_pdb=base / "ligand_bound_MD" / "preflight" / "ligand_pose_allatom_named" / LIGAND_ID / system_name / f"{RESNAME}_{LIGAND_ID}_{system_name}_allatom_named.pdb",
        seed=19120 if system_name == "20chol" else 19140,
    )


def scripts_dir(config: PilotConfig) -> Path:
    return config.pilot_root / "scripts"


def logs_dir(config: PilotConfig) -> Path:
    return config.pilot_root / "logs"


def outputs_dir(config: PilotConfig) -> Path:
    return config.pilot_root / "outputs"


def reports_dir(config: PilotConfig) -> Path:
    return config.pilot_root / "reports"


def parameter_files(config: PilotConfig) -> list[Path]:
    return [
        config.toppar_dir / "top_all36_prot.rtf",
        config.toppar_dir / "top_all36_lipid.rtf",
        config.toppar_dir / "top_all36_cgenff.rtf",
        config.toppar_dir / "par_all36m_prot.prm",
        config.toppar_dir / "par_all36_lipid.prm",
        config.toppar_dir / "par_all36_cgenff.prm",
        config.toppar_dir / "toppar_water_ions.str",
        config.toppar_dir / "toppar_all36_lipid_cholesterol.str",
        config.ligand_str,
    ]


def pathset(config: PilotConfig) -> dict[str, Path]:
    out = outputs_dir(config)
    rep = reports_dir(config)
    scr = scripts_dir(config)
    return {
        "apo_psfgen_ready_pdb": out / f"TRKB_{config.system_name}_psfgen_ready_apo.pdb",
        "ligand_psfgen_ready_pdb": out / f"{RESNAME}_{LIGAND_ID}_{config.system_name}_psfgen_ready.pdb",
        "bound_psf": out / f"TRKB_{config.system_name}_{RESNAME}_{LIGAND_ID}_bound.psf",
        "bound_pdb": out / f"TRKB_{config.system_name}_{RESNAME}_{LIGAND_ID}_bound.pdb",
        "cleaned_psf": out / f"TRKB_{config.system_name}_{RESNAME}_{LIGAND_ID}_bound_cleaned.psf",
        "cleaned_pdb": out / f"TRKB_{config.system_name}_{RESNAME}_{LIGAND_ID}_bound_cleaned.pdb",
        "minimized_pdb": out / f"TRKB_{config.system_name}_{RESNAME}_{LIGAND_ID}_R2_minimized.pdb",
        "probe_final_pdb": out / f"TRKB_{config.system_name}_{RESNAME}_{LIGAND_ID}_very_short_probe_final.pdb",
        "tenps_final_pdb": out / f"TRKB_{config.system_name}_{RESNAME}_{LIGAND_ID}_10ps_R2_final.pdb",
        "tenps_dcd": out / f"TRKB_{config.system_name}_{RESNAME}_{LIGAND_ID}_10ps_R2.dcd",
        "tenps_state_xml": out / f"TRKB_{config.system_name}_{RESNAME}_{LIGAND_ID}_10ps_R2_state.xml",
        "tenps_chk": out / f"TRKB_{config.system_name}_{RESNAME}_{LIGAND_ID}_10ps_R2.chk",
        "build_tcl": scr / f"build_{LIGAND_ID}_{config.system_name}.tcl",
        "cleanup_tcl": scr / f"cleanup_{LIGAND_ID}_{config.system_name}.tcl",
        "preflight_json": rep / "preflight_check.json",
        "preflight_md": rep / "preflight_check.md",
        "readtest_json": rep / "openmm_read_test.json",
        "readtest_md": rep / "openmm_read_test.md",
        "clash_csv": rep / "clash_diagnosis.csv",
        "clash_md": rep / "clash_diagnosis.md",
        "minimization_csv": rep / "R2_minimization_summary.csv",
        "minimization_md": rep / "R2_minimization_summary.md",
        "probe_csv": rep / "very_short_probe.csv",
        "probe_md": rep / "very_short_probe.md",
        "tenps_csv": rep / "pilot_10ps_monitor.csv",
        "tenps_md": rep / "pilot_10ps_summary.md",
        "gate_json": rep / "pilot_gate_summary.json",
        "gate_md": rep / "pilot_gate_summary.md",
        "build_log": logs_dir(config) / "psfgen_build.log",
        "readtest_log": logs_dir(config) / "openmm_read_test.log",
        "cleanup_log": logs_dir(config) / "cleanup.log",
        "minimization_log": logs_dir(config) / "R2_pbcaware_minimization.log",
        "probe_log": logs_dir(config) / "very_short_probe.log",
        "tenps_log": logs_dir(config) / "pilot_10ps.log",
    }


def ensure_dirs(config: PilotConfig) -> None:
    for path in [config.pilot_root, scripts_dir(config), logs_dir(config), outputs_dir(config), reports_dir(config)]:
        path.mkdir(parents=True, exist_ok=True)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def normalize_element(raw: str, atom_name: str = "") -> str:
    token = (raw or atom_name).strip()
    letters = "".join(ch for ch in token if ch.isalpha())
    if not letters:
        return ""
    upper = letters.upper()
    if upper.startswith("CL"):
        return "Cl"
    if upper.startswith("BR"):
        return "Br"
    if upper.startswith("NA"):
        return "Na"
    return upper[0]


def atom_lines(path: Path) -> list[str]:
    return [line.rstrip("\n") for line in path.read_text(encoding="utf-8", errors="replace").splitlines() if line.startswith(("ATOM", "HETATM"))]


def parse_pdb_xyz(line: str) -> tuple[float, float, float]:
    return float(line[30:38]), float(line[38:46]), float(line[46:54])


def parse_pdb_atom_line(line: str, index: int) -> PdbAtom | None:
    if not line.startswith(("ATOM", "HETATM")):
        return None
    try:
        return PdbAtom(
            index=index,
            serial=int(line[6:11]),
            name=line[12:16].strip(),
            resname=line[17:21].strip(),
            segid=line[72:76].strip() if len(line) >= 76 else "",
            resid=line[22:26].strip(),
            element=normalize_element(line[76:78].strip() if len(line) >= 78 else "", line[12:16].strip()),
            x=float(line[30:38]),
            y=float(line[38:46]),
            z=float(line[46:54]),
        )
    except ValueError:
        return None


def read_pdb_atoms(path: Path) -> list[PdbAtom]:
    atoms: list[PdbAtom] = []
    atom_index = 0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            atom = parse_pdb_atom_line(line, atom_index)
            atom_index += 1
            if atom is not None:
                atoms.append(atom)
    return atoms


def parse_psf_atom_count(path: Path) -> int:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    for line in lines:
        if "!NATOM" in line:
            return int(line.split()[0])
    raise ValueError(f"Could not find !NATOM section in {path}")


def parse_str_atom_names(path: Path, resname: str = RESNAME) -> list[str]:
    names: list[str] = []
    in_residue = False
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if line.startswith("RESI"):
            parts = line.split()
            in_residue = len(parts) > 1 and parts[1] == resname
            continue
        if in_residue and line.startswith("ATOM"):
            names.append(line.split()[1])
            continue
        if in_residue and line.startswith(("BOND", "DOUBLE", "IMPR", "CMAP", "DONOR", "ACCEPTOR")):
            break
    return names


def format_atom_name(name: str, element: str) -> str:
    clean = name[:4]
    if len(clean) == 4:
        return clean
    if len(element) == 1:
        return f" {clean:<3}"[:4]
    return f"{clean:<4}"[:4]


def format_pdb_line(record: str, serial: int, atom_name: str, resname: str, chain: str, resid: str, x: float, y: float, z: float, segid: str, element: str) -> str:
    atom_field = format_atom_name(atom_name, element)
    resid_field = int(resid) if str(resid).lstrip("-").isdigit() else 1
    return (
        f"{record:<6}{serial:5d} {atom_field} {resname:>4s}{chain[:1]:1s}{resid_field:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}      {segid:<4s}{element:>2s}"
    )


def to_float(value: Any, default: float = float("nan")) -> float:
    try:
        return float(value)
    except Exception:
        return default


def is_ligand(atom: Any) -> bool:
    return getattr(atom, "segid", "").upper() == SEGID or getattr(atom, "resname", "").upper() == RESNAME


def is_protein(atom: Any) -> bool:
    return getattr(atom, "resname", "").upper() in PROTEIN_RESNAMES or getattr(atom, "segid", "").upper().startswith("PRO")


def is_lipid(atom: Any) -> bool:
    return getattr(atom, "resname", "").upper() in LIPID_RESNAMES or getattr(atom, "segid", "").upper() == "MEMB"


def is_cholesterol(atom: Any) -> bool:
    return getattr(atom, "resname", "").upper() in CHOLESTEROL_RESNAMES


def is_water(atom: Any) -> bool:
    return getattr(atom, "resname", "").upper() in WATER_RESNAMES or getattr(atom, "segid", "").upper().startswith("WT")


def is_ion(atom: Any) -> bool:
    return getattr(atom, "resname", "").upper() in ION_RESNAMES or getattr(atom, "segid", "").upper().startswith("ION")


def is_heavy(atom: Any) -> bool:
    return normalize_element(getattr(atom, "element", ""), getattr(atom, "name", "")).upper() != "H"


def atom_category(atom: Any) -> str:
    if is_ligand(atom):
        return "ligand"
    if is_protein(atom):
        return "protein"
    if is_cholesterol(atom):
        return "cholesterol"
    if is_lipid(atom):
        return "lipid"
    if is_water(atom):
        return "water"
    if is_ion(atom):
        return "ion"
    return "unknown"


def r2_restraint_assignment(atoms: list[Any]) -> dict[int, float]:
    mapping: dict[int, float] = {}
    for index, atom in enumerate(atoms):
        k_value = 0.0
        if is_protein(atom) and getattr(atom, "name", "") in PROTEIN_BACKBONE_NAMES:
            k_value = PROTEIN_BACKBONE_K
        elif is_ligand(atom) and is_heavy(atom):
            k_value = LIGAND_HEAVY_K
        elif (is_lipid(atom) or is_cholesterol(atom)) and is_heavy(atom):
            k_value = LIPID_CHOL_HEAVY_K
        if k_value > 0.0:
            mapping[index] = k_value
    return mapping


def choose_cleanup_targets(contacts: list[ClashContact]) -> tuple[list[tuple[str, str, str, str]], str, Counter[tuple[str, str, str, str]]]:
    severe = [contact for contact in contacts if contact.contact_type == "severe_clash"]
    counts: Counter[tuple[str, str, str, str]] = Counter(
        (contact.other_category, contact.other_segid, contact.other_resid, contact.other_resname)
        for contact in severe
    )
    if not severe or not counts:
        return [], "no severe clashes; copied bound PSF/PDB as cleaned", counts
    top_key, top_count = counts.most_common(1)[0]
    top_category = top_key[0]
    dominant = top_count >= 3 and top_count / len(severe) >= 0.50
    if dominant and top_category in {"water", "ion"}:
        targets = sorted(key for key in counts if key[0] in {"water", "ion"})
        return targets, "deleted severe-clashing water/ion residues from dominant solvent/ion clash source", counts
    if dominant and top_category in {"lipid", "cholesterol"}:
        return [top_key], f"deleted dominant severe-clashing {top_category} molecule {top_key[3]} {top_key[1]}:{top_key[2]}", counts
    return [], "no dominant deletable severe clash source; copied bound PSF/PDB as cleaned", counts


def read_gro_box_vectors(path: Path) -> Any:
    values = [float(value) for value in path.read_text(encoding="utf-8", errors="replace").splitlines()[-1].split()]
    if len(values) == 3:
        a, b, c = values
        return (openmm.Vec3(a, 0.0, 0.0) * nanometer, openmm.Vec3(0.0, b, 0.0) * nanometer, openmm.Vec3(0.0, 0.0, c) * nanometer)
    if len(values) == 9:
        v1x, v2y, v3z, v1y, v1z, v2x, v2z, v3x, v3y = values
        return (
            openmm.Vec3(v1x, v1y, v1z) * nanometer,
            openmm.Vec3(v2x, v2y, v2z) * nanometer,
            openmm.Vec3(v3x, v3y, v3z) * nanometer,
        )
    raise ValueError(f"Unsupported GRO box vector line with {len(values)} values")


def read_gro_box_lengths_a(path: Path) -> tuple[float, float, float]:
    values = [float(value) for value in path.read_text(encoding="utf-8", errors="replace").splitlines()[-1].split()]
    return values[0] * 10.0, values[1] * 10.0, values[2] * 10.0


def set_psf_box_from_gro(psf: Any, gro_path: Path) -> None:
    box_vectors = read_gro_box_vectors(gro_path)
    if hasattr(psf, "setBoxVectors"):
        psf.setBoxVectors(*box_vectors)
    else:
        psf.setBox(
            box_vectors[0].value_in_unit(nanometer)[0] * nanometer,
            box_vectors[1].value_in_unit(nanometer)[1] * nanometer,
            box_vectors[2].value_in_unit(nanometer)[2] * nanometer,
        )


def choose_opencl_platform() -> Any:
    names = [Platform.getPlatform(index).getName() for index in range(Platform.getNumPlatforms())]
    if "OpenCL" not in names:
        raise RuntimeError(f"OpenCL platform required; available platforms: {names}")
    return Platform.getPlatformByName("OpenCL")


def build_system(config: PilotConfig, psf_path: Path, pdb_path: Path, constraints: Any = None, rigid_water: bool = False) -> tuple[Any, Any, Any]:
    psf = CharmmPsfFile(str(psf_path))
    pdb = PDBFile(str(pdb_path))
    set_psf_box_from_gro(psf, config.gro_box)
    params = CharmmParameterSet(*[str(path) for path in parameter_files(config)])
    kwargs = {
        "nonbondedMethod": PME,
        "nonbondedCutoff": 1.2 * nanometer,
        "switchDistance": 1.0 * nanometer,
        "constraints": constraints,
    }
    if rigid_water:
        kwargs["rigidWater"] = True
    system = psf.createSystem(params, **kwargs)
    return psf, pdb, system


def positions_to_angstrom(positions: Any) -> list[tuple[float, float, float]]:
    return [(float(p.x), float(p.y), float(p.z)) for p in positions.value_in_unit(angstrom)]


def distance_a(a: tuple[float, float, float], b: tuple[float, float, float], box: tuple[float, float, float] | None = None) -> float:
    total = 0.0
    for avalue, bvalue, length in zip(a, b, box or (0.0, 0.0, 0.0)):
        delta = avalue - bvalue
        if box is not None and length > 0:
            delta -= round(delta / length) * length
        total += delta * delta
    return math.sqrt(total)


def direct_rmsd(indices: list[int], reference: list[tuple[float, float, float]], current: list[tuple[float, float, float]]) -> float:
    if not indices:
        return float("nan")
    return math.sqrt(sum(distance_a(reference[index], current[index]) ** 2 for index in indices) / len(indices))


def pbc_restraint_force(restraint_map: dict[int, float], positions: Any) -> Any:
    force = openmm.CustomExternalForce(PBC_RESTRAINT_EXPRESSION)
    force.addPerParticleParameter("k")
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    for index, k_value in restraint_map.items():
        pos = positions[index].value_in_unit(nanometer)
        force.addParticle(index, [k_value, pos.x, pos.y, pos.z])
    return force


def coords_have_nan(coords: list[tuple[float, float, float]]) -> bool:
    return any(not math.isfinite(value) for xyz in coords for value in xyz)


def force_rows(atoms: list[PdbAtom], forces: Any) -> tuple[list[dict[str, Any]], bool]:
    vals = forces.value_in_unit(kilojoules_per_mole / nanometer)
    rows: list[dict[str, Any]] = []
    has_nan = False
    for index, (atom, force) in enumerate(zip(atoms, vals)):
        fx, fy, fz = float(force.x), float(force.y), float(force.z)
        mag = math.sqrt(fx * fx + fy * fy + fz * fz)
        finite = all(math.isfinite(value) for value in [fx, fy, fz, mag])
        has_nan = has_nan or not finite
        rows.append(
            {
                "index": index,
                "atom_name": atom.name,
                "resname": atom.resname,
                "segid": atom.segid,
                "resid": atom.resid,
                "force_magnitude_kJ_mol_nm": mag,
            }
        )
    rows.sort(key=lambda row: row["force_magnitude_kJ_mol_nm"], reverse=True)
    return rows, has_nan


def key_distance_metrics(atoms: list[PdbAtom], coords: list[tuple[float, float, float]], box: tuple[float, float, float]) -> dict[str, float]:
    ligand_heavy = [index for index, atom in enumerate(atoms) if is_ligand(atom) and is_heavy(atom)]
    out: dict[str, float] = {}
    for resid, label in KEY_RESIDUES.items():
        key_indices = [index for index, atom in enumerate(atoms) if atom.segid in {"PROA", "PROB"} and atom.resid == resid and is_heavy(atom)]
        out[label] = min(distance_a(coords[i], coords[j], box) for i in ligand_heavy for j in key_indices)
    out["key_min"] = min(out.values())
    return out


def near_chain_count(atoms: list[PdbAtom], coords: list[tuple[float, float, float]], box: tuple[float, float, float]) -> int:
    ligand_heavy = [index for index, atom in enumerate(atoms) if is_ligand(atom) and is_heavy(atom)]
    count = 0
    for segid in sorted({atom.segid for atom in atoms if is_protein(atom)}):
        indices = [index for index, atom in enumerate(atoms) if atom.segid == segid and is_protein(atom) and is_heavy(atom)]
        min_dist = min(distance_a(coords[i], coords[j], box) for i in ligand_heavy for j in indices)
        if min_dist <= RETENTION_CUTOFF_A:
            count += 1
    return count


def clash_metrics(atoms: list[PdbAtom], coords: list[tuple[float, float, float]], box: tuple[float, float, float]) -> dict[str, int]:
    ligand_heavy = [index for index, atom in enumerate(atoms) if is_ligand(atom) and is_heavy(atom)]
    other_heavy = [index for index, atom in enumerate(atoms) if not is_ligand(atom) and is_heavy(atom)]
    severe = 0
    close = 0
    for i in ligand_heavy:
        for j in other_heavy:
            dist = distance_a(coords[i], coords[j], box)
            if dist < CLOSE_CUTOFF_A:
                close += 1
                if dist < SEVERE_CUTOFF_A:
                    severe += 1
    return {"severe_lt_1A": severe, "close_lt_1p5A": close}


def estimate_temperature_k(kinetic_energy_kj_mol: float, system: Any) -> float:
    dof = 3 * system.getNumParticles() - system.getNumConstraints()
    for index in range(system.getNumForces()):
        if system.getForce(index).__class__.__name__ == "CMMotionRemover":
            dof -= 3
            break
    return 2 * kinetic_energy_kj_mol / (dof * 0.00831446261815324)


def write_pdb(template: Path, coords: list[tuple[float, float, float]], output: Path, remarks: list[str]) -> None:
    coord_iter = iter(coords)
    lines = remarks[:]
    for line in template.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            x, y, z = next(coord_iter)
            lines.append(line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[54:])
    lines.append("END")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_command(command: list[str], log_path: Path, cwd: Path | None = None) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(command, cwd=str(cwd) if cwd else None, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    log_path.write_text(completed.stdout, encoding="utf-8", errors="replace")
    return completed.returncode


def prepare_apo_pdb(config: PilotConfig) -> Path:
    paths = pathset(config)
    template_atoms = atom_lines(config.apo_template_pdb)
    short_atoms = atom_lines(config.apo_shortmd_pdb)
    psf_count = parse_psf_atom_count(config.apo_psf)
    if len(template_atoms) != len(short_atoms) or len(template_atoms) != psf_count:
        raise ValueError(f"Atom count mismatch for apo prep: template={len(template_atoms)} short_md={len(short_atoms)} psf={psf_count}")
    lines = [
        f"REMARK psfgen-ready apo coordinate PDB for {LIGAND_ID} / {RESNAME} / {config.system_name}",
        "REMARK atom identity/segid/resid from step5_assembly.pdb; coordinates from short_md_final.pdb",
    ]
    for template_line, coord_line in zip(template_atoms, short_atoms):
        x, y, z = parse_pdb_xyz(coord_line)
        lines.append(template_line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + template_line[54:])
    lines.append("END")
    paths["apo_psfgen_ready_pdb"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    return paths["apo_psfgen_ready_pdb"]


def prepare_ligand_pdb(config: PilotConfig) -> Path:
    paths = pathset(config)
    pose_lines = atom_lines(config.input_ligand_pdb)
    str_names = parse_str_atom_names(config.ligand_str, RESNAME)
    if len(pose_lines) != len(str_names):
        raise ValueError(f"Atom count mismatch for ligand prep: pose={len(pose_lines)} str={len(str_names)}")
    lines = [
        f"REMARK psfgen-ready ligand coordinate PDB for {LIGAND_ID} / {RESNAME} / {config.system_name}",
        "REMARK atom names copied from ligand .str atom order; coordinates from all-atom docked ligand pose",
    ]
    for serial, (line, str_name) in enumerate(zip(pose_lines, str_names), start=1):
        x, y, z = parse_pdb_xyz(line)
        element = normalize_element(line[76:78].strip() if len(line) >= 78 else "", line[12:16].strip())
        lines.append(format_pdb_line("HETATM", serial, str_name, RESNAME, CHAIN_ID, LIGAND_RESID, x, y, z, SEGID, element))
    lines.append("END")
    paths["ligand_psfgen_ready_pdb"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    return paths["ligand_psfgen_ready_pdb"]


def write_build_tcl(config: PilotConfig) -> Path:
    paths = pathset(config)
    lines = [
        f"# Pilot VMD/psfgen build script for ligand {LIGAND_ID} / {RESNAME} / {config.system_name}.",
        "package require psfgen",
        "resetpsf",
        "",
        f"set toppar_dir [file normalize \"{str(config.toppar_dir).replace(chr(92), '/')}\" ]",
        f"set ligand_str [file normalize \"{str(config.ligand_str).replace(chr(92), '/')}\" ]",
        f"set apo_psf [file normalize \"{str(config.apo_psf).replace(chr(92), '/')}\" ]",
        f"set apo_psfgen_pdb [file normalize \"{str(paths['apo_psfgen_ready_pdb']).replace(chr(92), '/')}\" ]",
        f"set ligand_pdb [file normalize \"{str(paths['ligand_psfgen_ready_pdb']).replace(chr(92), '/')}\" ]",
        f"set out_psf [file normalize \"{str(paths['bound_psf']).replace(chr(92), '/')}\" ]",
        f"set out_pdb [file normalize \"{str(paths['bound_pdb']).replace(chr(92), '/')}\" ]",
        "",
        "topology [file join $toppar_dir \"top_all36_prot.rtf\"]",
        "topology [file join $toppar_dir \"top_all36_lipid.rtf\"]",
        "topology [file join $toppar_dir \"top_all36_cgenff.rtf\"]",
        "topology [file join $toppar_dir \"toppar_water_ions.str\"]",
        "topology [file join $toppar_dir \"toppar_all36_lipid_cholesterol.str\"]",
        "topology $ligand_str",
        "",
        "readpsf $apo_psf pdb $apo_psfgen_pdb",
        "",
        "segment LIG {",
        "    first none",
        "    last none",
        f"    residue 1 {RESNAME}",
        "}",
        "coordpdb $ligand_pdb LIG",
        "",
        "guesscoord",
        "writepsf $out_psf",
        "writepdb $out_pdb",
        "puts \"Wrote ligand-bound PSF: $out_psf\"",
        "puts \"Wrote ligand-bound PDB: $out_pdb\"",
        "quit",
        "",
    ]
    paths["build_tcl"].write_text("\n".join(lines), encoding="utf-8")
    return paths["build_tcl"]


def run_vmd_script(script_path: Path, log_path: Path) -> None:
    if not VMD_EXE.exists():
        raise FileNotFoundError(f"VMD executable not found: {VMD_EXE}")
    code = run_command([str(VMD_EXE), "-dispdev", "text", "-e", str(script_path)], log_path, cwd=script_path.parent)
    if code != 0:
        raise RuntimeError(f"VMD script failed: {script_path}")


def run_preflight(config: PilotConfig) -> dict[str, Any]:
    ensure_dirs(config)
    paths = pathset(config)
    str_names = parse_str_atom_names(config.ligand_str, RESNAME)
    pose_atoms_20 = read_pdb_atoms(config.input_ligand_pdb)
    ligand_resnames = sorted({atom.resname for atom in pose_atoms_20})
    params_ok = False
    param_message = ""
    try:
        CharmmParameterSet(*[str(path) for path in parameter_files(config)])
        params_ok = True
        param_message = "CharmmParameterSet loaded successfully"
    except Exception as exc:
        param_message = f"CharmmParameterSet load failed: {exc!r}"
    atom_names_match = len(str_names) == len(pose_atoms_20)
    payload = {
        "ligand_id": LIGAND_ID,
        "resname": RESNAME,
        "system": config.system_name,
        "ligand_str_exists": config.ligand_str.exists(),
        "ligand_pose_exists": config.input_ligand_pdb.exists(),
        "ligand_pose_atom_count": len(pose_atoms_20),
        "ligand_str_atom_count": len(str_names),
        "ligand_pose_resnames": ligand_resnames,
        "ligand_pose_first10_atom_names": [atom.name for atom in pose_atoms_20[:10]],
        "ligand_str_first10_atom_names": str_names[:10],
        "atom_count_matches": atom_names_match,
        "parameter_load_ok": params_ok,
        "parameter_load_message": param_message,
    }
    lines = [
        f"# Preflight Check: {LIGAND_ID} / {RESNAME} / {config.system_name}",
        "",
        f"- Ligand STR exists: {config.ligand_str.exists()}",
        f"- All-atom ligand pose exists: {config.input_ligand_pdb.exists()}",
        f"- Ligand pose atom count: {len(pose_atoms_20)}",
        f"- Ligand STR atom count: {len(str_names)}",
        f"- Ligand pose residue names: {', '.join(ligand_resnames)}",
        f"- Atom count matches STR topology: {atom_names_match}",
        f"- CharmmParameterSet readable: {params_ok}",
        f"- Parameter load note: {param_message}",
    ]
    paths["preflight_md"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    write_json(paths["preflight_json"], payload)
    return payload


def run_openmm_read_test(config: PilotConfig) -> dict[str, Any]:
    paths = pathset(config)
    missing = [str(path) for path in [paths["bound_psf"], paths["bound_pdb"], config.gro_box, *parameter_files(config)] if not path.exists()]
    result = {"success": False, "initial_energy_kJ_mol": float("nan"), "platform": "", "missing": missing}
    lines = [
        f"OpenMM ligand-bound read test for {LIGAND_ID} / {RESNAME} / {config.system_name}",
        "OpenCL only. No minimization and no MD will be run.",
    ]
    if missing:
        lines.append("Missing required files:")
        lines.extend([f"- {path}" for path in missing])
        paths["readtest_log"].write_text("\n".join(lines) + "\n", encoding="utf-8")
        paths["readtest_md"].write_text("\n".join(lines) + "\n", encoding="utf-8")
        write_json(paths["readtest_json"], result)
        return result
    try:
        psf = CharmmPsfFile(str(paths["bound_psf"]))
        pdb = PDBFile(str(paths["bound_pdb"]))
        set_psf_box_from_gro(psf, config.gro_box)
        params = CharmmParameterSet(*[str(path) for path in parameter_files(config)])
        system = psf.createSystem(params, nonbondedMethod=PME, nonbondedCutoff=1.2 * nanometer, switchDistance=1.0 * nanometer, constraints=None)
        platform = choose_opencl_platform()
        integrator = VerletIntegrator(0.001 * picoseconds)
        context = openmm.Context(system, integrator, platform)
        context.setPositions(pdb.positions)
        state = context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        result.update({"success": math.isfinite(energy), "initial_energy_kJ_mol": energy, "platform": platform.getName()})
        lines.extend([
            "CharmmParameterSet loaded.",
            f"Selected platform: {platform.getName()}",
            f"Initial potential energy (kJ/mol): {energy}",
            f"Initial energy finite: {math.isfinite(energy)}",
        ])
        del context
        del integrator
    except Exception as exc:
        result["error"] = repr(exc)
        lines.append(f"FAILED: {exc!r}")
    paths["readtest_log"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    paths["readtest_md"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    write_json(paths["readtest_json"], result)
    return result


def find_contacts(bound_pdb: Path) -> list[ClashContact]:
    atoms = read_pdb_atoms(bound_pdb)
    ligand_heavy = [atom for atom in atoms if is_ligand(atom) and is_heavy(atom)]
    other_heavy = [atom for atom in atoms if not is_ligand(atom) and is_heavy(atom)]
    contacts: list[ClashContact] = []
    close_sq = CLOSE_CUTOFF_A * CLOSE_CUTOFF_A
    severe_sq = SEVERE_CUTOFF_A * SEVERE_CUTOFF_A
    for ligand_atom in ligand_heavy:
        for other_atom in other_heavy:
            dx = ligand_atom.x - other_atom.x
            dy = ligand_atom.y - other_atom.y
            dz = ligand_atom.z - other_atom.z
            d2 = dx * dx + dy * dy + dz * dz
            if d2 < close_sq:
                contacts.append(
                    ClashContact(
                        contact_type="severe_clash" if d2 < severe_sq else "close_contact",
                        distance=math.sqrt(d2),
                        other_category=atom_category(other_atom),
                        other_segid=other_atom.segid,
                        other_resid=other_atom.resid,
                        other_resname=other_atom.resname,
                        ligand_atom_name=ligand_atom.name,
                        ligand_serial=ligand_atom.serial,
                        other_atom_name=other_atom.name,
                        other_serial=other_atom.serial,
                    )
                )
    contacts.sort(key=lambda contact: contact.distance)
    return contacts


def write_clash_reports(config: PilotConfig, contacts: list[ClashContact], targets: list[tuple[str, str, str, str]], cleanup_reason: str, molecule_counts: Counter[tuple[str, str, str, str]]) -> None:
    paths = pathset(config)
    rows = []
    for contact in contacts:
        rows.append(
            {
                "contact_type": contact.contact_type,
                "distance_A": f"{contact.distance:.3f}",
                "category": contact.other_category,
                "ligand_atom_name": contact.ligand_atom_name,
                "ligand_serial": contact.ligand_serial,
                "other_atom_name": contact.other_atom_name,
                "other_serial": contact.other_serial,
                "other_resname": contact.other_resname,
                "other_segid": contact.other_segid,
                "other_resid": contact.other_resid,
            }
        )
    write_csv(paths["clash_csv"], rows)
    severe = [contact for contact in contacts if contact.contact_type == "severe_clash"]
    close = [contact for contact in contacts if contact.distance < CLOSE_CUTOFF_A]
    lines = [
        f"# Clash Diagnosis: {LIGAND_ID} / {RESNAME} / {config.system_name}",
        "",
        f"- Severe clashes <1.0 A: {len(severe)}",
        f"- Close contacts <1.5 A: {len(close)}",
        f"- Cleanup decision: {cleanup_reason}",
    ]
    for target in targets:
        lines.append(f"- Delete target: {target[0]} {target[3]} {target[1]}:{target[2]}")
    if molecule_counts:
        lines.extend(["", "## Dominant Severe Clash Counts", "", "| Category | Molecule | Severe count |", "|---|---|---:|"])
        for (category, segid, resid, resname), count in molecule_counts.most_common(12):
            lines.append(f"| {category} | {resname} {segid}:{resid} | {count} |")
    paths["clash_md"].write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_selection(targets: list[tuple[str, str, str, str]]) -> str:
    return " or ".join(f'(segid "{segid}" and resid "{resid}" and resname "{resname}")' for _category, segid, resid, resname in targets)


def write_cleanup_tcl(config: PilotConfig, targets: list[tuple[str, str, str, str]]) -> Path:
    paths = pathset(config)
    lines = [
        f"# Cleanup Tcl for {LIGAND_ID} / {RESNAME} / {config.system_name}",
        f"set input_psf \"{str(paths['bound_psf']).replace(chr(92), '/')}\"",
        f"set input_pdb \"{str(paths['bound_pdb']).replace(chr(92), '/')}\"",
        f"set output_psf \"{str(paths['cleaned_psf']).replace(chr(92), '/')}\"",
        f"set output_pdb \"{str(paths['cleaned_pdb']).replace(chr(92), '/')}\"",
        "",
    ]
    if not targets:
        lines.extend(
            [
                'puts "No targeted cleanup deletion selected. Copying bound PSF/PDB to cleaned outputs."',
                "file copy -force $input_psf $output_psf",
                "file copy -force $input_pdb $output_pdb",
                "quit",
                "",
            ]
        )
    else:
        selection = make_selection(targets)
        lines.extend(
            [
                "mol new $input_psf type psf waitfor all",
                "mol addfile $input_pdb type pdb waitfor all",
                f"set keep [atomselect top {{not ({selection})}}]",
                "$keep writepsf $output_psf",
                "$keep writepdb $output_pdb",
                "$keep delete",
                "quit",
                "",
            ]
        )
    paths["cleanup_tcl"].write_text("\n".join(lines), encoding="utf-8")
    return paths["cleanup_tcl"]


def run_cleanup(config: PilotConfig, targets: list[tuple[str, str, str, str]]) -> None:
    paths = pathset(config)
    if not targets:
        shutil.copyfile(paths["bound_psf"], paths["cleaned_psf"])
        shutil.copyfile(paths["bound_pdb"], paths["cleaned_pdb"])
        paths["cleanup_log"].write_text("No cleanup deletion selected; copied bound PSF/PDB to cleaned outputs.\n", encoding="utf-8")
        return
    cleanup_tcl = write_cleanup_tcl(config, targets)
    run_vmd_script(cleanup_tcl, paths["cleanup_log"])


def run_minimization(config: PilotConfig) -> dict[str, Any]:
    paths = pathset(config)
    atoms = read_pdb_atoms(paths["cleaned_pdb"])
    psf, pdb, system = build_system(config, paths["cleaned_psf"], paths["cleaned_pdb"], constraints=HBonds, rigid_water=True)
    restraints = r2_restraint_assignment(atoms)
    restraint_force = pbc_restraint_force(restraints, pdb.positions)
    restraint_force.setForceGroup(system.getNumForces())
    system.addForce(restraint_force)
    platform = choose_opencl_platform()
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = openmm.Context(system, integrator, platform)
    context.setPositions(pdb.positions)
    ref_coords = positions_to_angstrom(pdb.positions)
    initial_state = context.getState(getEnergy=True)
    initial_energy = initial_state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    LocalEnergyMinimizer.minimize(context, tolerance=1.0 * kilojoules_per_mole / nanometer, maxIterations=MAX_MINIMIZATION_ITERATIONS)
    state = context.getState(getEnergy=True, getForces=True, getPositions=True)
    final_energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    coords = positions_to_angstrom(state.getPositions(asNumpy=False))
    force_table, force_nan = force_rows(atoms, state.getForces(asNumpy=False))
    box = read_gro_box_lengths_a(config.gro_box)
    key_dist = key_distance_metrics(atoms, coords, box)
    near_chains = near_chain_count(atoms, coords, box)
    clashes = clash_metrics(atoms, coords, box)
    ligand_heavy = [index for index, atom in enumerate(atoms) if is_ligand(atom) and is_heavy(atom)]
    ligand_rmsd = direct_rmsd(ligand_heavy, ref_coords, coords)
    coord_nan = coords_have_nan(coords)
    energy_finite = math.isfinite(final_energy)
    probe_gate_passed = energy_finite and not force_nan and not coord_nan and clashes["severe_lt_1A"] == 0 and near_chains == 2 and all(value <= KEY_DISTANCE_GATE_A for value in key_dist.values())
    write_pdb(paths["cleaned_pdb"], coords, paths["minimized_pdb"], [f"REMARK {LIGAND_ID} {config.system_name} R2 minimized", f"REMARK final energy {final_energy:.12g} kJ/mol"])
    result = {
        "success": energy_finite and not force_nan and not coord_nan,
        "platform": platform.getName(),
        "initial_energy_kJ_mol": initial_energy,
        "final_energy_kJ_mol": final_energy,
        "max_force_kJ_mol_nm": force_table[0]["force_magnitude_kJ_mol_nm"],
        "top_force_atom": f"{force_table[0]['atom_name']} {force_table[0]['resname']} {force_table[0]['segid']}:{force_table[0]['resid']}",
        "coordinate_nan": coord_nan,
        "force_nan": force_nan,
        "ligand_distance_TYR13_A": key_dist["TYR13"],
        "ligand_distance_VAL17_A": key_dist["VAL17"],
        "ligand_distance_SER20_A": key_dist["SER20"],
        "ligand_key_min_distance_A": key_dist["key_min"],
        "ligand_near_chains_count": near_chains,
        "ligand_rmsd_A": ligand_rmsd,
        "severe_lt_1A": clashes["severe_lt_1A"],
        "close_lt_1p5A": clashes["close_lt_1p5A"],
        "probe_gate_passed": probe_gate_passed,
    }
    write_csv(paths["minimization_csv"], [{k: (f"{v:.12g}" if isinstance(v, float) else v) for k, v in result.items()}])
    lines = [
        f"# R2 Minimization Summary: {LIGAND_ID} / {RESNAME} / {config.system_name}",
        "",
        f"- Platform: {platform.getName()}",
        f"- Initial restrained energy: {initial_energy:.12g} kJ/mol",
        f"- Final restrained energy: {final_energy:.12g} kJ/mol",
        f"- Severe/close clashes after minimization: {clashes['severe_lt_1A']}/{clashes['close_lt_1p5A']}",
        f"- Ligand distances TYR13/VAL17/SER20: {key_dist['TYR13']:.6g}/{key_dist['VAL17']:.6g}/{key_dist['SER20']:.6g} A",
        f"- Ligand near chains count: {near_chains}",
        f"- Probe gate passed: {probe_gate_passed}",
        "- R2 periodicdistance restraints used: protein backbone medium, ligand heavy weak, lipid/chol heavy extremely weak, water/ions unrestrained.",
    ]
    paths["minimization_md"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    log_lines = [
        f"Selected platform: {platform.getName()}",
        f"Initial energy: {initial_energy:.12g}",
        f"Final energy: {final_energy:.12g}",
        f"Probe gate passed: {probe_gate_passed}",
    ]
    paths["minimization_log"].write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    del context
    del integrator
    return result


def run_probe(config: PilotConfig) -> dict[str, Any]:
    paths = pathset(config)
    atoms = read_pdb_atoms(paths["minimized_pdb"])
    psf, pdb, system = build_system(config, paths["cleaned_psf"], paths["minimized_pdb"], constraints=HBonds, rigid_water=True)
    system.addForce(pbc_restraint_force(r2_restraint_assignment(atoms), pdb.positions))
    platform = choose_opencl_platform()
    integrator = LangevinMiddleIntegrator(50 * kelvin, 10.0 / picoseconds, PROBE_TIMESTEP_PS * picoseconds)
    integrator.setRandomNumberSeed(config.seed)
    simulation = Simulation(psf.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    simulation.context.setVelocitiesToTemperature(50 * kelvin, config.seed)
    ref_coords = positions_to_angstrom(pdb.positions)
    rows = []
    last_positions = pdb.positions
    stopped = False
    stop_reason = ""
    step = 0
    while step <= PROBE_TOTAL_STEPS:
        state = simulation.context.getState(getEnergy=True, getForces=True, getPositions=True)
        pe = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        ke = state.getKineticEnergy().value_in_unit(kilojoules_per_mole)
        temp = estimate_temperature_k(ke, system)
        coords = positions_to_angstrom(state.getPositions(asNumpy=False))
        box = read_gro_box_lengths_a(config.gro_box)
        key_dist = key_distance_metrics(atoms, coords, box)
        near_chains = near_chain_count(atoms, coords, box)
        clashes = clash_metrics(atoms, coords, box)
        ligand_heavy = [index for index, atom in enumerate(atoms) if is_ligand(atom) and is_heavy(atom)]
        rmsd = direct_rmsd(ligand_heavy, ref_coords, coords)
        statuses = []
        if not all(math.isfinite(value) for value in [pe, ke, temp]) or coords_have_nan(coords):
            statuses.append("nonfinite_state")
        if math.isfinite(temp) and temp > PROBE_ABORT_TEMPERATURE_K:
            statuses.append("temperature_above_200K")
        if near_chains == 0:
            statuses.append("ligand_near_chains_0")
        if not statuses:
            statuses.append("ok")
        row = {
            "step": step,
            "time_ps": f"{step * PROBE_TIMESTEP_PS:.6g}",
            "temperature_K": f"{temp:.6g}",
            "potential_energy_kJ_mol": f"{pe:.12g}",
            "kinetic_energy_kJ_mol": f"{ke:.12g}",
            "ligand_key_min_distance_A": f"{key_dist['key_min']:.6g}",
            "ligand_near_chains_count": near_chains,
            "ligand_rmsd_A": f"{rmsd:.6g}",
            "severe_lt_1A": clashes["severe_lt_1A"],
            "close_lt_1p5A": clashes["close_lt_1p5A"],
            "status": ";".join(statuses),
        }
        if step == 0 or step % 100 == 0 or step == PROBE_TOTAL_STEPS or statuses[0] != "ok":
            rows.append(row)
        last_positions = state.getPositions(asNumpy=False)
        if statuses[0] != "ok":
            stopped = True
            stop_reason = row["status"]
            break
        if step == PROBE_TOTAL_STEPS:
            break
        chunk = min(10, PROBE_TOTAL_STEPS - step)
        simulation.step(chunk)
        step += chunk
    write_csv(paths["probe_csv"], rows)
    write_pdb(paths["minimized_pdb"], positions_to_angstrom(last_positions), paths["probe_final_pdb"], [f"REMARK {LIGAND_ID} {config.system_name} very-short probe final", f"REMARK stop reason {stop_reason or 'none'}"])
    result = {
        "passed": (not stopped) and step == PROBE_TOTAL_STEPS,
        "last_step": step,
        "stop_reason": stop_reason or "none",
        "final_temperature_K": to_float(rows[-1]["temperature_K"]),
        "final_potential_energy_kJ_mol": to_float(rows[-1]["potential_energy_kJ_mol"]),
    }
    lines = [
        f"# Very-short Probe Summary: {LIGAND_ID} / {RESNAME} / {config.system_name}",
        "",
        "- This is a non-production, 50 K, <=0.2 ps probe under R2 PBC-aware restraints.",
        f"- Probe passed: {result['passed']}",
        f"- Stop reason: {result['stop_reason']}",
        f"- Final temperature: {result['final_temperature_K']:.6g} K",
        f"- Final potential energy: {result['final_potential_energy_kJ_mol']:.12g} kJ/mol",
    ]
    paths["probe_md"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    paths["probe_log"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    return result


def run_10ps(config: PilotConfig) -> dict[str, Any]:
    paths = pathset(config)
    atoms = read_pdb_atoms(paths["probe_final_pdb"])
    psf, pdb, system = build_system(config, paths["cleaned_psf"], paths["probe_final_pdb"], constraints=HBonds, rigid_water=True)
    system.addForce(pbc_restraint_force(r2_restraint_assignment(atoms), pdb.positions))
    platform = choose_opencl_platform()
    integrator = LangevinMiddleIntegrator(STAGES_10PS[0].target_temperature_k * kelvin, 10.0 / picoseconds, TIMESTEP_10PS * picoseconds)
    integrator.setRandomNumberSeed(config.seed + 7)
    simulation = Simulation(psf.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    simulation.context.setVelocitiesToTemperature(STAGES_10PS[0].target_temperature_k * kelvin, config.seed + 7)
    simulation.reporters.append(DCDReporter(str(paths["tenps_dcd"]), 1000))
    ref_coords = positions_to_angstrom(pdb.positions)
    box = read_gro_box_lengths_a(config.gro_box)
    rows = []
    last_positions = pdb.positions
    global_step = 0
    stopped = False
    stop_reason = "none"
    current_stage = STAGES_10PS[0]
    stage_offset = 0
    for stage in STAGES_10PS:
        current_stage = stage
        integrator.setTemperature(stage.target_temperature_k * kelvin)
        stage_offset = 0
        while stage_offset < stage.steps:
            chunk = min(50, stage.steps - stage_offset)
            simulation.step(chunk)
            global_step += chunk
            stage_offset += chunk
            if global_step % 100 != 0 and global_step != TOTAL_STEPS_10PS:
                continue
            state = simulation.context.getState(getEnergy=True, getPositions=True, getVelocities=True)
            pe = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
            ke = state.getKineticEnergy().value_in_unit(kilojoules_per_mole)
            temp = estimate_temperature_k(ke, system)
            coords = positions_to_angstrom(state.getPositions(asNumpy=False))
            key_dist = key_distance_metrics(atoms, coords, box)
            near_chains = near_chain_count(atoms, coords, box)
            clashes = clash_metrics(atoms, coords, box)
            ligand_heavy = [index for index, atom in enumerate(atoms) if is_ligand(atom) and is_heavy(atom)]
            rmsd = direct_rmsd(ligand_heavy, ref_coords, coords)
            statuses = []
            if not all(math.isfinite(value) for value in [pe, ke, temp]) or coords_have_nan(coords):
                statuses.append("nonfinite_state")
            if math.isfinite(temp) and temp > TENPS_ABORT_TEMPERATURE_K:
                statuses.append("temperature_above_500K")
            if near_chains == 0:
                statuses.append("ligand_near_chains_0")
            if key_dist["VAL17"] > 10.0:
                statuses.append("warning_VAL17_gt_10A")
            if key_dist["SER20"] > 10.0:
                statuses.append("warning_SER20_gt_10A")
            if clashes["severe_lt_1A"] > 0:
                statuses.append("warning_severe_clash_gt_0")
            if not statuses:
                statuses.append("ok")
            rows.append(
                {
                    "step": global_step,
                    "time_ps": f"{global_step * TIMESTEP_10PS:.6g}",
                    "stage": stage.name,
                    "target_temperature_K": f"{stage.target_temperature_k:.6g}",
                    "temperature_K": f"{temp:.6g}",
                    "potential_energy_kJ_mol": f"{pe:.12g}",
                    "ligand_distance_TYR13_A": f"{key_dist['TYR13']:.6g}",
                    "ligand_distance_VAL17_A": f"{key_dist['VAL17']:.6g}",
                    "ligand_distance_SER20_A": f"{key_dist['SER20']:.6g}",
                    "ligand_key_min_distance_A": f"{key_dist['key_min']:.6g}",
                    "ligand_near_chains_count": near_chains,
                    "ligand_heavy_rmsd_A": f"{rmsd:.6g}",
                    "severe_lt_1A": clashes["severe_lt_1A"],
                    "close_lt_1p5A": clashes["close_lt_1p5A"],
                    "status": ";".join(statuses),
                }
            )
            last_positions = state.getPositions(asNumpy=False)
            if any(token in rows[-1]["status"] for token in ["nonfinite_state", "temperature_above_500K", "ligand_near_chains_0"]):
                stopped = True
                stop_reason = rows[-1]["status"]
                break
        if stopped:
            break
    state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
    paths["tenps_state_xml"].write_text(XmlSerializer.serialize(state), encoding="utf-8")
    simulation.saveCheckpoint(str(paths["tenps_chk"]))
    write_csv(paths["tenps_csv"], rows)
    write_pdb(paths["probe_final_pdb"], positions_to_angstrom(last_positions), paths["tenps_final_pdb"], [f"REMARK {LIGAND_ID} {config.system_name} 10 ps R2 final", f"REMARK stop reason {stop_reason}"])
    final = rows[-1]
    completed = (not stopped) and global_step == TOTAL_STEPS_10PS
    recommend_50ps = completed and final["status"] == "ok" and to_float(final["ligand_key_min_distance_A"]) <= RETENTION_CUTOFF_A and int(final["ligand_near_chains_count"]) == 2 and int(final["severe_lt_1A"]) == 0 and int(final["close_lt_1p5A"]) == 0
    result = {
        "completed_10ps": completed,
        "last_step": global_step,
        "stop_reason": stop_reason,
        "final_temperature_K": to_float(final["temperature_K"]),
        "final_potential_energy_kJ_mol": to_float(final["potential_energy_kJ_mol"]),
        "ligand_distance_TYR13_A": to_float(final["ligand_distance_TYR13_A"]),
        "ligand_distance_VAL17_A": to_float(final["ligand_distance_VAL17_A"]),
        "ligand_distance_SER20_A": to_float(final["ligand_distance_SER20_A"]),
        "ligand_key_min_distance_A": to_float(final["ligand_key_min_distance_A"]),
        "ligand_near_chains_count": int(final["ligand_near_chains_count"]),
        "ligand_rmsd_A": to_float(final["ligand_heavy_rmsd_A"]),
        "severe_lt_1A": int(final["severe_lt_1A"]),
        "close_lt_1p5A": int(final["close_lt_1p5A"]),
        "has_nan": any("nonfinite_state" in row["status"] for row in rows),
        "recommend_50ps": recommend_50ps,
        "platform": platform.getName(),
    }
    summary_lines = [
        f"# 10 ps Restrained MD Pilot Summary: {LIGAND_ID} / {RESNAME} / {config.system_name}",
        "",
        "- OpenCL only.",
        "- R2 PBC-aware restraints.",
        "- Staged heating: 50 K 2 ps, 100 K 2 ps, 200 K 2 ps, 310 K 4 ps.",
        "- This is not production MD.",
        f"- Completed 10 ps: {completed}",
        f"- Stop reason: {stop_reason}",
        f"- Final temperature: {result['final_temperature_K']:.6g} K",
        f"- Final potential energy: {result['final_potential_energy_kJ_mol']:.12g} kJ/mol",
        f"- Final TYR13/VAL17/SER20 distances: {result['ligand_distance_TYR13_A']:.6g}/{result['ligand_distance_VAL17_A']:.6g}/{result['ligand_distance_SER20_A']:.6g} A",
        f"- Ligand near chains count: {result['ligand_near_chains_count']}",
        f"- Ligand RMSD: {result['ligand_rmsd_A']:.6g} A",
        f"- Final severe/close contacts: {result['severe_lt_1A']}/{result['close_lt_1p5A']}",
        f"- Recommend 50 ps: {recommend_50ps}",
    ]
    paths["tenps_md"].write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    paths["tenps_log"].write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    return result


def run_pipeline(config: PilotConfig) -> dict[str, Any]:
    ensure_dirs(config)
    paths = pathset(config)
    result: dict[str, Any] = {
        "system": config.system_name,
        "psf_pdb_generated": False,
        "read_test_succeeded": False,
        "initial_energy_kJ_mol": float("nan"),
        "severe_before_cleanup": 0,
        "close_before_cleanup": 0,
        "cleanup_executed": False,
        "deleted_molecules": [],
        "minimization_succeeded": False,
        "very_short_probe_passed": False,
        "completed_10ps": False,
        "recommend_50ps": False,
        "error": "",
    }
    try:
        preflight = run_preflight(config)
        if not (preflight["parameter_load_ok"] and preflight["atom_count_matches"] and preflight["ligand_pose_exists"] and preflight["ligand_str_exists"]):
            raise RuntimeError("Preflight failed")
        prepare_apo_pdb(config)
        prepare_ligand_pdb(config)
        write_build_tcl(config)
        run_vmd_script(paths["build_tcl"], paths["build_log"])
        result["psf_pdb_generated"] = paths["bound_psf"].exists() and paths["bound_pdb"].exists()
        read_test = run_openmm_read_test(config)
        result["read_test_succeeded"] = bool(read_test["success"])
        result["initial_energy_kJ_mol"] = read_test["initial_energy_kJ_mol"]
        contacts = find_contacts(paths["bound_pdb"])
        severe_before = sum(1 for contact in contacts if contact.contact_type == "severe_clash")
        close_before = sum(1 for contact in contacts if contact.distance < CLOSE_CUTOFF_A)
        targets, cleanup_reason, molecule_counts = choose_cleanup_targets(contacts)
        result["severe_before_cleanup"] = severe_before
        result["close_before_cleanup"] = close_before
        result["cleanup_executed"] = bool(targets)
        result["cleanup_reason"] = cleanup_reason
        result["deleted_molecules"] = [{"category": cat, "segid": segid, "resid": resid, "resname": resname} for cat, segid, resid, resname in targets]
        write_clash_reports(config, contacts, targets, cleanup_reason, molecule_counts)
        run_cleanup(config, targets)
        minim = run_minimization(config)
        result["minimization_succeeded"] = bool(minim["success"])
        result["minimization_final_energy_kJ_mol"] = minim["final_energy_kJ_mol"]
        result["minimization_final_severe_lt_1A"] = minim["severe_lt_1A"]
        result["minimization_final_close_lt_1p5A"] = minim["close_lt_1p5A"]
        if not minim["probe_gate_passed"]:
            raise RuntimeError("Minimization did not pass very-short probe gate")
        probe = run_probe(config)
        result["very_short_probe_passed"] = bool(probe["passed"])
        if not probe["passed"]:
            raise RuntimeError("Very-short probe did not pass")
        tenps = run_10ps(config)
        result.update(tenps)
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
        result["traceback"] = traceback.format_exc()
    write_json(paths["gate_json"], result)
    gate_lines = [
        f"# Pilot Gate Summary: {LIGAND_ID} / {RESNAME} / {config.system_name}",
        "",
        f"- PSF/PDB generated: {result['psf_pdb_generated']}",
        f"- OpenMM read test succeeded: {result['read_test_succeeded']}",
        f"- Initial energy: {result['initial_energy_kJ_mol']}",
        f"- Severe/close before cleanup: {result['severe_before_cleanup']}/{result['close_before_cleanup']}",
        f"- Cleanup executed: {result['cleanup_executed']}",
        f"- Deleted molecules: {result['deleted_molecules']}",
        f"- Minimization succeeded: {result['minimization_succeeded']}",
        f"- Very-short probe passed: {result['very_short_probe_passed']}",
        f"- 10 ps completed: {result['completed_10ps']}",
        f"- Final temperature: {result.get('final_temperature_K', float('nan'))}",
        f"- Final potential energy: {result.get('final_potential_energy_kJ_mol', float('nan'))}",
        f"- Final TYR13/VAL17/SER20 distances: {result.get('ligand_distance_TYR13_A', float('nan'))}/{result.get('ligand_distance_VAL17_A', float('nan'))}/{result.get('ligand_distance_SER20_A', float('nan'))}",
        f"- Near chains count: {result.get('ligand_near_chains_count', 'nan')}",
        f"- Ligand RMSD: {result.get('ligand_rmsd_A', float('nan'))}",
        f"- Final severe/close: {result.get('severe_lt_1A', 'nan')}/{result.get('close_lt_1p5A', 'nan')}",
        f"- Recommend 50 ps: {result['recommend_50ps']}",
    ]
    if result["error"]:
        gate_lines.extend(["", "## Error", "", result["error"]])
    paths["gate_md"].write_text("\n".join(gate_lines) + "\n", encoding="utf-8")
    return result
