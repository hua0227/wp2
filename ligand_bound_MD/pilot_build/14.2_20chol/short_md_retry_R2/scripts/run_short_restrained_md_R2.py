from __future__ import annotations

import csv
import math
import sys
import traceback
from collections import OrderedDict
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass
from pathlib import Path
from typing import Any, NamedTuple


PILOT_ROOT = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol")
ROOT = PILOT_ROOT / "short_md_retry_R2"
SCRIPTS_DIR = ROOT / "scripts"
LOGS_DIR = ROOT / "logs"
OUTPUTS_DIR = ROOT / "outputs"
REPORTS_DIR = ROOT / "reports"

INPUT_PSF = PILOT_ROOT / "targeted_clash_cleanup" / "outputs" / "branch_B_remove_POPC94" / "branch_B_cleaned.psf"
INPUT_PDB = PILOT_ROOT / "restraint_strategy_fix" / "outputs" / "R2_pbcaware_minimized.pdb"
GRO_BOX = Path(r"C:\TRKB_WP2\TRKB_20CHOL\gromacs\step5_input.gro")
TOPPAR_DIR = Path(r"C:\TRKB_WP2\TRKB_20CHOL\toppar")
LIGAND_STR = Path(r"C:\TRKB_WP2\ligand_bound_MD\cgenff_parameterization\returned_str\L002_14.2.str")

LOG_PATH = LOGS_DIR / "short_md_retry_R2.log"
FINAL_PDB = OUTPUTS_DIR / "TRKB_20chol_L002_14.2_shortMD_R2_final.pdb"
DCD_PATH = OUTPUTS_DIR / "TRKB_20chol_L002_14.2_shortMD_R2.dcd"
MONITOR_CSV = REPORTS_DIR / "short_md_retry_R2_monitor.csv"
SUMMARY_MD = REPORTS_DIR / "short_md_retry_R2_summary.md"

PARAMETER_FILES = [
    TOPPAR_DIR / "top_all36_prot.rtf",
    TOPPAR_DIR / "top_all36_lipid.rtf",
    TOPPAR_DIR / "top_all36_cgenff.rtf",
    TOPPAR_DIR / "par_all36m_prot.prm",
    TOPPAR_DIR / "par_all36_lipid.prm",
    TOPPAR_DIR / "par_all36_cgenff.prm",
    TOPPAR_DIR / "toppar_water_ions.str",
    TOPPAR_DIR / "toppar_all36_lipid_cholesterol.str",
    LIGAND_STR,
]

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

PROTEIN_BACKBONE_K = 500.0
LIGAND_HEAVY_K = 25.0
LIPID_CHOL_HEAVY_K = 5.0

TIMESTEP_PS = 0.0005
TOTAL_STEPS = 20000  # 10 ps at 0.5 fs
REPORT_INTERVAL_STEPS = 100
SAFETY_INTERVAL_STEPS = 10
DCD_INTERVAL_STEPS = 1000
RETENTION_CUTOFF_A = 6.0
KEY_DISTANCE_WARNING_A = 10.0
TEMPERATURE_ABORT_K = 500.0
CLASH_SAMPLE_INTERVAL_STEPS = 100


class Stage(NamedTuple):
    name: str
    target_temperature_k: float
    steps: int


STAGES = [
    Stage("50K_2ps", 50.0, 4000),
    Stage("100K_2ps", 100.0, 4000),
    Stage("200K_2ps", 200.0, 4000),
    Stage("310K_4ps", 310.0, 8000),
]


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


class Tee:
    def __init__(self, *streams: Any) -> None:
        self.streams = streams

    def write(self, data: str) -> int:
        for stream in self.streams:
            stream.write(data)
            stream.flush()
        return len(data)

    def flush(self) -> None:
        for stream in self.streams:
            stream.flush()


def ensure_dirs() -> None:
    for path in [ROOT, SCRIPTS_DIR, LOGS_DIR, OUTPUTS_DIR, REPORTS_DIR]:
        path.mkdir(parents=True, exist_ok=True)


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


def is_ligand(atom: PdbAtom) -> bool:
    return atom.segid == "LIG" or atom.resname == "L002"


def is_hydrogen(atom: PdbAtom) -> bool:
    return atom.element.upper() == "H"


def is_heavy(atom: PdbAtom) -> bool:
    return not is_hydrogen(atom)


def is_protein(atom: PdbAtom) -> bool:
    return atom.resname in PROTEIN_RESNAMES or atom.segid.startswith("PRO")


def is_lipid(atom: PdbAtom) -> bool:
    return atom.resname in LIPID_RESNAMES or atom.segid == "MEMB"


def is_cholesterol(atom: PdbAtom) -> bool:
    return atom.resname in CHOLESTEROL_RESNAMES


def is_water_or_ion(atom: PdbAtom) -> bool:
    return atom.resname in WATER_RESNAMES or atom.resname in ION_RESNAMES or atom.segid.startswith("WT") or atom.segid.startswith("ION")


def read_gro_box_lengths_a(path: Path = GRO_BOX) -> tuple[float, float, float]:
    values = [float(value) for value in path.read_text(encoding="utf-8", errors="replace").splitlines()[-1].split()]
    if len(values) < 3:
        raise ValueError(f"Unsupported GRO box vector line with {len(values)} values")
    return values[0] * 10.0, values[1] * 10.0, values[2] * 10.0


def read_gro_box_vectors(path: Path = GRO_BOX) -> Any:
    from openmm import Vec3
    from openmm.unit import nanometer

    values = [float(value) for value in path.read_text(encoding="utf-8", errors="replace").splitlines()[-1].split()]
    if len(values) == 3:
        a, b, c = values
        return (Vec3(a, 0.0, 0.0) * nanometer, Vec3(0.0, b, 0.0) * nanometer, Vec3(0.0, 0.0, c) * nanometer)
    if len(values) == 9:
        v1x, v2y, v3z, v1y, v1z, v2x, v2z, v3x, v3y = values
        return (Vec3(v1x, v1y, v1z) * nanometer, Vec3(v2x, v2y, v2z) * nanometer, Vec3(v3x, v3y, v3z) * nanometer)
    raise ValueError(f"Unsupported GRO box vector line with {len(values)} values")


def set_psf_box_from_gro(psf: Any) -> None:
    from openmm.unit import nanometer

    box_vectors = read_gro_box_vectors()
    if hasattr(psf, "setBoxVectors"):
        psf.setBoxVectors(*box_vectors)
    else:
        psf.setBox(
            box_vectors[0].value_in_unit(nanometer)[0] * nanometer,
            box_vectors[1].value_in_unit(nanometer)[1] * nanometer,
            box_vectors[2].value_in_unit(nanometer)[2] * nanometer,
        )


def choose_opencl_platform(platform_cls: Any) -> Any:
    names = [platform_cls.getPlatform(index).getName() for index in range(platform_cls.getNumPlatforms())]
    if "OpenCL" not in names:
        raise RuntimeError(f"OpenCL platform is required; available platforms: {names}")
    return platform_cls.getPlatformByName("OpenCL")


def build_system() -> tuple[Any, Any, Any]:
    from openmm.app import CharmmParameterSet, CharmmPsfFile, HBonds, PDBFile, PME
    from openmm.unit import nanometer

    psf = CharmmPsfFile(str(INPUT_PSF))
    pdb = PDBFile(str(INPUT_PDB))
    set_psf_box_from_gro(psf)
    params = CharmmParameterSet(*[str(path) for path in PARAMETER_FILES])
    system = psf.createSystem(
        params,
        nonbondedMethod=PME,
        nonbondedCutoff=1.2 * nanometer,
        switchDistance=1.0 * nanometer,
        constraints=HBonds,
        rigidWater=True,
    )
    return psf, pdb, system


def restraint_map(atoms: list[PdbAtom]) -> dict[int, float]:
    mapping: dict[int, float] = {}
    for index, atom in enumerate(atoms):
        k = 0.0
        if is_protein(atom) and atom.name in PROTEIN_BACKBONE_NAMES:
            k = PROTEIN_BACKBONE_K
        elif is_ligand(atom) and is_heavy(atom):
            k = LIGAND_HEAVY_K
        elif (is_lipid(atom) or is_cholesterol(atom)) and is_heavy(atom):
            k = LIPID_CHOL_HEAVY_K
        if k > 0.0:
            mapping[index] = k
    return mapping


def add_r2_pbc_restraints(system: Any, atoms: list[PdbAtom], positions: Any) -> dict[str, int]:
    from openmm import CustomExternalForce
    from openmm.unit import nanometer

    force = CustomExternalForce("0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2")
    force.addPerParticleParameter("k")
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    counts = {
        "protein_backbone_medium": 0,
        "ligand_heavy_weak": 0,
        "lipid_cholesterol_heavy_extremely_weak": 0,
        "water_ion_unrestrained": 0,
        "other_unrestrained": 0,
    }
    for index, atom in enumerate(atoms):
        k = restraint_map([atom]).get(0, 0.0)
        if k == 0.0:
            if is_water_or_ion(atom):
                counts["water_ion_unrestrained"] += 1
            else:
                counts["other_unrestrained"] += 1
            continue
        pos = positions[index].value_in_unit(nanometer)
        force.addParticle(index, [k, pos.x, pos.y, pos.z])
        if is_ligand(atom):
            counts["ligand_heavy_weak"] += 1
        elif is_protein(atom):
            counts["protein_backbone_medium"] += 1
        else:
            counts["lipid_cholesterol_heavy_extremely_weak"] += 1
    system.addForce(force)
    return counts


def positions_to_angstrom(positions: Any) -> list[tuple[float, float, float]]:
    from openmm.unit import angstrom

    return [(pos.x, pos.y, pos.z) for pos in positions.value_in_unit(angstrom)]


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


def residue_order_by_protein_segment(atoms: list[PdbAtom]) -> dict[str, list[tuple[str, str]]]:
    ordered: dict[str, OrderedDict[str, str]] = {}
    for atom in atoms:
        if not is_protein(atom):
            continue
        ordered.setdefault(atom.segid, OrderedDict())
        ordered[atom.segid].setdefault(atom.resid, atom.resname)
    return {segid: list(residues.items()) for segid, residues in ordered.items()}


def key_residue_keys(atoms: list[PdbAtom]) -> list[tuple[str, str, str]]:
    wanted = {13: ("TYR", "TYR13"), 17: ("VAL", "VAL17"), 20: ("SER", "SER20")}
    keys: list[tuple[str, str, str]] = []
    for segid, residues in sorted(residue_order_by_protein_segment(atoms).items()):
        for ordinal, (expected_resname, label) in wanted.items():
            if ordinal - 1 >= len(residues):
                continue
            resid, actual_resname = residues[ordinal - 1]
            if actual_resname == expected_resname:
                keys.append((segid, resid, label))
    return keys


def retention_metrics(
    atoms: list[PdbAtom],
    reference_coords: list[tuple[float, float, float]],
    current_coords: list[tuple[float, float, float]],
    box: tuple[float, float, float],
) -> dict[str, Any]:
    ligand_heavy = [index for index, atom in enumerate(atoms) if is_ligand(atom) and is_heavy(atom)]
    key_rows: list[dict[str, Any]] = []
    for segid, resid, label in key_residue_keys(atoms):
        residue_indices = [
            index for index, atom in enumerate(atoms)
            if atom.segid == segid and atom.resid == resid and is_heavy(atom)
        ]
        min_dist = min(distance_a(current_coords[i], current_coords[j], box) for i in ligand_heavy for j in residue_indices)
        key_rows.append({"label": label, "segid": segid, "resid": resid, "distance_A": min_dist, "within_6A": min_dist <= RETENTION_CUTOFF_A})
    chain_distances: dict[str, float] = {}
    for segid in sorted({atom.segid for atom in atoms if is_protein(atom)}):
        chain_indices = [index for index, atom in enumerate(atoms) if atom.segid == segid and is_protein(atom) and is_heavy(atom)]
        chain_distances[segid] = min(distance_a(current_coords[i], current_coords[j], box) for i in ligand_heavy for j in chain_indices)
    key_min = min(row["distance_A"] for row in key_rows)
    near_chains = sum(1 for value in chain_distances.values() if value <= RETENTION_CUTOFF_A)
    return {
        "key_min_distance_A": key_min,
        "near_chains": near_chains,
        "ligand_rmsd_A": direct_rmsd(ligand_heavy, reference_coords, current_coords),
        "key_rows": key_rows,
        "chain_distances": chain_distances,
        "key_retained": any(row["within_6A"] for row in key_rows),
        "near_two_chains": near_chains >= 2,
    }


def clash_metrics(atoms: list[PdbAtom], current_coords: list[tuple[float, float, float]], box: tuple[float, float, float]) -> dict[str, int]:
    import numpy as np

    ligand_indices = np.array([index for index, atom in enumerate(atoms) if is_ligand(atom)], dtype=np.int64)
    other_indices = np.array([index for index, atom in enumerate(atoms) if not is_ligand(atom)], dtype=np.int64)
    coords = np.asarray(current_coords, dtype=np.float64)
    other_coords = coords[other_indices]
    box_array = np.asarray(box, dtype=np.float64)
    severe = 0
    close = 0
    for ligand_index in ligand_indices:
        delta = coords[ligand_index] - other_coords
        delta -= np.round(delta / box_array) * box_array
        d2 = np.einsum("ij,ij->i", delta, delta)
        close += int(np.count_nonzero(d2 < (1.5 * 1.5)))
        severe += int(np.count_nonzero(d2 < (1.0 * 1.0)))
    return {"severe_lt_1A": severe, "close_lt_1p5A": close}


def estimate_temperature_k(kinetic_energy_kj_mol: float, system: Any) -> float:
    dof = 3 * system.getNumParticles() - system.getNumConstraints()
    for force_index in range(system.getNumForces()):
        if system.getForce(force_index).__class__.__name__ == "CMMotionRemover":
            dof -= 3
            break
    return 2.0 * kinetic_energy_kj_mol / (dof * 0.00831446261815324) if dof > 0 else float("nan")


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_pdb(template_pdb: Path, coords: list[tuple[float, float, float]], output_pdb: Path, remarks: list[str]) -> None:
    atom_iter = iter(coords)
    lines = remarks[:]
    for line in template_pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            x, y, z = next(atom_iter)
            lines.append(line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[54:])
    lines.append("END")
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    output_pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")


def fmt(value: Any, digits: int = 6) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "nan"
    return f"{number:.{digits}g}" if math.isfinite(number) else "nan"


MONITOR_FIELDS = [
    "step",
    "stage",
    "target_temperature_K",
    "actual_temperature_K",
    "potential_energy_kJ_mol",
    "kinetic_energy_kJ_mol",
    "total_energy_kJ_mol",
    "ligand_key_min_distance_A",
    "ligand_near_chains_count",
    "ligand_heavy_rmsd_A",
    "severe_lt_1A",
    "close_lt_1p5A",
    "status",
]


def make_monitor_row(
    step: int,
    stage_name: str,
    target_temperature: float,
    actual_temperature: float,
    potential: float,
    kinetic: float,
    retention: dict[str, Any],
    clashes: dict[str, int],
    status: str,
) -> dict[str, Any]:
    return {
        "step": step,
        "stage": stage_name,
        "target_temperature_K": f"{target_temperature:.6g}",
        "actual_temperature_K": f"{actual_temperature:.6g}" if math.isfinite(actual_temperature) else "nan",
        "potential_energy_kJ_mol": f"{potential:.12g}" if math.isfinite(potential) else "nan",
        "kinetic_energy_kJ_mol": f"{kinetic:.12g}" if math.isfinite(kinetic) else "nan",
        "total_energy_kJ_mol": f"{potential + kinetic:.12g}" if math.isfinite(potential + kinetic) else "nan",
        "ligand_key_min_distance_A": f"{retention['key_min_distance_A']:.6g}" if math.isfinite(retention["key_min_distance_A"]) else "nan",
        "ligand_near_chains_count": retention["near_chains"],
        "ligand_heavy_rmsd_A": f"{retention['ligand_rmsd_A']:.6g}" if math.isfinite(retention["ligand_rmsd_A"]) else "nan",
        "severe_lt_1A": clashes["severe_lt_1A"],
        "close_lt_1p5A": clashes["close_lt_1p5A"],
        "status": status,
    }


def to_int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def to_float(value: Any, default: float = float("nan")) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def write_summary(result: dict[str, Any]) -> None:
    final_row = result.get("final_row") or {}
    lines = [
        "# Short Restrained MD Retry R2 Summary",
        "",
        "## Inputs",
        "",
        f"- PSF: `{INPUT_PSF}`",
        f"- R2 minimized PDB: `{INPUT_PDB}`",
        f"- GRO box: `{GRO_BOX}`",
        f"- CHARMM toppar: `{TOPPAR_DIR}`",
        f"- Ligand stream: `{LIGAND_STR}`",
        "",
        "## R2 Restraint Strategy",
        "",
        "- Protein backbone: medium PBC-aware restraint using `periodicdistance(x,y,z,x0,y0,z0)`.",
        "- Ligand heavy atoms: weak PBC-aware restraint using `periodicdistance(x,y,z,x0,y0,z0)`.",
        "- Lipid/cholesterol heavy atoms: extremely weak PBC-aware restraint using `periodicdistance(x,y,z,x0,y0,z0)`.",
        "- Water/ions: unrestrained.",
        f"- Force constants: protein backbone {PROTEIN_BACKBONE_K:g}, ligand heavy {LIGAND_HEAVY_K:g}, lipid/chol heavy {LIPID_CHOL_HEAVY_K:g} kJ/mol/nm^2.",
        "- The old absolute Cartesian lipid/cholesterol restraint expression was not used.",
        "",
        "## Settings",
        "",
        f"- OpenMM platform used: {result.get('platform', 'not reached')}",
        "- CUDA used: NO",
        "- Nonbonded method: PME",
        "- Nonbonded cutoff: 1.2 nm",
        "- constraints: HBonds",
        "- rigidWater: True",
        f"- Timestep: {TIMESTEP_PS * 1000.0:.3f} fs",
        "- Integrator: LangevinMiddleIntegrator",
        "- Schedule: 50 K 2 ps, 100 K 2 ps, 200 K 2 ps, 310 K 4 ps.",
        "",
        "## Outcome",
        "",
        f"- Completed 10 ps: {'YES' if result.get('completed_10ps') else 'NO'}",
        f"- Last completed step: {result.get('last_step', 0)} / {TOTAL_STEPS}",
        f"- Stop reason: {result.get('stop_reason') or 'none'}",
        f"- NaN/nonfinite detected: {'YES' if result.get('has_nan') else 'NO'}",
        f"- Temperature stable below {TEMPERATURE_ABORT_K:.0f} K: {'YES' if result.get('temperature_stable') else 'NO'}",
        f"- Maximum sampled temperature: {fmt(result.get('max_temperature_K'))} K",
        f"- Final potential energy: {final_row.get('potential_energy_kJ_mol', 'nan')} kJ/mol",
        f"- Ligand min distance to TYR13/VAL17/SER20: {final_row.get('ligand_key_min_distance_A', 'nan')} A",
        f"- Ligand near two chains: {'YES' if to_int(final_row.get('ligand_near_chains_count'), 0) >= 2 else 'NO'}",
        f"- Ligand heavy atom RMSD vs R2 input pose: {final_row.get('ligand_heavy_rmsd_A', 'nan')} A",
        f"- Final severe/close contacts: {final_row.get('severe_lt_1A', 'nan')} / {final_row.get('close_lt_1p5A', 'nan')}",
        f"- Recommend next 50 ps restrained MD pilot: {'YES' if result.get('recommend_50ps') else 'NO'}",
        "",
        "## Output Files",
        "",
        f"- Final PDB: `{FINAL_PDB}`",
        f"- DCD: `{DCD_PATH}`",
        f"- Monitor CSV: `{MONITOR_CSV}`",
        f"- Log: `{LOG_PATH}`",
        "",
        "This stage is a short restrained MD pilot retry only, not production MD.",
    ]
    if not result.get("recommend_50ps"):
        lines.extend(
            [
                "",
                "## Next Step",
                "",
                "- Do not move to the 50 ps pilot until the stop reason or final monitor warnings are reviewed.",
            ]
        )
    else:
        lines.extend(
            [
                "",
                "## Next Step",
                "",
                "- The monitored criteria support a 50 ps restrained MD pilot reattempt using the same R2 PBC-aware restraints; this should remain non-production.",
            ]
        )
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_short_md() -> dict[str, Any]:
    import openmm
    from openmm import LangevinMiddleIntegrator, Platform
    from openmm.app import DCDReporter, Simulation
    from openmm.unit import kelvin, kilojoules_per_mole, picoseconds

    ensure_dirs()
    atoms = read_pdb_atoms(INPUT_PDB)
    psf, pdb, system = build_system()
    if len(atoms) != len(pdb.positions):
        raise ValueError(f"Atom count mismatch: parsed {len(atoms)}, OpenMM positions {len(pdb.positions)}")
    restraint_counts = add_r2_pbc_restraints(system, atoms, pdb.positions)
    box = read_gro_box_lengths_a()
    reference_coords = positions_to_angstrom(pdb.positions)

    print("Short restrained MD retry R2 for 14.2 / L002 / 20chol")
    print("This is not production MD. No CUDA will be used.")
    print("R2 PBC-aware restraint expression: 0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2")
    print(f"Input PSF: {INPUT_PSF}")
    print(f"Input PDB: {INPUT_PDB}")
    print(f"OpenMM version: {openmm.version.version}")
    print("Available OpenMM platforms:")
    for index in range(Platform.getNumPlatforms()):
        print(f"- {Platform.getPlatform(index).getName()}")
    platform = choose_opencl_platform(Platform)
    print(f"Selected OpenMM platform: {platform.getName()}")
    print("R2 restraint counts:")
    for key, value in restraint_counts.items():
        print(f"- {key}: {value}")

    integrator = LangevinMiddleIntegrator(50.0 * kelvin, 10.0 / picoseconds, TIMESTEP_PS * picoseconds)
    integrator.setRandomNumberSeed(142310)
    simulation = Simulation(psf.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    simulation.context.setVelocitiesToTemperature(50.0 * kelvin, 142310)
    simulation.reporters.append(DCDReporter(str(DCD_PATH), DCD_INTERVAL_STEPS))

    rows: list[dict[str, Any]] = []
    last_positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=False)
    last_row: dict[str, Any] | None = None
    max_temperature = float("-inf")
    has_nan = False
    stopped = False
    stop_reason = ""
    key_distance_warning_streak = 0
    global_step = 0

    def sample(stage_name: str, target_temperature: float, force_write: bool) -> tuple[dict[str, Any], bool]:
        nonlocal last_positions, last_row, max_temperature, has_nan, key_distance_warning_streak, stop_reason
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        potential = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        kinetic = state.getKineticEnergy().value_in_unit(kilojoules_per_mole)
        actual_temp = estimate_temperature_k(kinetic, system)
        positions = state.getPositions(asNumpy=False)
        coords = positions_to_angstrom(positions)
        coords_finite = not any(not math.isfinite(value) for xyz in coords for value in xyz)
        values_finite = all(math.isfinite(value) for value in [potential, kinetic, actual_temp]) and coords_finite
        if math.isfinite(actual_temp):
            max_temperature = max(max_temperature, actual_temp)

        if values_finite:
            retention = retention_metrics(atoms, reference_coords, coords, box)
            clashes = clash_metrics(atoms, coords, box)
            last_positions = positions
        else:
            retention = {"key_min_distance_A": float("nan"), "near_chains": 0, "ligand_rmsd_A": float("nan")}
            clashes = {"severe_lt_1A": -1, "close_lt_1p5A": -1}

        statuses: list[str] = []
        abort = False
        if not values_finite:
            statuses.append("nonfinite_energy_temperature_or_coordinates")
            has_nan = True
            abort = True
        if math.isfinite(actual_temp) and actual_temp > TEMPERATURE_ABORT_K:
            statuses.append("temperature_above_500K")
            abort = True
        if retention["near_chains"] == 0:
            statuses.append("ligand_near_chains_0")
            abort = True
        if math.isfinite(retention["key_min_distance_A"]) and retention["key_min_distance_A"] > KEY_DISTANCE_WARNING_A:
            statuses.append("warning_key_distance_gt_10A")
            key_distance_warning_streak += 1
            if key_distance_warning_streak >= 5:
                statuses.append("key_distance_warning_persistent")
                abort = True
        else:
            key_distance_warning_streak = 0
        if not statuses:
            statuses.append("ok")
        row = make_monitor_row(
            global_step,
            stage_name,
            target_temperature,
            actual_temp,
            potential,
            kinetic,
            retention,
            clashes,
            ";".join(statuses),
        )
        last_row = row
        if force_write or abort:
            rows.append(row)
            print(
                "step {step} stage {stage} target {target:.1f}K temp {temp:.3f}K PE {pe:.6g} KE {ke:.6g} "
                "key_min {key:.3f}A chains {chains} RMSD {rmsd:.3f}A severe/close {sev}/{close} status {status}".format(
                    step=global_step,
                    stage=stage_name,
                    target=target_temperature,
                    temp=actual_temp,
                    pe=potential,
                    ke=kinetic,
                    key=retention["key_min_distance_A"],
                    chains=retention["near_chains"],
                    rmsd=retention["ligand_rmsd_A"],
                    sev=clashes["severe_lt_1A"],
                    close=clashes["close_lt_1p5A"],
                    status=row["status"],
                )
            )
        if abort:
            stop_reason = row["status"]
        return row, abort

    sample("initial_50K", 50.0, True)
    for stage in STAGES:
        if stopped:
            break
        integrator.setTemperature(stage.target_temperature_k * kelvin)
        print(f"Starting stage {stage.name}: {stage.target_temperature_k:.1f} K for {stage.steps} steps")
        remaining = stage.steps
        while remaining > 0:
            chunk = min(SAFETY_INTERVAL_STEPS, remaining)
            simulation.step(chunk)
            global_step += chunk
            remaining -= chunk
            force_write = global_step % REPORT_INTERVAL_STEPS == 0 or remaining == 0
            _row, abort = sample(stage.name, stage.target_temperature_k, force_write)
            if abort:
                stopped = True
                print(f"Safety stop triggered at step {global_step}: {stop_reason}")
                break

    completed = (not stopped) and global_step == TOTAL_STEPS
    final_coords = positions_to_angstrom(last_positions)
    write_pdb(
        INPUT_PDB,
        final_coords,
        FINAL_PDB,
        [
            "REMARK short restrained MD retry R2 final safe coordinates",
            "REMARK not production MD",
            f"REMARK completed_10ps {completed}",
            f"REMARK last_step {global_step}",
            f"REMARK stop_reason {stop_reason or 'none'}",
        ],
    )
    write_csv(MONITOR_CSV, rows, MONITOR_FIELDS)
    final_row = last_row or {}
    recommend_50ps = (
        completed
        and not has_nan
        and math.isfinite(max_temperature)
        and max_temperature <= TEMPERATURE_ABORT_K
        and to_int(final_row.get("ligand_near_chains_count"), 0) >= 2
        and to_float(final_row.get("ligand_key_min_distance_A")) <= RETENTION_CUTOFF_A
        and to_int(final_row.get("severe_lt_1A"), 999) == 0
        and to_int(final_row.get("close_lt_1p5A"), 999) == 0
        and to_float(final_row.get("ligand_heavy_rmsd_A")) <= 5.0
    )
    result = {
        "completed_10ps": completed,
        "last_step": global_step,
        "stop_reason": stop_reason,
        "has_nan": has_nan,
        "temperature_stable": math.isfinite(max_temperature) and max_temperature <= TEMPERATURE_ABORT_K,
        "max_temperature_K": max_temperature if math.isfinite(max_temperature) else float("nan"),
        "platform": platform.getName(),
        "final_row": final_row,
        "recommend_50ps": recommend_50ps,
    }
    write_summary(result)
    print("")
    print("Terminal summary")
    print(f"completed_10ps: {completed}")
    print(f"last_step: {global_step}/{TOTAL_STEPS}")
    print(f"last_temperature_K: {final_row.get('actual_temperature_K', 'nan')}")
    print(f"last_potential_energy_kJ_mol: {final_row.get('potential_energy_kJ_mol', 'nan')}")
    print(f"nan_detected: {has_nan}")
    print(f"ligand_final_min_distance_to_key_residues_A: {final_row.get('ligand_key_min_distance_A', 'nan')}")
    print(f"ligand_final_near_chains_count: {final_row.get('ligand_near_chains_count', 'nan')}")
    print(f"ligand_final_RMSD_A: {final_row.get('ligand_heavy_rmsd_A', 'nan')}")
    print(f"recommend_50ps_restrained_MD_pilot: {recommend_50ps}")
    print(f"final PDB: {FINAL_PDB}")
    print(f"DCD: {DCD_PATH}")
    print(f"monitor CSV: {MONITOR_CSV}")
    print(f"summary report: {SUMMARY_MD}")
    del simulation
    del integrator
    return result


def failure_result(error: str) -> dict[str, Any]:
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    result = {
        "completed_10ps": False,
        "last_step": 0,
        "stop_reason": error,
        "has_nan": False,
        "temperature_stable": False,
        "max_temperature_K": float("nan"),
        "platform": "not reached",
        "final_row": {},
        "recommend_50ps": False,
    }
    if not MONITOR_CSV.exists():
        write_csv(MONITOR_CSV, [], MONITOR_FIELDS)
    write_summary(result)
    return result


def main() -> int:
    ensure_dirs()
    with LOG_PATH.open("w", encoding="utf-8") as log_file:
        tee_out = Tee(sys.__stdout__, log_file)
        tee_err = Tee(sys.__stderr__, log_file)
        with redirect_stdout(tee_out), redirect_stderr(tee_err):
            try:
                result = run_short_md()
            except Exception as exc:
                error = repr(exc)
                print("FAILED: short restrained MD retry R2 failed:", error)
                traceback.print_exc()
                result = failure_result(error)
    return 0 if result.get("completed_10ps") and not result.get("has_nan") else 1


if __name__ == "__main__":
    raise SystemExit(main())
