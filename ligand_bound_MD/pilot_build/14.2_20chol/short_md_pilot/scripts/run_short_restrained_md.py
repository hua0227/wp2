from __future__ import annotations

import csv
import math
import sys
import traceback
from collections import OrderedDict
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from typing import Any, NamedTuple


PILOT_ROOT = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol")
CLASH_OUTPUTS_DIR = PILOT_ROOT / "clash_resolution" / "outputs"
SHORT_ROOT = PILOT_ROOT / "short_md_pilot"
SCRIPTS_DIR = SHORT_ROOT / "scripts"
LOGS_DIR = SHORT_ROOT / "logs"
OUTPUTS_DIR = SHORT_ROOT / "outputs"
REPORTS_DIR = SHORT_ROOT / "reports"

CLEANED_PSF = CLASH_OUTPUTS_DIR / "TRKB_20chol_L002_14.2_bound_cleaned.psf"
MINIMIZED_PDB = CLASH_OUTPUTS_DIR / "TRKB_20chol_L002_14.2_minimized.pdb"
GRO_BOX = Path(r"C:\TRKB_WP2\TRKB_20CHOL\gromacs\step5_input.gro")
TOPPAR_DIR = Path(r"C:\TRKB_WP2\TRKB_20CHOL\toppar")
LIGAND_STR = Path(r"C:\TRKB_WP2\ligand_bound_MD\cgenff_parameterization\returned_str\L002_14.2.str")

LOG_PATH = LOGS_DIR / "short_restrained_md.log"
FINAL_PDB = OUTPUTS_DIR / "TRKB_20chol_L002_14.2_shortMD_final.pdb"
DCD_PATH = OUTPUTS_DIR / "TRKB_20chol_L002_14.2_shortMD.dcd"
RETENTION_CSV = REPORTS_DIR / "short_md_retention_summary.csv"
SUMMARY_MD = REPORTS_DIR / "short_md_pilot_summary.md"

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
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "ASH",
    "CYS",
    "CYM",
    "CYX",
    "GLN",
    "GLU",
    "GLH",
    "GLY",
    "HIS",
    "HSD",
    "HSE",
    "HSP",
    "ILE",
    "LEU",
    "LYS",
    "LYN",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}
LIPID_RESNAMES = {"POP", "POPC", "POPE", "POPS", "POPG", "DOPC", "DPPC", "DMPC"}
CHOLESTEROL_RESNAMES = {"CHL", "CHL1", "CHOL", "CLR"}
WATER_RESNAMES = {"TIP3", "HOH", "WAT", "SOL"}
ION_RESNAMES = {"SOD", "CLA", "POT", "NA", "CL", "K", "CAL", "MG", "ZN"}
PROTEIN_BACKBONE_NAMES = {"N", "CA", "C", "O", "OT1", "OT2"}

STRONG_RESTRAINT_K = 1000.0  # kJ/mol/nm^2
LIGAND_RESTRAINT_K = 150.0  # kJ/mol/nm^2
RETENTION_CUTOFF_A = 6.0
KEY_DISTANCE_WARNING_A = 10.0
TEMPERATURE_ABORT_K = 500.0
TIMESTEP_PS = 0.0005
REPORT_INTERVAL_STEPS = 100
SAFETY_CHECK_INTERVAL_STEPS = 10
DCD_INTERVAL_STEPS = 1000


class HeatingStage(NamedTuple):
    name: str
    temperature_k: float
    steps: int


class PdbAtom(NamedTuple):
    serial: int
    name: str
    resname: str
    segid: str
    resid: str
    element: str
    x: float
    y: float
    z: float


HEATING_STAGES = [
    HeatingStage("50K_2ps", 50.0, 4000),
    HeatingStage("100K_2ps", 100.0, 4000),
    HeatingStage("200K_2ps", 200.0, 4000),
    HeatingStage("310K_4ps", 310.0, 8000),
]


RETENTION_FIELDS = [
    "stage",
    "step",
    "temperature_K",
    "potential_energy_kJ_mol",
    "total_energy_kJ_mol",
    "ligand_min_distance_to_key_residues_A",
    "ligand_near_chains",
    "ligand_rmsd_A",
    "status",
]


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data: str) -> int:
        for stream in self.streams:
            stream.write(data)
            stream.flush()
        return len(data)

    def flush(self) -> None:
        for stream in self.streams:
            stream.flush()


def select_platform_name(available_platform_names: list[str]) -> str:
    if "OpenCL" not in available_platform_names:
        raise RuntimeError(f"OpenCL platform is required; available platforms: {available_platform_names}")
    return "OpenCL"


def recommend_longer_restrained_md(
    completed: bool,
    has_nan: bool,
    max_temperature_k: float,
    final_key_min_distance_a: float,
    final_near_chains: int,
    final_ligand_rmsd_a: float,
) -> bool:
    return (
        completed
        and not has_nan
        and math.isfinite(max_temperature_k)
        and max_temperature_k <= TEMPERATURE_ABORT_K
        and math.isfinite(final_key_min_distance_a)
        and final_key_min_distance_a <= RETENTION_CUTOFF_A
        and final_near_chains >= 2
        and math.isfinite(final_ligand_rmsd_a)
        and final_ligand_rmsd_a <= 5.0
    )


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


def parse_pdb_atom_line(line: str) -> PdbAtom | None:
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None
    try:
        serial = int(line[6:11])
        name = line[12:16].strip()
        resname = line[17:21].strip()
        resid = line[22:26].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        return None
    segid = line[72:76].strip() if len(line) >= 76 else ""
    element = normalize_element(line[76:78].strip() if len(line) >= 78 else "", name)
    return PdbAtom(serial, name, resname, segid, resid, element, x, y, z)


def read_pdb_atoms(path: Path) -> list[PdbAtom]:
    atoms: list[PdbAtom] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        atom = parse_pdb_atom_line(line)
        if atom is not None:
            atoms.append(atom)
    return atoms


def is_ligand(atom: PdbAtom) -> bool:
    return atom.segid == "LIG" or atom.resname == "L002"


def is_heavy(atom: PdbAtom) -> bool:
    return atom.element != "H"


def is_protein(atom: PdbAtom) -> bool:
    return atom.resname in PROTEIN_RESNAMES or atom.segid.startswith("PRO")


def is_lipid_or_cholesterol(atom: PdbAtom) -> bool:
    return (
        atom.resname in LIPID_RESNAMES
        or atom.resname in CHOLESTEROL_RESNAMES
        or atom.segid == "MEMB"
    )


def is_water_or_ion(atom: PdbAtom) -> bool:
    return (
        atom.resname in WATER_RESNAMES
        or atom.resname in ION_RESNAMES
        or atom.segid.startswith("WT")
        or atom.segid.startswith("ION")
    )


def residue_order_by_protein_segment(atoms: list[PdbAtom]) -> dict[str, list[tuple[str, str]]]:
    ordered: dict[str, OrderedDict[str, str]] = {}
    for atom in atoms:
        if not is_protein(atom):
            continue
        ordered.setdefault(atom.segid, OrderedDict())
        ordered[atom.segid].setdefault(atom.resid, atom.resname)
    return {segid: list(residues.items()) for segid, residues in ordered.items()}


def key_residue_keys(residue_order: dict[str, list[tuple[str, str]]]) -> list[tuple[str, str, str]]:
    wanted = {13: ("TYR", "TYR13"), 17: ("VAL", "VAL17"), 20: ("SER", "SER20")}
    keys: list[tuple[str, str, str]] = []
    for segid, residues in sorted(residue_order.items()):
        for ordinal, (expected_resname, label) in wanted.items():
            if ordinal - 1 >= len(residues):
                continue
            resid, actual_resname = residues[ordinal - 1]
            if actual_resname == expected_resname:
                keys.append((segid, resid, label))
    return keys


def read_gro_box_vectors(path: Path):
    from openmm import Vec3
    from openmm.unit import nanometer

    values = [float(value) for value in path.read_text(encoding="utf-8", errors="replace").splitlines()[-1].split()]
    if len(values) == 3:
        ax, by, cz = values
        return (
            Vec3(ax, 0.0, 0.0) * nanometer,
            Vec3(0.0, by, 0.0) * nanometer,
            Vec3(0.0, 0.0, cz) * nanometer,
        )
    if len(values) == 9:
        v1x, v2y, v3z, v1y, v1z, v2x, v2z, v3x, v3y = values
        return (
            Vec3(v1x, v1y, v1z) * nanometer,
            Vec3(v2x, v2y, v2z) * nanometer,
            Vec3(v3x, v3y, v3z) * nanometer,
        )
    raise ValueError(f"Unsupported GRO box vector line with {len(values)} values")


def restraint_k_for_atom(atom: PdbAtom) -> float | None:
    if is_ligand(atom) and is_heavy(atom):
        return LIGAND_RESTRAINT_K
    if is_lipid_or_cholesterol(atom) and is_heavy(atom):
        return STRONG_RESTRAINT_K
    if is_protein(atom) and atom.name in PROTEIN_BACKBONE_NAMES:
        return STRONG_RESTRAINT_K
    return None


def add_position_restraints(system: Any, atoms: list[PdbAtom], positions: Any) -> tuple[int, int, int]:
    from openmm import CustomExternalForce
    from openmm.unit import nanometer

    force = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addPerParticleParameter("k")
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    strong_count = 0
    ligand_count = 0
    unrestrained_count = 0
    for index, atom in enumerate(atoms):
        k = restraint_k_for_atom(atom)
        if k is None:
            unrestrained_count += 1
            continue
        pos = positions[index].value_in_unit(nanometer)
        force.addParticle(index, [k, pos.x, pos.y, pos.z])
        if is_ligand(atom):
            ligand_count += 1
        else:
            strong_count += 1
    system.addForce(force)
    return strong_count, ligand_count, unrestrained_count


def distance_a(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def positions_to_angstrom(positions: Any) -> list[tuple[float, float, float]]:
    from openmm.unit import angstrom

    return [(pos.x, pos.y, pos.z) for pos in positions.value_in_unit(angstrom)]


def ligand_heavy_indices(atoms: list[PdbAtom]) -> list[int]:
    return [index for index, atom in enumerate(atoms) if is_ligand(atom) and is_heavy(atom)]


def direct_rmsd_a(
    indices: list[int],
    before: list[tuple[float, float, float]],
    after: list[tuple[float, float, float]],
) -> float:
    if not indices:
        return float("nan")
    total = 0.0
    for index in indices:
        total += distance_a(before[index], after[index]) ** 2
    return math.sqrt(total / len(indices))


def min_distance_for_indices(
    ligand_indices: list[int],
    other_indices: list[int],
    coords_a: list[tuple[float, float, float]],
) -> float:
    if not ligand_indices or not other_indices:
        return float("nan")
    best = float("inf")
    for ligand_index in ligand_indices:
        ligand_coord = coords_a[ligand_index]
        for other_index in other_indices:
            best = min(best, distance_a(ligand_coord, coords_a[other_index]))
    return best


def retention_metrics(
    atoms: list[PdbAtom],
    reference_coords_a: list[tuple[float, float, float]],
    current_coords_a: list[tuple[float, float, float]],
) -> dict[str, Any]:
    ligand_indices = ligand_heavy_indices(atoms)
    residue_order = residue_order_by_protein_segment(atoms)
    key_keys = key_residue_keys(residue_order)
    key_distances: list[dict[str, Any]] = []
    for segid, resid, label in key_keys:
        residue_indices = [
            index
            for index, atom in enumerate(atoms)
            if atom.segid == segid and atom.resid == resid and is_heavy(atom)
        ]
        min_dist = min_distance_for_indices(ligand_indices, residue_indices, current_coords_a)
        key_distances.append({"label": label, "segid": segid, "resid": resid, "distance_a": min_dist})

    protein_chain_distances: dict[str, float] = {}
    for segid in sorted({atom.segid for atom in atoms if is_protein(atom)}):
        chain_indices = [
            index
            for index, atom in enumerate(atoms)
            if atom.segid == segid and is_protein(atom) and is_heavy(atom)
        ]
        protein_chain_distances[segid] = min_distance_for_indices(ligand_indices, chain_indices, current_coords_a)

    finite_key_distances = [row["distance_a"] for row in key_distances if math.isfinite(row["distance_a"])]
    key_min_distance = min(finite_key_distances) if finite_key_distances else float("nan")
    near_chain_count = sum(1 for value in protein_chain_distances.values() if math.isfinite(value) and value <= RETENTION_CUTOFF_A)
    ligand_rmsd = direct_rmsd_a(ligand_indices, reference_coords_a, current_coords_a)
    return {
        "key_min_distance_a": key_min_distance,
        "near_chain_count": near_chain_count,
        "ligand_rmsd_a": ligand_rmsd,
        "key_distances": key_distances,
        "protein_chain_distances": protein_chain_distances,
    }


def estimate_temperature_k(kinetic_energy_kj_mol: float, system: Any) -> float:
    dof = 3 * system.getNumParticles() - system.getNumConstraints()
    for force_index in range(system.getNumForces()):
        if system.getForce(force_index).__class__.__name__ == "CMMotionRemover":
            dof -= 3
            break
    if dof <= 0:
        return float("nan")
    gas_constant_kj_mol_k = 0.00831446261815324
    return 2.0 * kinetic_energy_kj_mol / (dof * gas_constant_kj_mol_k)


def status_for_sample(
    potential_energy: float,
    kinetic_energy: float,
    total_energy: float,
    temperature_k: float,
    key_min_distance_a: float,
    near_chain_count: int,
) -> tuple[str, bool, bool]:
    statuses: list[str] = []
    nonfinite = not sample_values_are_finite(potential_energy, kinetic_energy, total_energy, temperature_k)
    abort = False
    if nonfinite:
        statuses.append("nonfinite_energy_or_temperature")
        abort = True
    if math.isfinite(temperature_k) and temperature_k > TEMPERATURE_ABORT_K:
        statuses.append("temperature_above_500K")
        abort = True
    if math.isfinite(key_min_distance_a) and key_min_distance_a > KEY_DISTANCE_WARNING_A:
        statuses.append("warning_key_distance_gt_10A")
    if near_chain_count == 0:
        statuses.append("warning_ligand_near_chains_0")
    if not statuses:
        statuses.append("ok")
    return ";".join(statuses), abort, nonfinite


def sample_values_are_finite(
    potential_energy: float,
    kinetic_energy: float,
    total_energy: float,
    temperature_k: float,
) -> bool:
    return all(math.isfinite(value) for value in [potential_energy, kinetic_energy, total_energy, temperature_k])


def write_pdb_with_positions(
    template_pdb: Path,
    positions: Any,
    output_pdb: Path,
    extra_remarks: list[str] | None = None,
) -> None:
    from openmm.unit import angstrom

    coords = positions.value_in_unit(angstrom)
    atom_line_iter = iter(coords)
    lines: list[str] = [
        "REMARK short restrained MD pilot output for 14.2 / L002 / 20chol",
        "REMARK not production MD",
    ]
    if extra_remarks:
        lines.extend(extra_remarks)
    for line in template_pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            pos = next(atom_line_iter)
            lines.append(line[:30] + f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}" + line[54:])
    lines.append("END")
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    output_pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_summary(result: dict[str, Any]) -> None:
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    recommendation = result.get("recommend_longer_restrained_md", False)
    final_metrics = result.get("final_metrics") or {}
    key_rows = final_metrics.get("key_distances") or []
    chain_rows = final_metrics.get("protein_chain_distances") or {}
    error = result.get("error")
    lines = [
        "# Short Restrained MD Pilot Summary",
        "",
        "## Purpose",
        "",
        "This pilot tests whether the minimized ligand-bound 14.2 / L002 / 20chol system can tolerate a very short restrained NVT run while retaining the docked ligand near the expected TRKB pocket.",
        "",
        "## Inputs",
        "",
        f"- Cleaned PSF: `{CLEANED_PSF}`",
        f"- Minimized PDB: `{MINIMIZED_PDB}`",
        f"- GRO box: `{GRO_BOX}`",
        f"- CHARMM toppar: `{TOPPAR_DIR}`",
        f"- Ligand stream: `{LIGAND_STR}`",
        "",
        "## Settings",
        "",
        f"- OpenMM platform requested: OpenCL",
        f"- OpenMM platform used: {result.get('platform_name', 'not reached')}",
        f"- Timestep: {TIMESTEP_PS * 1000.0:.3f} fs",
        f"- Nonbonded method: PME",
        f"- Nonbonded cutoff: 1.2 nm",
        f"- Constraints: HBonds",
        f"- rigidWater: True",
        f"- Protein backbone restraint: {STRONG_RESTRAINT_K:.1f} kJ/mol/nm^2",
        f"- Lipid/cholesterol heavy atom restraint: {STRONG_RESTRAINT_K:.1f} kJ/mol/nm^2",
        f"- Ligand heavy atom restraint: {LIGAND_RESTRAINT_K:.1f} kJ/mol/nm^2",
        f"- Water/ion restraints: none",
        "",
        "## Outcome",
        "",
        f"- Completed 10 ps: {'YES' if result.get('completed') else 'NO'}",
        f"- Last completed step: {result.get('last_step', 0)} / {sum(stage.steps for stage in HEATING_STAGES)}",
        f"- Last finite snapshot step: {result.get('last_finite_step', 'not available')}",
        f"- Safety stop reason: {result.get('stop_reason') or 'none'}",
        f"- NaN/nonfinite detected: {'YES' if result.get('has_nan') else 'NO'}",
        f"- Maximum sampled temperature: {format_float(result.get('max_temperature_k'))} K",
        f"- Temperature stable below {TEMPERATURE_ABORT_K:.0f} K: {'YES' if result.get('temperature_stable') else 'NO'}",
        f"- Potential energy remained finite: {'YES' if result.get('potential_energy_finite') else 'NO'}",
        f"- Last potential energy: {format_float(result.get('last_potential_energy_kj_mol'))} kJ/mol",
        f"- Last total energy: {format_float(result.get('last_total_energy_kj_mol'))} kJ/mol",
        f"- Last finite potential energy: {format_float(result.get('last_finite_potential_energy_kj_mol'))} kJ/mol",
        f"- Ligand min distance to TYR13/VAL17/SER20: {format_float(final_metrics.get('key_min_distance_a'))} A",
        f"- Ligand near two chains: {'YES' if final_metrics.get('near_chain_count', 0) >= 2 else 'NO'}",
        f"- Ligand heavy atom RMSD relative to minimized pose: {format_float(final_metrics.get('ligand_rmsd_a'))} A",
        f"- Recommend next longer restrained MD: {'YES' if recommendation else 'NO'}",
        "",
        "## Final Key-Residue Distances",
        "",
        "| Key residue | Segid | Resid | Min distance A | Within 6 A |",
        "|---|---|---|---:|---:|",
    ]
    if key_rows:
        for row in key_rows:
            dist = row["distance_a"]
            lines.append(
                f"| {row['label']} | {row['segid']} | {row['resid']} | {format_float(dist)} | {int(math.isfinite(dist) and dist <= RETENTION_CUTOFF_A)} |"
            )
    else:
        lines.append("| not available |  |  | nan | 0 |")
    lines.extend(["", "## Final Chain Distances", "", "| Chain | Min distance A | Within 6 A |", "|---|---:|---:|"])
    if chain_rows:
        for segid, dist in chain_rows.items():
            lines.append(f"| {segid} | {format_float(dist)} | {int(math.isfinite(dist) and dist <= RETENTION_CUTOFF_A)} |")
    else:
        lines.append("| not available | nan | 0 |")
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- This stage is a short restrained MD pilot only, not production MD.",
            "- No production MD was run.",
            "- No original apo files were modified.",
            "- The ligand `.str` file was not modified.",
        ]
    )
    if error:
        lines.extend(["", "## Failure / Next Step", "", f"- Key error: `{error}`"])
    elif not recommendation:
        lines.extend(
            [
                "",
                "## Next Step",
                "",
                "- Review the warnings in the retention CSV and log before deciding whether to repeat this pilot with adjusted restraints or shorter heating increments.",
            ]
        )
    else:
        lines.extend(
            [
                "",
                "## Next Step",
                "",
                "- The pilot criteria support moving to the next longer restrained MD pilot, still with careful retention and energy monitoring.",
            ]
        )
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def format_float(value: Any) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "nan"
    if not math.isfinite(number):
        return "nan"
    return f"{number:.6g}"


def validate_inputs() -> None:
    missing = [path for path in [CLEANED_PSF, MINIMIZED_PDB, GRO_BOX, LIGAND_STR, *PARAMETER_FILES] if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required input files: " + "; ".join(str(path) for path in missing))


def choose_opencl_platform(platform_cls: Any) -> Any:
    names = [platform_cls.getPlatform(index).getName() for index in range(platform_cls.getNumPlatforms())]
    selected = select_platform_name(names)
    return platform_cls.getPlatformByName(selected)


def set_psf_box_from_gro(psf: Any) -> None:
    from openmm.unit import nanometer

    box_vectors = read_gro_box_vectors(GRO_BOX)
    if hasattr(psf, "setBoxVectors"):
        psf.setBoxVectors(*box_vectors)
    else:
        a = box_vectors[0].value_in_unit(nanometer)[0]
        b = box_vectors[1].value_in_unit(nanometer)[1]
        c = box_vectors[2].value_in_unit(nanometer)[2]
        psf.setBox(a * nanometer, b * nanometer, c * nanometer)


def run_short_md() -> dict[str, Any]:
    import openmm
    from openmm import LangevinMiddleIntegrator, Platform
    from openmm.app import CharmmParameterSet, CharmmPsfFile, DCDReporter, HBonds, PDBFile, PME, Simulation
    from openmm.unit import kelvin, kilojoules_per_mole, nanometer, picoseconds

    validate_inputs()
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)

    print("Short restrained MD pilot for 14.2 / L002 / 20chol")
    print("This is not production MD. CUDA will not be used.")
    print("Available OpenMM platforms:")
    platform_names = [Platform.getPlatform(index).getName() for index in range(Platform.getNumPlatforms())]
    for name in platform_names:
        print(f"- {name}")
    platform = choose_opencl_platform(Platform)
    print(f"Selected OpenMM platform: {platform.getName()}")
    print(f"OpenMM version: {openmm.version.version}")

    atoms = read_pdb_atoms(MINIMIZED_PDB)
    psf = CharmmPsfFile(str(CLEANED_PSF))
    pdb = PDBFile(str(MINIMIZED_PDB))
    if len(atoms) != len(pdb.positions):
        raise ValueError(f"PDB atom count mismatch: parsed {len(atoms)} atoms, OpenMM positions {len(pdb.positions)}")
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
    strong_count, ligand_restraint_count, unrestrained_count = add_position_restraints(system, atoms, pdb.positions)
    print(f"Strong protein/lipid/cholesterol restraints: {strong_count}")
    print(f"Ligand heavy-atom restraints: {ligand_restraint_count}")
    print(f"Unrestrained atoms: {unrestrained_count}")
    print(f"Timestep: {TIMESTEP_PS * 1000.0:.3f} fs")
    print("Heating schedule:")
    for stage in HEATING_STAGES:
        print(f"- {stage.name}: {stage.temperature_k:.1f} K for {stage.steps} steps")

    integrator = LangevinMiddleIntegrator(50.0 * kelvin, 1.0 / picoseconds, TIMESTEP_PS * picoseconds)
    integrator.setRandomNumberSeed(140220)
    simulation = Simulation(psf.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    simulation.context.applyConstraints(1e-8)
    simulation.context.setVelocitiesToTemperature(50.0 * kelvin, 140220)
    simulation.reporters.append(DCDReporter(str(DCD_PATH), DCD_INTERVAL_STEPS))

    reference_coords_a = positions_to_angstrom(pdb.positions)
    last_finite_positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=False)
    last_finite_step = 0
    max_temperature = float("-inf")
    has_nan = False
    potential_energy_finite = True
    stopped_early = False
    stop_reason = ""
    last_row: dict[str, Any] | None = None
    last_finite_row: dict[str, Any] | None = None

    def empty_metrics() -> dict[str, Any]:
        return {
            "key_min_distance_a": float("nan"),
            "near_chain_count": 0,
            "ligand_rmsd_a": float("nan"),
            "key_distances": [],
            "protein_chain_distances": {},
        }

    def make_row(
        stage_name: str,
        step: int,
        potential: float,
        kinetic: float,
        total: float,
        temperature: float,
        metrics: dict[str, Any],
        status: str,
    ) -> dict[str, Any]:
        return {
            "stage": stage_name,
            "step": step,
            "temperature_K": f"{temperature:.6g}" if math.isfinite(temperature) else "nan",
            "potential_energy_kJ_mol": f"{potential:.12g}" if math.isfinite(potential) else "nan",
            "total_energy_kJ_mol": f"{total:.12g}" if math.isfinite(total) else "nan",
            "ligand_min_distance_to_key_residues_A": f"{metrics['key_min_distance_a']:.6g}" if math.isfinite(metrics["key_min_distance_a"]) else "nan",
            "ligand_near_chains": metrics["near_chain_count"],
            "ligand_rmsd_A": f"{metrics['ligand_rmsd_a']:.6g}" if math.isfinite(metrics["ligand_rmsd_a"]) else "nan",
            "status": status,
        }

    def record_sample(
        writer: csv.DictWriter,
        stage_name: str,
        step: int,
        potential: float,
        kinetic: float,
        total: float,
        temperature: float,
        positions: Any,
        force_write: bool,
    ) -> tuple[dict[str, Any], bool, bool]:
        nonlocal last_finite_positions, last_finite_step, last_finite_row, max_temperature, has_nan, potential_energy_finite

        max_temperature = max(max_temperature, temperature if math.isfinite(temperature) else float("-inf"))
        potential_energy_finite = potential_energy_finite and math.isfinite(potential)
        current_coords_a: list[tuple[float, float, float]] = []
        positions_finite = False
        if sample_values_are_finite(potential, kinetic, total, temperature):
            current_coords_a = positions_to_angstrom(positions)
            positions_finite = all(math.isfinite(value) for xyz in current_coords_a for value in xyz)

        if positions_finite:
            metrics = retention_metrics(atoms, reference_coords_a, current_coords_a)
            last_finite_positions = positions
            last_finite_step = step
        else:
            metrics = empty_metrics()

        status, abort, nonfinite = status_for_sample(
            potential,
            kinetic,
            total,
            temperature,
            metrics["key_min_distance_a"],
            metrics["near_chain_count"],
        )
        nonfinite = nonfinite or not positions_finite
        has_nan = has_nan or nonfinite
        if nonfinite and "nonfinite_energy_or_temperature" not in status:
            status = status + ";nonfinite_positions"
            abort = True

        row = make_row(stage_name, step, potential, kinetic, total, temperature, metrics, status)
        if positions_finite:
            last_finite_row = row
        if force_write or abort:
            writer.writerow(row)
            print(
                "step {step} stage {stage} PE {pe:.6g} KE {ke:.6g} Total {total:.6g} "
                "Temp {temp:.3f} K key_min {key:.3f} A near_chains {chains} ligand_RMSD {rmsd:.3f} A status {status}".format(
                    step=step,
                    stage=stage_name,
                    pe=potential,
                    ke=kinetic,
                    total=total,
                    temp=temperature,
                    key=metrics["key_min_distance_a"],
                    chains=metrics["near_chain_count"],
                    rmsd=metrics["ligand_rmsd_a"],
                    status=status,
                )
            )
        return row, abort, nonfinite

    with RETENTION_CSV.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=RETENTION_FIELDS)
        writer.writeheader()

        global_step = 0
        initial_state = simulation.context.getState(getEnergy=True, getPositions=True)
        initial_potential = initial_state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        initial_kinetic = initial_state.getKineticEnergy().value_in_unit(kilojoules_per_mole)
        initial_total = initial_potential + initial_kinetic
        initial_temperature = estimate_temperature_k(initial_kinetic, system)
        last_row, abort, _ = record_sample(
            writer,
            "initial_50K",
            0,
            initial_potential,
            initial_kinetic,
            initial_total,
            initial_temperature,
            initial_state.getPositions(asNumpy=False),
            True,
        )
        f.flush()
        if abort:
            print(f"Initial sample warning before MD steps: {last_row['status']}")

        for stage in HEATING_STAGES:
            if stopped_early:
                break
            integrator.setTemperature(stage.temperature_k * kelvin)
            print(f"Starting stage {stage.name}")
            remaining = stage.steps
            while remaining > 0:
                chunk_size = 1 if global_step < REPORT_INTERVAL_STEPS else SAFETY_CHECK_INTERVAL_STEPS
                chunk = min(chunk_size, remaining)
                simulation.step(chunk)
                global_step += chunk
                remaining -= chunk
                state = simulation.context.getState(getEnergy=True, getPositions=True)
                potential = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
                kinetic = state.getKineticEnergy().value_in_unit(kilojoules_per_mole)
                total = potential + kinetic
                temperature = estimate_temperature_k(kinetic, system)
                current_positions = state.getPositions(asNumpy=False)
                force_write = global_step % REPORT_INTERVAL_STEPS == 0 or remaining == 0
                last_row, abort, _ = record_sample(
                    writer,
                    stage.name,
                    global_step,
                    potential,
                    kinetic,
                    total,
                    temperature,
                    current_positions,
                    force_write,
                )
                f.flush()
                if abort:
                    stopped_early = True
                    stop_reason = last_row["status"]
                    print(f"Safety stop triggered at step {global_step}: {stop_reason}")
                    break
            if stopped_early:
                break

    write_pdb_with_positions(
        MINIMIZED_PDB,
        last_finite_positions,
        FINAL_PDB,
        [
            f"REMARK last finite snapshot step {last_finite_step}",
            f"REMARK safety stop reason {stop_reason or 'none'}",
        ],
    )
    final_coords_a = positions_to_angstrom(last_finite_positions)
    final_metrics = retention_metrics(atoms, reference_coords_a, final_coords_a)
    completed = (not stopped_early) and last_row is not None and int(last_row["step"]) == sum(stage.steps for stage in HEATING_STAGES)
    max_temperature = max_temperature if math.isfinite(max_temperature) else float("nan")
    recommendation = recommend_longer_restrained_md(
        completed=completed,
        has_nan=has_nan,
        max_temperature_k=max_temperature,
        final_key_min_distance_a=final_metrics["key_min_distance_a"],
        final_near_chains=final_metrics["near_chain_count"],
        final_ligand_rmsd_a=final_metrics["ligand_rmsd_a"],
    )
    result = {
        "completed": completed,
        "last_step": int(last_row["step"]) if last_row is not None else 0,
        "has_nan": has_nan,
        "temperature_stable": math.isfinite(max_temperature) and max_temperature <= TEMPERATURE_ABORT_K,
        "max_temperature_k": max_temperature,
        "potential_energy_finite": potential_energy_finite,
        "last_potential_energy_kj_mol": float(last_row["potential_energy_kJ_mol"]) if last_row is not None else float("nan"),
        "last_total_energy_kj_mol": float(last_row["total_energy_kJ_mol"]) if last_row is not None else float("nan"),
        "last_temperature_k": float(last_row["temperature_K"]) if last_row is not None else float("nan"),
        "last_finite_step": last_finite_step,
        "last_finite_potential_energy_kj_mol": float(last_finite_row["potential_energy_kJ_mol"]) if last_finite_row is not None else float("nan"),
        "final_metrics": final_metrics,
        "recommend_longer_restrained_md": recommendation,
        "platform_name": platform.getName(),
        "stopped_early": stopped_early,
        "stop_reason": stop_reason,
        "final_pdb": str(FINAL_PDB),
        "dcd": str(DCD_PATH),
        "retention_csv": str(RETENTION_CSV),
        "summary_md": str(SUMMARY_MD),
    }
    write_summary(result)
    print("Short MD completed:", completed)
    print("Final PDB:", FINAL_PDB)
    print("DCD:", DCD_PATH)
    print("Retention CSV:", RETENTION_CSV)
    print("Summary Markdown:", SUMMARY_MD)
    del simulation
    del integrator
    return result


def failure_result(error: str) -> dict[str, Any]:
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    if not RETENTION_CSV.exists():
        with RETENTION_CSV.open("w", newline="", encoding="utf-8-sig") as f:
            csv.DictWriter(f, fieldnames=RETENTION_FIELDS).writeheader()
    return {
        "completed": False,
        "last_step": 0,
        "has_nan": False,
        "temperature_stable": False,
        "max_temperature_k": float("nan"),
        "potential_energy_finite": False,
        "last_potential_energy_kj_mol": float("nan"),
        "last_total_energy_kj_mol": float("nan"),
        "last_temperature_k": float("nan"),
        "final_metrics": {},
        "recommend_longer_restrained_md": False,
        "platform_name": "not reached",
        "error": error,
    }


def print_terminal_summary(result: dict[str, Any]) -> None:
    final_metrics = result.get("final_metrics") or {}
    print("")
    print("Terminal summary")
    print("short MD completed:", result.get("completed"))
    print("last potential energy (kJ/mol):", format_float(result.get("last_potential_energy_kj_mol")))
    print("last temperature (K):", format_float(result.get("last_temperature_k")))
    print("ligand final min distance to key residues (A):", format_float(final_metrics.get("key_min_distance_a")))
    print("ligand final near chains count:", final_metrics.get("near_chain_count", "nan"))
    print("ligand RMSD (A):", format_float(final_metrics.get("ligand_rmsd_a")))
    print("NaN/nonfinite detected:", result.get("has_nan"))
    print("recommend next longer restrained MD:", result.get("recommend_longer_restrained_md"))
    print("report path:", SUMMARY_MD)


def main() -> int:
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    with LOG_PATH.open("w", encoding="utf-8") as log_file:
        tee_out = Tee(sys.__stdout__, log_file)
        tee_err = Tee(sys.__stderr__, log_file)
        with redirect_stdout(tee_out), redirect_stderr(tee_err):
            try:
                result = run_short_md()
            except Exception as exc:
                error = repr(exc)
                print("FAILED: short restrained MD pilot failed:", error)
                traceback.print_exc()
                result = failure_result(error)
                write_summary(result)
            print_terminal_summary(result)
    return 0 if result.get("completed") and not result.get("has_nan") else 1


if __name__ == "__main__":
    raise SystemExit(main())
