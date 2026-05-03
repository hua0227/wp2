from __future__ import annotations

import csv
import json
import math
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Any


PILOT_ROOT = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol")
BRANCH_B_DIR = PILOT_ROOT / "targeted_clash_cleanup" / "outputs" / "branch_B_remove_POPC94"
RELAX_ROOT = PILOT_ROOT / "branchB_relaxation"
DIAG_ROOT = PILOT_ROOT / "residual_force_diagnostics"
SCRIPTS_DIR = DIAG_ROOT / "scripts"
LOGS_DIR = DIAG_ROOT / "logs"
OUTPUTS_DIR = DIAG_ROOT / "outputs"
REPORTS_DIR = DIAG_ROOT / "reports"

CLEANED_PSF = BRANCH_B_DIR / "branch_B_cleaned.psf"
STAGE4_PDB = RELAX_ROOT / "outputs" / "branchB_stage4_minimized.pdb"
BRANCH_B_INITIAL_MINIMIZED_PDB = BRANCH_B_DIR / "branch_B_minimized.pdb"
ORIGINAL_APO_PSF = Path(r"C:\TRKB_WP2\TRKB_20CHOL\step5_assembly.psf")
ORIGINAL_APO_PDB = Path(r"C:\TRKB_WP2\TRKB_20CHOL\openmm_short_md_output\short_md_final.pdb")
GRO_BOX = Path(r"C:\TRKB_WP2\TRKB_20CHOL\gromacs\step5_input.gro")
TOPPAR_DIR = Path(r"C:\TRKB_WP2\TRKB_20CHOL\toppar")
LIGAND_STR = Path(r"C:\TRKB_WP2\ligand_bound_MD\cgenff_parameterization\returned_str\L002_14.2.str")

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

RETENTION_CUTOFF_A = 6.0
SEVERE_CUTOFF_A = 1.0
CLOSE_CUTOFF_A = 1.5
MEMB49_KEY = ("MEMB", "49", "POPC")


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


def ensure_dirs() -> None:
    for path in [DIAG_ROOT, SCRIPTS_DIR, LOGS_DIR, OUTPUTS_DIR, REPORTS_DIR]:
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
    return PdbAtom(index, serial, name, resname, segid, resid, element, x, y, z)


def read_pdb_atoms(path: Path) -> list[PdbAtom]:
    atoms: list[PdbAtom] = []
    atom_index = 0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        atom = parse_pdb_atom_line(line, atom_index)
        atom_index += 1
        if atom is not None:
            atoms.append(atom)
    return atoms


def is_ligand(atom: PdbAtom) -> bool:
    return atom.segid == "LIG" or atom.resname == "L002"


def is_memb49(atom: PdbAtom) -> bool:
    return (atom.segid, atom.resid, atom.resname) == MEMB49_KEY


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
    return (
        atom.resname in WATER_RESNAMES
        or atom.resname in ION_RESNAMES
        or atom.segid.startswith("WT")
        or atom.segid.startswith("ION")
    )


def atom_category(atom: PdbAtom) -> str:
    if is_ligand(atom):
        return "ligand"
    if is_protein(atom):
        return "protein"
    if is_cholesterol(atom):
        return "cholesterol"
    if is_lipid(atom):
        return "lipid"
    if is_water_or_ion(atom):
        return "water/ion"
    return "other"


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


def choose_opencl_platform(platform_cls: Any) -> Any:
    names = [platform_cls.getPlatform(index).getName() for index in range(platform_cls.getNumPlatforms())]
    if "OpenCL" not in names:
        raise RuntimeError(f"OpenCL platform is required; available platforms: {names}")
    return platform_cls.getPlatformByName("OpenCL")


def build_charmm_system(psf_path: Path, pdb_path: Path, constraints: Any = None, rigid_water: bool = False) -> tuple[Any, Any, Any, Any]:
    from openmm.app import CharmmParameterSet, CharmmPsfFile, PDBFile, PME
    from openmm.unit import nanometer

    psf = CharmmPsfFile(str(psf_path))
    pdb = PDBFile(str(pdb_path))
    set_psf_box_from_gro(psf)
    params = CharmmParameterSet(*[str(path) for path in PARAMETER_FILES])
    kwargs: dict[str, Any] = {
        "nonbondedMethod": PME,
        "nonbondedCutoff": 1.2 * nanometer,
        "switchDistance": 1.0 * nanometer,
        "constraints": constraints,
    }
    if rigid_water:
        kwargs["rigidWater"] = True
    system = psf.createSystem(params, **kwargs)
    return psf, pdb, params, system


def positions_to_angstrom(positions: Any) -> list[tuple[float, float, float]]:
    from openmm.unit import angstrom

    return [(pos.x, pos.y, pos.z) for pos in positions.value_in_unit(angstrom)]


def distance_a(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    box_lengths_a: tuple[float, float, float] | None = None,
) -> float:
    total = 0.0
    for avalue, bvalue, box_length in zip(a, b, box_lengths_a or (0.0, 0.0, 0.0)):
        delta = avalue - bvalue
        if box_lengths_a is not None and box_length > 0.0:
            delta -= round(delta / box_length) * box_length
        total += delta * delta
    return math.sqrt(total)


def coords_have_nan(coords_a: list[tuple[float, float, float]]) -> bool:
    return any(not math.isfinite(value) for xyz in coords_a for value in xyz)


def force_rows(atoms: list[PdbAtom], forces: Any) -> tuple[list[dict[str, Any]], bool]:
    from openmm.unit import kilojoules_per_mole, nanometer

    values = forces.value_in_unit(kilojoules_per_mole / nanometer)
    rows: list[dict[str, Any]] = []
    has_nan = False
    for index, (atom, force) in enumerate(zip(atoms, values)):
        fx, fy, fz = float(force.x), float(force.y), float(force.z)
        mag = math.sqrt(fx * fx + fy * fy + fz * fz)
        finite = math.isfinite(fx) and math.isfinite(fy) and math.isfinite(fz) and math.isfinite(mag)
        has_nan = has_nan or not finite
        rows.append(
            {
                "index": index,
                "serial": atom.serial,
                "atom_name": atom.name,
                "resname": atom.resname,
                "segid": atom.segid,
                "resid": atom.resid,
                "element": atom.element,
                "category": atom_category(atom),
                "force_x_kJ_mol_nm": fx,
                "force_y_kJ_mol_nm": fy,
                "force_z_kJ_mol_nm": fz,
                "force_magnitude_kJ_mol_nm": mag,
                "is_MEMB49": int(is_memb49(atom)),
                "is_ligand": int(is_ligand(atom)),
                "is_ligand_H": int(is_ligand(atom) and is_hydrogen(atom)),
            }
        )
    rows.sort(key=lambda row: row["force_magnitude_kJ_mol_nm"], reverse=True)
    return rows, has_nan


def stage4_restraint_map(atoms: list[PdbAtom]) -> dict[int, float]:
    # Same schedule as Branch B staged minimization Stage 4.
    mapping: dict[int, float] = {}
    for index, atom in enumerate(atoms):
        k = 0.0
        if is_protein(atom) and atom.name in PROTEIN_BACKBONE_NAMES:
            k = 1000.0
        elif (is_lipid(atom) or is_cholesterol(atom)) and is_heavy(atom):
            k = 1000.0
        elif is_ligand(atom) and is_heavy(atom):
            k = 1.0
        if k > 0.0:
            mapping[index] = k
    return mapping


def add_position_restraints(
    system: Any,
    positions: Any,
    restraint_k_by_index: dict[int, float],
    force_group: int | None = None,
) -> int:
    from openmm import CustomExternalForce
    from openmm.unit import nanometer

    force = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addPerParticleParameter("k")
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    if force_group is not None:
        force.setForceGroup(force_group)
    count = 0
    for index, k in restraint_k_by_index.items():
        pos = positions[index].value_in_unit(nanometer)
        force.addParticle(index, [k, pos.x, pos.y, pos.z])
        count += 1
    system.addForce(force)
    return count


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, data: dict[str, Any]) -> None:
    def json_safe(value: Any) -> Any:
        if isinstance(value, dict):
            return {key: json_safe(subvalue) for key, subvalue in value.items()}
        if isinstance(value, list):
            return [json_safe(item) for item in value]
        if isinstance(value, tuple):
            return [json_safe(item) for item in value]
        if isinstance(value, float):
            return value if math.isfinite(value) else None
        return value

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_safe(data), indent=2) + "\n", encoding="utf-8")


def fmt(value: Any, digits: int = 6) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "nan"
    if not math.isfinite(number):
        return "nan"
    return f"{number:.{digits}g}"
