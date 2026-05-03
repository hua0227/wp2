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
SCRIPTS_DIR = RELAX_ROOT / "scripts"
LOGS_DIR = RELAX_ROOT / "logs"
OUTPUTS_DIR = RELAX_ROOT / "outputs"
REPORTS_DIR = RELAX_ROOT / "reports"

CLEANED_PSF = BRANCH_B_DIR / "branch_B_cleaned.psf"
CLEANED_PDB = BRANCH_B_DIR / "branch_B_cleaned.pdb"
BRANCH_B_MINIMIZED_PDB = BRANCH_B_DIR / "branch_B_minimized.pdb"
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
    for path in [RELAX_ROOT, SCRIPTS_DIR, LOGS_DIR, OUTPUTS_DIR, REPORTS_DIR]:
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


def read_pdb_atom_lines(path: Path) -> list[str]:
    return [
        line
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith("ATOM") or line.startswith("HETATM")
    ]


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
    if atom.resname == "POPC":
        return "POPC"
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


def direct_rmsd_a(indices: list[int], ref: list[tuple[float, float, float]], cur: list[tuple[float, float, float]]) -> float:
    if not indices:
        return float("nan")
    return math.sqrt(sum(distance_a(ref[index], cur[index]) ** 2 for index in indices) / len(indices))


def min_distance_for_indices(
    ligand_indices: list[int],
    other_indices: list[int],
    coords_a: list[tuple[float, float, float]],
    box_lengths_a: tuple[float, float, float] | None = None,
) -> float:
    if not ligand_indices or not other_indices:
        return float("nan")
    best = float("inf")
    for ligand_index in ligand_indices:
        ligand_coord = coords_a[ligand_index]
        for other_index in other_indices:
            best = min(best, distance_a(ligand_coord, coords_a[other_index], box_lengths_a))
    return best


def ligand_heavy_indices(atoms: list[PdbAtom]) -> list[int]:
    return [index for index, atom in enumerate(atoms) if is_ligand(atom) and is_heavy(atom)]


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


def retention_metrics(
    atoms: list[PdbAtom],
    reference_coords_a: list[tuple[float, float, float]],
    current_coords_a: list[tuple[float, float, float]],
    box_lengths_a: tuple[float, float, float] | None = None,
) -> dict[str, Any]:
    ligand_indices = ligand_heavy_indices(atoms)
    key_rows: list[dict[str, Any]] = []
    for segid, resid, label in key_residue_keys(residue_order_by_protein_segment(atoms)):
        residue_indices = [
            index
            for index, atom in enumerate(atoms)
            if atom.segid == segid and atom.resid == resid and is_heavy(atom)
        ]
        min_dist = min_distance_for_indices(ligand_indices, residue_indices, current_coords_a, box_lengths_a)
        key_rows.append(
            {
                "label": label,
                "segid": segid,
                "resid": resid,
                "distance_a": min_dist,
                "within_6A": math.isfinite(min_dist) and min_dist <= RETENTION_CUTOFF_A,
            }
        )

    chain_distances: dict[str, float] = {}
    for segid in sorted({atom.segid for atom in atoms if is_protein(atom)}):
        chain_indices = [
            index
            for index, atom in enumerate(atoms)
            if atom.segid == segid and is_protein(atom) and is_heavy(atom)
        ]
        chain_distances[segid] = min_distance_for_indices(ligand_indices, chain_indices, current_coords_a, box_lengths_a)

    finite_key_distances = [row["distance_a"] for row in key_rows if math.isfinite(row["distance_a"])]
    key_min_distance = min(finite_key_distances) if finite_key_distances else float("nan")
    near_chain_count = sum(1 for value in chain_distances.values() if math.isfinite(value) and value <= RETENTION_CUTOFF_A)
    ligand_rmsd = direct_rmsd_a(ligand_indices, reference_coords_a, current_coords_a)
    return {
        "key_min_distance_a": key_min_distance,
        "near_chain_count": near_chain_count,
        "ligand_rmsd_a": ligand_rmsd,
        "key_rows": key_rows,
        "chain_distances": chain_distances,
        "key_retained": any(row["within_6A"] for row in key_rows),
        "near_two_chains": near_chain_count >= 2,
    }


def clash_metrics(
    atoms: list[PdbAtom],
    coords_a: list[tuple[float, float, float]],
    box_lengths_a: tuple[float, float, float] | None = None,
) -> dict[str, Any]:
    ligand_indices = [index for index, atom in enumerate(atoms) if is_ligand(atom)]
    other_indices = [index for index, atom in enumerate(atoms) if not is_ligand(atom)]
    severe = 0
    close = 0
    h9_h11_popc_severe = 0
    h9_h11_popc_close = 0
    severe_sq = SEVERE_CUTOFF_A * SEVERE_CUTOFF_A
    close_sq = CLOSE_CUTOFF_A * CLOSE_CUTOFF_A
    for ligand_index in ligand_indices:
        ligand = atoms[ligand_index]
        lcoord = coords_a[ligand_index]
        for other_index in other_indices:
            other = atoms[other_index]
            dist = distance_a(lcoord, coords_a[other_index], box_lengths_a)
            d2 = dist * dist
            if d2 >= close_sq:
                continue
            close += 1
            is_severe = d2 < severe_sq
            if is_severe:
                severe += 1
            if ligand.name in {"H9", "H11"} and other.resname == "POPC":
                h9_h11_popc_close += 1
                if is_severe:
                    h9_h11_popc_severe += 1
    return {
        "severe_lt_1A": severe,
        "close_lt_1p5A": close,
        "h9_h11_popc_severe_lt_1A": h9_h11_popc_severe,
        "h9_h11_popc_close_lt_1p5A": h9_h11_popc_close,
    }


def ligand_image_translation_a(
    atoms: list[PdbAtom],
    reference_coords_a: list[tuple[float, float, float]],
    current_coords_a: list[tuple[float, float, float]],
    box_lengths_a: tuple[float, float, float],
) -> tuple[float, float, float]:
    ligand_indices = ligand_heavy_indices(atoms)
    if not ligand_indices:
        return (0.0, 0.0, 0.0)
    best_rmsd = float("inf")
    best_translation = (0.0, 0.0, 0.0)
    for ix in [-1, 0, 1]:
        for iy in [-1, 0, 1]:
            for iz in [-1, 0, 1]:
                translation = (ix * box_lengths_a[0], iy * box_lengths_a[1], iz * box_lengths_a[2])
                total = 0.0
                for index in ligand_indices:
                    shifted = (
                        current_coords_a[index][0] + translation[0],
                        current_coords_a[index][1] + translation[1],
                        current_coords_a[index][2] + translation[2],
                    )
                    total += distance_a(reference_coords_a[index], shifted) ** 2
                rmsd = math.sqrt(total / len(ligand_indices))
                if rmsd < best_rmsd:
                    best_rmsd = rmsd
                    best_translation = translation
    return best_translation


def image_ligand_coords(
    atoms: list[PdbAtom],
    coords_a: list[tuple[float, float, float]],
    translation_a: tuple[float, float, float],
) -> list[tuple[float, float, float]]:
    imaged: list[tuple[float, float, float]] = []
    for atom, coord in zip(atoms, coords_a):
        if is_ligand(atom):
            imaged.append((coord[0] + translation_a[0], coord[1] + translation_a[1], coord[2] + translation_a[2]))
        else:
            imaged.append(coord)
    return imaged


def write_pdb_with_coords(template_pdb: Path, coords_a: list[tuple[float, float, float]], output_pdb: Path, remarks: list[str]) -> None:
    atom_iter = iter(coords_a)
    lines = [*remarks]
    for line in template_pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            x, y, z = next(atom_iter)
            lines.append(line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[54:])
    lines.append("END")
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    output_pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")


def add_position_restraints(
    system: Any,
    atoms: list[PdbAtom],
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


def restraint_map_for_stage(atoms: list[PdbAtom], stage: int) -> dict[int, float]:
    # kJ/mol/nm^2. "Strong" is intentionally larger than ligand restraints used in the failed probe.
    settings = {
        1: {"protein": 5000.0, "environment": 5000.0, "lig_heavy": 5000.0, "lig_h": 5000.0},
        2: {"protein": 5000.0, "environment": 5000.0, "lig_heavy": 1000.0, "lig_h": 100.0},
        3: {"protein": 1000.0, "environment": 1000.0, "lig_heavy": 100.0, "lig_h": 1.0},
        4: {"protein": 1000.0, "environment": 1000.0, "lig_heavy": 1.0, "lig_h": 0.0},
    }[stage]
    mapping: dict[int, float] = {}
    for index, atom in enumerate(atoms):
        k = 0.0
        if is_protein(atom) and atom.name in PROTEIN_BACKBONE_NAMES:
            k = settings["protein"]
        elif (is_lipid(atom) or is_cholesterol(atom)) and is_heavy(atom):
            k = settings["environment"]
        elif is_ligand(atom):
            k = settings["lig_h"] if is_hydrogen(atom) else settings["lig_heavy"]
        if k > 0.0:
            mapping[index] = k
    return mapping


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
                "is_ligand_H": int(is_ligand(atom) and is_hydrogen(atom)),
            }
        )
    rows.sort(key=lambda row: row["force_magnitude_kJ_mol_nm"], reverse=True)
    return rows, has_nan


def max_force_summary(atoms: list[PdbAtom], forces: Any) -> tuple[dict[str, Any], bool]:
    rows, has_nan = force_rows(atoms, forces)
    return rows[0], has_nan


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


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


def write_json(path: Path, data: dict[str, Any]) -> None:
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
