from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Any


PILOT_ROOT = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol")
CLEANUP_ROOT = PILOT_ROOT / "targeted_clash_cleanup"
REPORTS_DIR = CLEANUP_ROOT / "reports"
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
STRONG_RESTRAINT_K = 5000.0  # kJ/mol/nm^2, protein backbone + lipid/cholesterol heavy atoms
LIGAND_HEAVY_RESTRAINT_K = 1000.0  # kJ/mol/nm^2, moderate relative to the 5000 heavy environment restraint
MINIMIZATION_MAX_ITERATIONS = 10000
MINIMIZATION_TOLERANCE = 1.0
RETENTION_CUTOFF_A = 6.0
SEVERE_CUTOFF_A = 1.0
CLOSE_CUTOFF_A = 1.5
TARGET_LIGAND_ATOMS = {"H9", "H11"}


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


def restraint_k_for_atom(atom: PdbAtom) -> float | None:
    if is_ligand(atom) and is_heavy(atom):
        return LIGAND_HEAVY_RESTRAINT_K
    if (is_lipid(atom) or is_cholesterol(atom)) and is_heavy(atom):
        return STRONG_RESTRAINT_K
    if is_protein(atom) and atom.name in PROTEIN_BACKBONE_NAMES:
        return STRONG_RESTRAINT_K
    return None


def add_position_restraints(system: Any, atoms: list[PdbAtom], positions: Any) -> dict[str, int]:
    from openmm import CustomExternalForce
    from openmm.unit import nanometer

    force = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addPerParticleParameter("k")
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    counts = {
        "protein_backbone_restrained": 0,
        "lipid_chol_heavy_restrained": 0,
        "ligand_heavy_restrained": 0,
        "water_ion_unrestrained": 0,
        "other_unrestrained": 0,
    }
    for index, atom in enumerate(atoms):
        k = restraint_k_for_atom(atom)
        if k is None:
            if is_water_or_ion(atom):
                counts["water_ion_unrestrained"] += 1
            else:
                counts["other_unrestrained"] += 1
            continue
        pos = positions[index].value_in_unit(nanometer)
        force.addParticle(index, [k, pos.x, pos.y, pos.z])
        if is_ligand(atom):
            counts["ligand_heavy_restrained"] += 1
        elif is_protein(atom):
            counts["protein_backbone_restrained"] += 1
        else:
            counts["lipid_chol_heavy_restrained"] += 1
    system.addForce(force)
    return counts


def positions_to_angstrom(positions: Any) -> list[tuple[float, float, float]]:
    from openmm.unit import angstrom

    return [(pos.x, pos.y, pos.z) for pos in positions.value_in_unit(angstrom)]


def ligand_heavy_indices(atoms: list[PdbAtom]) -> list[int]:
    return [index for index, atom in enumerate(atoms) if is_ligand(atom) and is_heavy(atom)]


def direct_rmsd_a(indices: list[int], before: list[tuple[float, float, float]], after: list[tuple[float, float, float]]) -> float:
    if not indices:
        return float("nan")
    return math.sqrt(sum(distance_a(before[index], after[index]) ** 2 for index in indices) / len(indices))


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
    residue_order = residue_order_by_protein_segment(atoms)
    key_rows: list[dict[str, Any]] = []
    for segid, resid, label in key_residue_keys(residue_order):
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
    h9_h11_popc_details: list[dict[str, Any]] = []
    severe_sq = SEVERE_CUTOFF_A * SEVERE_CUTOFF_A
    close_sq = CLOSE_CUTOFF_A * CLOSE_CUTOFF_A
    for ligand_index in ligand_indices:
        ligand = atoms[ligand_index]
        lcoord = coords_a[ligand_index]
        for other_index in other_indices:
            other = atoms[other_index]
            ocoord = coords_a[other_index]
            dist = distance_a(lcoord, ocoord, box_lengths_a)
            d2 = dist * dist
            if d2 >= close_sq:
                continue
            close += 1
            is_severe = d2 < severe_sq
            if is_severe:
                severe += 1
            if ligand.name in TARGET_LIGAND_ATOMS and other.resname == "POPC":
                h9_h11_popc_close += 1
                if is_severe:
                    h9_h11_popc_severe += 1
                h9_h11_popc_details.append(
                    {
                        "ligand_atom": ligand.name,
                        "other_atom": other.name,
                        "other_segid": other.segid,
                        "other_resid": other.resid,
                        "distance_a": dist,
                        "severe": is_severe,
                    }
                )
    h9_h11_popc_details.sort(key=lambda row: row["distance_a"])
    return {
        "severe_lt_1A": severe,
        "close_lt_1p5A": close,
        "h9_h11_popc_severe_lt_1A": h9_h11_popc_severe,
        "h9_h11_popc_close_lt_1p5A": h9_h11_popc_close,
        "h9_h11_popc_details": h9_h11_popc_details[:20],
    }


def read_gro_box_vectors(path: Path) -> Any:
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


def read_gro_box_lengths_a(path: Path) -> tuple[float, float, float]:
    values = [float(value) for value in path.read_text(encoding="utf-8", errors="replace").splitlines()[-1].split()]
    if len(values) < 3:
        raise ValueError(f"Unsupported GRO box vector line with {len(values)} values")
    return values[0] * 10.0, values[1] * 10.0, values[2] * 10.0


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
    for image in itertools.product([-1, 0, 1], repeat=3):
        translation = (
            image[0] * box_lengths_a[0],
            image[1] * box_lengths_a[1],
            image[2] * box_lengths_a[2],
        )
        total = 0.0
        for index in ligand_indices:
            current = current_coords_a[index]
            shifted = (
                current[0] + translation[0],
                current[1] + translation[1],
                current[2] + translation[2],
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
            imaged.append(
                (
                    coord[0] + translation_a[0],
                    coord[1] + translation_a[1],
                    coord[2] + translation_a[2],
                )
            )
        else:
            imaged.append(coord)
    return imaged


def write_pdb_with_coords(template_pdb: Path, coords_a: list[tuple[float, float, float]], output_pdb: Path, remarks: list[str]) -> None:
    atom_iter = iter(coords_a)
    lines = ["REMARK deeper restrained minimization output for 14.2 / L002 / 20chol", *remarks]
    for line in template_pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            x, y, z = next(atom_iter)
            lines.append(line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[54:])
    lines.append("END")
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    output_pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: json_safe(subvalue) for key, subvalue in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    if isinstance(value, tuple):
        return [json_safe(item) for item in value]
    if isinstance(value, float):
        if math.isfinite(value):
            return value
        return None
    return value


def write_key_distance_csv(path: Path, branch_id: str, metrics: dict[str, Any]) -> None:
    fieldnames = ["branch", "target", "segid", "resid", "min_distance_A", "within_6A"]
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in metrics["key_rows"]:
            writer.writerow(
                {
                    "branch": branch_id,
                    "target": row["label"],
                    "segid": row["segid"],
                    "resid": row["resid"],
                    "min_distance_A": f"{row['distance_a']:.4f}" if math.isfinite(row["distance_a"]) else "nan",
                    "within_6A": int(row["within_6A"]),
                }
            )
        for segid, dist in metrics["chain_distances"].items():
            writer.writerow(
                {
                    "branch": branch_id,
                    "target": "nearest_protein_chain",
                    "segid": segid,
                    "resid": "",
                    "min_distance_A": f"{dist:.4f}" if math.isfinite(dist) else "nan",
                    "within_6A": int(math.isfinite(dist) and dist <= RETENTION_CUTOFF_A),
                }
            )


def run_minimization(branch_id: str, input_psf: Path, input_pdb: Path, output_pdb: Path, summary_json: Path, key_csv: Path) -> dict[str, Any]:
    import openmm
    from openmm import LocalEnergyMinimizer, Platform, VerletIntegrator
    from openmm.app import CharmmParameterSet, CharmmPsfFile, HBonds, PDBFile, PME
    from openmm.unit import kilojoules_per_mole, nanometer, picoseconds

    print(f"Deeper restrained minimization branch: {branch_id}")
    print("No MD will be run.")
    print(f"Input PSF: {input_psf}")
    print(f"Input PDB: {input_pdb}")
    print(f"Output PDB: {output_pdb}")
    print(f"Max iterations: {MINIMIZATION_MAX_ITERATIONS}")
    print(f"Strong restraints protein backbone/lipid/cholesterol heavy atoms: {STRONG_RESTRAINT_K} kJ/mol/nm^2")
    print(f"Moderate ligand heavy restraints: {LIGAND_HEAVY_RESTRAINT_K} kJ/mol/nm^2")
    print("Water/ion restraints: none")
    print("OpenMM constraints: HBonds, rigidWater=True")
    print("Available OpenMM platforms:")
    for index in range(Platform.getNumPlatforms()):
        print(f"- {Platform.getPlatform(index).getName()}")
    platform = choose_opencl_platform(Platform)
    print(f"Selected OpenMM platform: {platform.getName()}")
    print(f"OpenMM version: {openmm.version.version}")

    atoms = read_pdb_atoms(input_pdb)
    psf = CharmmPsfFile(str(input_psf))
    pdb = PDBFile(str(input_pdb))
    if len(atoms) != len(pdb.positions):
        raise ValueError(f"Atom count mismatch: parsed {len(atoms)} atoms, OpenMM positions {len(pdb.positions)}")
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
    restraint_counts = add_position_restraints(system, atoms, pdb.positions)
    print("Restraint counts:")
    for key, value in restraint_counts.items():
        print(f"- {key}: {value}")

    integrator = VerletIntegrator(0.001 * picoseconds)
    context = openmm.Context(system, integrator, platform)
    context.setPositions(pdb.positions)
    input_coords_a = positions_to_angstrom(pdb.positions)
    box_lengths_a = read_gro_box_lengths_a(GRO_BOX)
    initial_clashes = clash_metrics(atoms, input_coords_a, box_lengths_a)
    raw_initial_state = context.getState(getEnergy=True)
    raw_initial_energy = raw_initial_state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    print(f"Initial input potential energy before constraint projection (kJ/mol): {raw_initial_energy:.12g}")
    context.applyConstraints(1e-8)
    constrained_state = context.getState(getEnergy=True, getPositions=True)
    constrained_energy = constrained_state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    constrained_coords_a = positions_to_angstrom(constrained_state.getPositions(asNumpy=False))
    constrained_clashes = clash_metrics(atoms, constrained_coords_a, box_lengths_a)
    print(f"Initial potential energy after constraint projection (kJ/mol): {constrained_energy:.12g}")
    print(
        "Initial clashes from input coords: severe_lt_1A={severe} close_lt_1p5A={close} H9/H11-POPC severe={hsev} close={hclose}".format(
            severe=initial_clashes["severe_lt_1A"],
            close=initial_clashes["close_lt_1p5A"],
            hsev=initial_clashes["h9_h11_popc_severe_lt_1A"],
            hclose=initial_clashes["h9_h11_popc_close_lt_1p5A"],
        )
    )
    print(
        "Initial clashes after constraints: severe_lt_1A={severe} close_lt_1p5A={close} H9/H11-POPC severe={hsev} close={hclose}".format(
            severe=constrained_clashes["severe_lt_1A"],
            close=constrained_clashes["close_lt_1p5A"],
            hsev=constrained_clashes["h9_h11_popc_severe_lt_1A"],
            hclose=constrained_clashes["h9_h11_popc_close_lt_1p5A"],
        )
    )

    LocalEnergyMinimizer.minimize(
        context,
        tolerance=MINIMIZATION_TOLERANCE * kilojoules_per_mole / nanometer,
        maxIterations=MINIMIZATION_MAX_ITERATIONS,
    )
    final_state = context.getState(getEnergy=True, getPositions=True)
    final_energy = final_state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    final_positions = final_state.getPositions(asNumpy=False)
    final_coords_a = positions_to_angstrom(final_positions)
    image_translation = ligand_image_translation_a(atoms, input_coords_a, final_coords_a, box_lengths_a)
    final_imaged_coords_a = image_ligand_coords(atoms, final_coords_a, image_translation)
    final_clashes = clash_metrics(atoms, final_imaged_coords_a, box_lengths_a)
    retention = retention_metrics(atoms, input_coords_a, final_imaged_coords_a, box_lengths_a)
    finite_energy = math.isfinite(raw_initial_energy) and math.isfinite(constrained_energy) and math.isfinite(final_energy)
    print(f"Final potential energy (kJ/mol): {final_energy:.12g}")
    print(
        "Ligand image translation applied for output/checks (A): "
        f"{image_translation[0]:.4f}, {image_translation[1]:.4f}, {image_translation[2]:.4f}"
    )
    print(f"Energy finite: {finite_energy}")
    print(f"Energy decreased vs constrained initial: {final_energy < constrained_energy}")
    print(
        "Final clashes: severe_lt_1A={severe} close_lt_1p5A={close} H9/H11-POPC severe={hsev} close={hclose}".format(
            severe=final_clashes["severe_lt_1A"],
            close=final_clashes["close_lt_1p5A"],
            hsev=final_clashes["h9_h11_popc_severe_lt_1A"],
            hclose=final_clashes["h9_h11_popc_close_lt_1p5A"],
        )
    )
    print(f"Ligand min distance to TYR13/VAL17/SER20 (A): {retention['key_min_distance_a']:.4f}")
    print(f"Ligand near chains count: {retention['near_chain_count']}")
    print(f"Ligand RMSD relative to pre-minimized pose (A): {retention['ligand_rmsd_a']:.4f}")
    print(f"Ligand retained near key residues: {retention['key_retained']}")
    print(f"Ligand retained near two chains: {retention['near_two_chains']}")

    write_pdb_with_coords(
        input_pdb,
        final_imaged_coords_a,
        output_pdb,
        [
            f"REMARK branch {branch_id}",
            "REMARK no production MD was run",
            f"REMARK final potential energy {final_energy:.12g} kJ/mol",
            "REMARK ligand coordinates imaged to nearest periodic copy of the pre-minimized ligand pose",
            f"REMARK ligand image translation A {image_translation[0]:.4f} {image_translation[1]:.4f} {image_translation[2]:.4f}",
        ],
    )
    write_key_distance_csv(key_csv, branch_id, retention)

    result = {
        "branch": branch_id,
        "input_psf": str(input_psf),
        "input_pdb": str(input_pdb),
        "output_pdb": str(output_pdb),
        "platform": platform.getName(),
        "max_iterations": MINIMIZATION_MAX_ITERATIONS,
        "constraints": "HBonds",
        "rigid_water": True,
        "strong_restraint_k_kj_mol_nm2": STRONG_RESTRAINT_K,
        "ligand_heavy_restraint_k_kj_mol_nm2": LIGAND_HEAVY_RESTRAINT_K,
        "restraint_counts": restraint_counts,
        "initial_input_energy_kj_mol": raw_initial_energy,
        "initial_energy_kj_mol": constrained_energy,
        "final_energy_kj_mol": final_energy,
        "energy_finite": finite_energy,
        "energy_decreased": final_energy < constrained_energy,
        "initial_input_clashes": initial_clashes,
        "initial_constrained_clashes": constrained_clashes,
        "final_clashes": final_clashes,
        "retention": retention,
        "ligand_image_translation_a": image_translation,
    }
    summary_json.parent.mkdir(parents=True, exist_ok=True)
    summary_json.write_text(json.dumps(json_safe(result), indent=2) + "\n", encoding="utf-8")
    print(f"Summary JSON: {summary_json}")
    print(f"Key-distance CSV: {key_csv}")
    print(f"Minimized PDB: {output_pdb}")
    del context
    del integrator
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run pilot-local deeper restrained minimization for one branch.")
    parser.add_argument("--branch-id", required=True)
    parser.add_argument("--input-psf", required=True, type=Path)
    parser.add_argument("--input-pdb", required=True, type=Path)
    parser.add_argument("--output-pdb", required=True, type=Path)
    parser.add_argument("--summary-json", required=True, type=Path)
    parser.add_argument("--key-csv", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = run_minimization(
        branch_id=args.branch_id,
        input_psf=args.input_psf,
        input_pdb=args.input_pdb,
        output_pdb=args.output_pdb,
        summary_json=args.summary_json,
        key_csv=args.key_csv,
    )
    return 0 if result["energy_finite"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
