from __future__ import annotations

import csv
import math
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Any


PILOT_ROOT = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol")
CLASH_ROOT = PILOT_ROOT / "clash_resolution"
OUTPUTS_DIR = CLASH_ROOT / "outputs"
LOGS_DIR = CLASH_ROOT / "logs"
REPORTS_DIR = CLASH_ROOT / "reports"

CLEANED_PSF = OUTPUTS_DIR / "TRKB_20chol_L002_14.2_bound_cleaned.psf"
CLEANED_PDB = OUTPUTS_DIR / "TRKB_20chol_L002_14.2_bound_cleaned.pdb"
MINIMIZED_PDB = OUTPUTS_DIR / "TRKB_20chol_L002_14.2_minimized.pdb"
RETENTION_CSV = REPORTS_DIR / "minimization_retention_check.csv"
RETENTION_MD = REPORTS_DIR / "minimization_retention_check.md"

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
STRONG_RESTRAINT_K = 1000.0  # kJ/mol/nm^2
WEAK_LIGAND_RESTRAINT_K = 100.0  # kJ/mol/nm^2
MINIMIZATION_MAX_ITERATIONS = 1000
RETENTION_CUTOFF_A = 6.0


@dataclass(frozen=True)
class PdbAtom:
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


def atom_lines(path: Path) -> list[str]:
    return [
        line.rstrip("\n")
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith("ATOM") or line.startswith("HETATM")
    ]


def is_ligand(atom: PdbAtom) -> bool:
    return atom.segid == "LIG" or atom.resname == "L002"


def is_heavy(atom: PdbAtom) -> bool:
    return atom.element != "H"


def is_protein(atom: PdbAtom) -> bool:
    return atom.resname in PROTEIN_RESNAMES or atom.segid.startswith("PRO")


def is_lipid_or_cholesterol(atom: PdbAtom) -> bool:
    return atom.resname in LIPID_RESNAMES or atom.resname in CHOLESTEROL_RESNAMES or atom.segid == "MEMB"


def is_water_or_ion(atom: PdbAtom) -> bool:
    return atom.resname in WATER_RESNAMES or atom.resname in ION_RESNAMES or atom.segid.startswith("WT") or atom.segid.startswith("ION")


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


def choose_platform(platform_cls):
    try:
        return platform_cls.getPlatformByName("OpenCL")
    except Exception:
        for index in range(platform_cls.getNumPlatforms()):
            platform = platform_cls.getPlatform(index)
            if platform.getName() != "CUDA":
                return platform
    raise RuntimeError("No non-CUDA OpenMM platform is available")


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
        return WEAK_LIGAND_RESTRAINT_K
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
    protein_lipid_count = 0
    ligand_count = 0
    other_count = 0
    for index, atom in enumerate(atoms):
        k = restraint_k_for_atom(atom)
        if k is None:
            other_count += 1
            continue
        pos = positions[index].value_in_unit(nanometer)
        force.addParticle(index, [k, pos.x, pos.y, pos.z])
        if is_ligand(atom):
            ligand_count += 1
        else:
            protein_lipid_count += 1
    system.addForce(force)
    return protein_lipid_count, ligand_count, other_count


def distance_a(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def positions_to_angstrom(positions: Any) -> list[tuple[float, float, float]]:
    from openmm.unit import angstrom

    return [(pos.x, pos.y, pos.z) for pos in positions.value_in_unit(angstrom)]


def ligand_heavy_indices(atoms: list[PdbAtom]) -> list[int]:
    return [index for index, atom in enumerate(atoms) if is_ligand(atom) and is_heavy(atom)]


def direct_rmsd_a(indices: list[int], before: list[tuple[float, float, float]], after: list[tuple[float, float, float]]) -> float:
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


def retention_rows(
    atoms: list[PdbAtom],
    before_coords_a: list[tuple[float, float, float]],
    after_coords_a: list[tuple[float, float, float]],
) -> tuple[list[dict[str, Any]], bool, bool, float]:
    ligand_indices = ligand_heavy_indices(atoms)
    residue_order = residue_order_by_protein_segment(atoms)
    key_keys = key_residue_keys(residue_order)
    rows: list[dict[str, Any]] = []
    for segid, resid, label in key_keys:
        residue_indices = [
            index
            for index, atom in enumerate(atoms)
            if atom.segid == segid and atom.resid == resid and is_heavy(atom)
        ]
        min_dist = min_distance_for_indices(ligand_indices, residue_indices, after_coords_a)
        rows.append(
            {
                "key_residue": label,
                "segid": segid,
                "resid": resid,
                "min_distance_A": f"{min_dist:.3f}",
                "within_6A": int(math.isfinite(min_dist) and min_dist <= RETENTION_CUTOFF_A),
            }
        )

    protein_chain_distances: dict[str, float] = {}
    for segid in sorted({atom.segid for atom in atoms if is_protein(atom)}):
        chain_indices = [
            index
            for index, atom in enumerate(atoms)
            if atom.segid == segid and is_protein(atom) and is_heavy(atom)
        ]
        protein_chain_distances[segid] = min_distance_for_indices(ligand_indices, chain_indices, after_coords_a)

    near_chain_count = sum(1 for value in protein_chain_distances.values() if math.isfinite(value) and value <= RETENTION_CUTOFF_A)
    key_retained = any(int(row["within_6A"]) == 1 for row in rows)
    near_two_chains = near_chain_count >= 2
    rmsd = direct_rmsd_a(ligand_indices, before_coords_a, after_coords_a)
    for segid, dist in protein_chain_distances.items():
        rows.append(
            {
                "key_residue": "nearest_protein_chain",
                "segid": segid,
                "resid": "",
                "min_distance_A": f"{dist:.3f}",
                "within_6A": int(math.isfinite(dist) and dist <= RETENTION_CUTOFF_A),
            }
        )
    rows.append(
        {
            "key_residue": "ligand_heavy_rmsd_before_after",
            "segid": "LIG",
            "resid": "1",
            "min_distance_A": f"{rmsd:.3f}",
            "within_6A": "",
        }
    )
    return rows, key_retained, near_two_chains, rmsd


def write_retention_reports(
    rows: list[dict[str, Any]],
    key_retained: bool,
    near_two_chains: bool,
    ligand_rmsd_a: float,
    initial_energy: float,
    final_energy: float,
    finite_energy: bool,
) -> None:
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = ["key_residue", "segid", "resid", "min_distance_A", "within_6A"]
    with RETENTION_CSV.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# Minimization Retention Check",
        "",
        f"- Initial energy: {initial_energy:.6g} kJ/mol",
        f"- Final energy: {final_energy:.6g} kJ/mol",
        f"- Energy finite: {'YES' if finite_energy else 'NO'}",
        f"- Energy decreased: {'YES' if final_energy < initial_energy else 'NO'}",
        f"- Ligand heavy atom RMSD before vs after minimization: {ligand_rmsd_a:.3f} A",
        f"- Ligand within {RETENTION_CUTOFF_A:.1f} A of at least one TYR13/VAL17/SER20 key residue: {'YES' if key_retained else 'NO'}",
        f"- Ligand within {RETENTION_CUTOFF_A:.1f} A of two protein chains: {'YES' if near_two_chains else 'NO'}",
        "",
        "| Target | Segid | Resid | Min distance A | Within 6 A |",
        "|---|---|---|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['key_residue']} | {row['segid']} | {row['resid']} | {row['min_distance_A']} | {row['within_6A']} |"
        )
    RETENTION_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_pdb_with_positions(template_pdb: Path, positions: Any, output_pdb: Path) -> None:
    from openmm.unit import angstrom

    coords = positions.value_in_unit(angstrom)
    atom_line_iter = iter(coords)
    lines: list[str] = [
        "REMARK restrained minimization output for pilot 14.2 / L002 / 20chol",
        "REMARK no MD trajectory was generated",
    ]
    for line in template_pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            pos = next(atom_line_iter)
            lines.append(line[:30] + f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}" + line[54:])
    lines.append("END")
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    output_pdb.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_minimization() -> dict[str, Any]:
    import openmm
    from openmm import LocalEnergyMinimizer, Platform, VerletIntegrator
    from openmm.app import CharmmParameterSet, CharmmPsfFile, PDBFile, PME
    from openmm.unit import kilojoules_per_mole, nanometer, picoseconds

    print("Restrained minimization pilot for 14.2 / L002 / 20chol")
    print("No MD will be run and no trajectory will be written.")
    atoms = read_pdb_atoms(CLEANED_PDB)
    psf = CharmmPsfFile(str(CLEANED_PSF))
    pdb = PDBFile(str(CLEANED_PDB))
    box_vectors = read_gro_box_vectors(GRO_BOX)
    if hasattr(psf, "setBoxVectors"):
        psf.setBoxVectors(*box_vectors)
    else:
        a = box_vectors[0].value_in_unit(nanometer)[0]
        b = box_vectors[1].value_in_unit(nanometer)[1]
        c = box_vectors[2].value_in_unit(nanometer)[2]
        psf.setBox(a * nanometer, b * nanometer, c * nanometer)
    params = CharmmParameterSet(*[str(path) for path in PARAMETER_FILES])
    system = psf.createSystem(
        params,
        nonbondedMethod=PME,
        nonbondedCutoff=1.2 * nanometer,
        switchDistance=1.0 * nanometer,
        constraints=None,
    )
    strong_count, ligand_restraint_count, unrestrained_count = add_position_restraints(system, atoms, pdb.positions)
    print(f"Strong protein/lipid/cholesterol restraints: {strong_count}")
    print(f"Weak ligand heavy-atom restraints: {ligand_restraint_count}")
    print(f"Unrestrained atoms: {unrestrained_count}")
    print("Available OpenMM platforms:")
    for index in range(Platform.getNumPlatforms()):
        print(f"- {Platform.getPlatform(index).getName()}")
    platform = choose_platform(Platform)
    print(f"Selected OpenMM platform: {platform.getName()}")

    integrator = VerletIntegrator(0.001 * picoseconds)
    context = openmm.Context(system, integrator, platform)
    context.setPositions(pdb.positions)
    before_coords_a = positions_to_angstrom(pdb.positions)

    initial_state = context.getState(getEnergy=True)
    initial_energy = initial_state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    print(f"Initial potential energy (kJ/mol): {initial_energy}")
    print(f"Initial potential energy finite: {math.isfinite(initial_energy)}")
    LocalEnergyMinimizer.minimize(
        context,
        tolerance=10.0 * kilojoules_per_mole / nanometer,
        maxIterations=MINIMIZATION_MAX_ITERATIONS,
    )
    final_state = context.getState(getEnergy=True, getPositions=True)
    final_energy = final_state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    final_positions = final_state.getPositions(asNumpy=False)
    after_coords_a = positions_to_angstrom(final_positions)
    finite_energy = math.isfinite(initial_energy) and math.isfinite(final_energy)
    print(f"Final potential energy (kJ/mol): {final_energy}")
    print(f"Final potential energy finite: {math.isfinite(final_energy)}")
    print(f"Energy decreased: {final_energy < initial_energy}")
    write_pdb_with_positions(CLEANED_PDB, final_positions, MINIMIZED_PDB)

    rows, key_retained, near_two_chains, ligand_rmsd_a = retention_rows(atoms, before_coords_a, after_coords_a)
    write_retention_reports(rows, key_retained, near_two_chains, ligand_rmsd_a, initial_energy, final_energy, finite_energy)
    print(f"Ligand heavy atom RMSD before/after (A): {ligand_rmsd_a:.3f}")
    print(f"Ligand retained near key residues within {RETENTION_CUTOFF_A:.1f} A: {key_retained}")
    print(f"Ligand near two protein chains within {RETENTION_CUTOFF_A:.1f} A: {near_two_chains}")
    print(f"Minimized PDB: {MINIMIZED_PDB}")
    print(f"Retention CSV: {RETENTION_CSV}")
    print(f"Retention Markdown: {RETENTION_MD}")
    del context
    del integrator
    return {
        "initial_energy": initial_energy,
        "final_energy": final_energy,
        "finite_energy": finite_energy,
        "energy_decreased": final_energy < initial_energy,
        "key_retained": key_retained,
        "near_two_chains": near_two_chains,
        "ligand_rmsd_a": ligand_rmsd_a,
    }


def main() -> int:
    try:
        result = run_minimization()
    except Exception as exc:
        print(f"FAILED: restrained minimization failed: {exc!r}")
        return 1
    return 0 if result["finite_energy"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
