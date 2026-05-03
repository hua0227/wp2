from __future__ import annotations

import csv
import math
import shutil
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any


PILOT_ROOT = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol")
SCRIPTS_DIR = PILOT_ROOT / "scripts"
LOGS_DIR = PILOT_ROOT / "logs"
OUTPUTS_DIR = PILOT_ROOT / "outputs"
REPORTS_DIR = PILOT_ROOT / "reports"

APO_PSF = Path(r"C:\TRKB_WP2\TRKB_40CHOL\step5_assembly.psf")
APO_TEMPLATE_PDB = Path(r"C:\TRKB_WP2\TRKB_40CHOL\step5_assembly.pdb")
APO_SHORTMD_PDB = Path(r"C:\TRKB_WP2\TRKB_40CHOL\openmm_short_md_output\short_md_final.pdb")
APO_PSFGEN_READY_PDB = OUTPUTS_DIR / "TRKB_40chol_psfgen_ready_apo.pdb"
GRO_BOX = Path(r"C:\TRKB_WP2\TRKB_40CHOL\gromacs\step5_input.gro")
TOPPAR_DIR = Path(r"C:\TRKB_WP2\TRKB_40CHOL\toppar")

LIGAND_STR = Path(r"C:\TRKB_WP2\ligand_bound_MD\cgenff_parameterization\returned_str\L002_14.2.str")
INPUT_LIGAND_PDB = Path(
    r"C:\TRKB_WP2\ligand_bound_MD\preflight\ligand_pose_allatom_named\14.2\40chol\L002_14.2_40chol_allatom_named.pdb"
)
LIGAND_PSFGEN_READY_PDB = OUTPUTS_DIR / "L002_14.2_40chol_psfgen_ready.pdb"

BOUND_PSF = OUTPUTS_DIR / "TRKB_40chol_L002_14.2_bound.psf"
BOUND_PDB = OUTPUTS_DIR / "TRKB_40chol_L002_14.2_bound.pdb"
CLEANED_PSF = OUTPUTS_DIR / "TRKB_40chol_L002_14.2_bound_cleaned.psf"
CLEANED_PDB = OUTPUTS_DIR / "TRKB_40chol_L002_14.2_bound_cleaned.pdb"
MINIMIZED_PDB = OUTPUTS_DIR / "TRKB_40chol_L002_14.2_R2_minimized.pdb"
PROBE_FINAL_PDB = OUTPUTS_DIR / "TRKB_40chol_L002_14.2_very_short_probe_final.pdb"

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

RESNAME = "L002"
SEGID = "LIG"
CHAIN_ID = "Z"
LIGAND_RESID = "1"

SEVERE_CUTOFF_A = 1.0
CLOSE_CUTOFF_A = 1.5
RETENTION_CUTOFF_A = 6.0
MAX_MINIMIZATION_ITERATIONS = 20000

PROTEIN_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYM", "CYX", "GLN", "GLU", "GLH", "GLY", "HIS", "HSD",
    "HSE", "HSP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
LIPID_RESNAMES = {"POP", "POPC", "POPE", "POPS", "POPG", "DOPC", "DPPC", "DMPC"}
CHOLESTEROL_RESNAMES = {"CHL", "CHL1", "CHOL", "CLR"}
WATER_RESNAMES = {"TIP3", "HOH", "WAT", "SOL"}
ION_RESNAMES = {"SOD", "CLA", "POT", "NA", "CL", "K", "CAL", "MG", "ZN"}
PROTEIN_BACKBONE_NAMES = {"N", "CA", "C", "O", "OT1", "OT2"}
KEY_RESIDUES = {"434": "TYR13", "438": "VAL17", "441": "SER20"}


@dataclass(frozen=True)
class PsfAtom:
    index: int
    segid: str
    resid: str
    resname: str
    name: str
    atom_type: str
    charge: float
    mass: float


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
class Contact:
    contact_type: str
    distance: float
    ligand_atom: PdbAtom
    other_atom: PdbAtom
    other_category: str


def ensure_dirs() -> None:
    for path in [SCRIPTS_DIR, LOGS_DIR, OUTPUTS_DIR, REPORTS_DIR]:
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


def infer_element_from_psf(atom_name: str, atom_type: str) -> str:
    return normalize_element(atom_name) or normalize_element(atom_type)


def atom_lines(path: Path) -> list[str]:
    return [
        line.rstrip("\n")
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith(("ATOM", "HETATM"))
    ]


def parse_psf_atoms(path: Path) -> list[PsfAtom]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    start = None
    count = None
    for i, line in enumerate(lines):
        if "!NATOM" in line:
            count = int(line.split()[0])
            start = i + 1
            break
    if start is None or count is None:
        raise ValueError(f"Could not find !NATOM section in {path}")
    atoms: list[PsfAtom] = []
    for raw in lines[start : start + count]:
        parts = raw.split()
        if len(parts) < 8:
            raise ValueError(f"Malformed PSF atom line: {raw}")
        atoms.append(
            PsfAtom(
                index=int(parts[0]),
                segid=parts[1],
                resid=parts[2],
                resname=parts[3],
                name=parts[4],
                atom_type=parts[5],
                charge=float(parts[6]),
                mass=float(parts[7]),
            )
        )
    return atoms


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


def parse_pdb_xyz(line: str) -> tuple[float, float, float]:
    return float(line[30:38]), float(line[38:46]), float(line[46:54])


def format_atom_name(name: str, element: str) -> str:
    clean = name[:4]
    if len(clean) == 4:
        return clean
    if len(element) == 1:
        return f" {clean:<3}"[:4]
    return f"{clean:<4}"[:4]


def format_pdb_line(
    record: str,
    serial: int,
    atom_name: str,
    resname: str,
    chain: str,
    resid: str,
    x: float,
    y: float,
    z: float,
    segid: str,
    element: str,
) -> str:
    atom_field = format_atom_name(atom_name, element)
    resid_field = int(resid) if resid.lstrip("-").isdigit() else 1
    return (
        f"{record:<6}{serial:5d} {atom_field} {resname:>4s}{chain[:1]:1s}{resid_field:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}      {segid:<4s}{element:>2s}"
    )


def parse_str_atom_names(path: Path, resname: str) -> list[str]:
    names: list[str] = []
    in_residue = False
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if line.startswith("RESI"):
            parts = line.split()
            in_residue = len(parts) > 1 and parts[1] == resname
            continue
        if in_residue and line.startswith("ATOM"):
            parts = line.split()
            if len(parts) >= 2:
                names.append(parts[1])
            continue
        if in_residue and line.startswith(("BOND", "DOUBLE", "IMPR", "CMAP", "DONOR", "ACCEPTOR")):
            break
    return names


def read_gro_box_vectors(path: Path = GRO_BOX) -> Any:
    from openmm import Vec3
    from openmm.unit import nanometer

    vals = [float(x) for x in path.read_text(encoding="utf-8", errors="replace").splitlines()[-1].split()]
    if len(vals) == 3:
        a, b, c = vals
        return (Vec3(a, 0, 0) * nanometer, Vec3(0, b, 0) * nanometer, Vec3(0, 0, c) * nanometer)
    if len(vals) == 9:
        v1x, v2y, v3z, v1y, v1z, v2x, v2z, v3x, v3y = vals
        return (Vec3(v1x, v1y, v1z) * nanometer, Vec3(v2x, v2y, v2z) * nanometer, Vec3(v3x, v3y, v3z) * nanometer)
    raise ValueError(f"Unsupported GRO box vector line with {len(vals)} values: {vals}")


def read_gro_box_lengths_a(path: Path = GRO_BOX) -> tuple[float, float, float]:
    vals = [float(x) for x in path.read_text(encoding="utf-8", errors="replace").splitlines()[-1].split()]
    return vals[0] * 10.0, vals[1] * 10.0, vals[2] * 10.0


def set_psf_box_from_gro(psf: Any, gro_path: Path = GRO_BOX) -> None:
    from openmm.unit import nanometer

    vecs = read_gro_box_vectors(gro_path)
    if hasattr(psf, "setBoxVectors"):
        psf.setBoxVectors(*vecs)
    else:
        psf.setBox(
            vecs[0].value_in_unit(nanometer)[0] * nanometer,
            vecs[1].value_in_unit(nanometer)[1] * nanometer,
            vecs[2].value_in_unit(nanometer)[2] * nanometer,
        )


def choose_opencl_platform(platform_cls: Any) -> Any:
    names = [platform_cls.getPlatform(i).getName() for i in range(platform_cls.getNumPlatforms())]
    if "OpenCL" not in names:
        raise RuntimeError(f"OpenCL platform required; available platforms: {names}")
    return platform_cls.getPlatformByName("OpenCL")


def build_system(psf_path: Path, pdb_path: Path, constraints: Any = None, rigid_water: bool = False) -> tuple[Any, Any, Any]:
    from openmm.app import CharmmParameterSet, CharmmPsfFile, PDBFile, PME
    from openmm.unit import nanometer

    psf = CharmmPsfFile(str(psf_path))
    pdb = PDBFile(str(pdb_path))
    set_psf_box_from_gro(psf)
    params = CharmmParameterSet(*[str(path) for path in PARAMETER_FILES])
    kwargs = {
        "nonbondedMethod": PME,
        "nonbondedCutoff": 1.2 * nanometer,
        "switchDistance": 1.0 * nanometer,
        "constraints": constraints,
    }
    if rigid_water:
        kwargs["rigidWater"] = True
    return psf, pdb, psf.createSystem(params, **kwargs)


def positions_to_angstrom(positions: Any) -> list[tuple[float, float, float]]:
    from openmm.unit import angstrom

    return [(float(p.x), float(p.y), float(p.z)) for p in positions.value_in_unit(angstrom)]


def distance_a(a: tuple[float, float, float], b: tuple[float, float, float], box: tuple[float, float, float] | None = None) -> float:
    total = 0.0
    for av, bv, length in zip(a, b, box or (0.0, 0.0, 0.0)):
        delta = av - bv
        if box is not None and length > 0:
            delta -= round(delta / length) * length
        total += delta * delta
    return math.sqrt(total)


def is_ligand(atom: PdbAtom) -> bool:
    return atom.segid == SEGID or atom.resname == RESNAME


def is_hydrogen(atom: PdbAtom) -> bool:
    return atom.element.upper() == "H" or atom.name.upper().startswith("H")


def is_heavy(atom: PdbAtom) -> bool:
    return not is_hydrogen(atom)


def is_protein(atom: PdbAtom) -> bool:
    return atom.resname.upper() in PROTEIN_RESNAMES or atom.segid.upper().startswith("PRO")


def is_cholesterol(atom: PdbAtom) -> bool:
    return atom.resname.upper() in CHOLESTEROL_RESNAMES


def is_lipid(atom: PdbAtom) -> bool:
    return atom.resname.upper() in LIPID_RESNAMES or atom.segid.upper() == "MEMB"


def is_water(atom: PdbAtom) -> bool:
    return atom.resname.upper() in WATER_RESNAMES or atom.segid.upper().startswith("WT")


def is_ion(atom: PdbAtom) -> bool:
    return atom.resname.upper() in ION_RESNAMES or atom.segid.upper().startswith("ION")


def atom_category(atom: PdbAtom) -> str:
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


def find_contacts(atoms: list[PdbAtom]) -> list[Contact]:
    ligand_heavy = [atom for atom in atoms if is_ligand(atom) and is_heavy(atom)]
    other_heavy = [atom for atom in atoms if not is_ligand(atom) and is_heavy(atom)]
    close_sq = CLOSE_CUTOFF_A * CLOSE_CUTOFF_A
    severe_sq = SEVERE_CUTOFF_A * SEVERE_CUTOFF_A
    contacts: list[Contact] = []
    for ligand_atom in ligand_heavy:
        for other_atom in other_heavy:
            dx = ligand_atom.x - other_atom.x
            dy = ligand_atom.y - other_atom.y
            dz = ligand_atom.z - other_atom.z
            d2 = dx * dx + dy * dy + dz * dz
            if d2 < close_sq:
                contacts.append(
                    Contact(
                        contact_type="severe_clash" if d2 < severe_sq else "close_contact",
                        distance=math.sqrt(d2),
                        ligand_atom=ligand_atom,
                        other_atom=other_atom,
                        other_category=atom_category(other_atom),
                    )
                )
    contacts.sort(key=lambda contact: contact.distance)
    return contacts


def molecule_key(atom: PdbAtom) -> tuple[str, str, str, str]:
    return atom_category(atom), atom.segid, atom.resid, atom.resname


def choose_cleanup_targets(contacts: list[Contact]) -> tuple[list[tuple[str, str, str, str]], str, Counter[tuple[str, str, str, str]]]:
    severe = [contact for contact in contacts if contact.contact_type == "severe_clash"]
    counts: Counter[tuple[str, str, str, str]] = Counter(molecule_key(contact.other_atom) for contact in severe)
    if not severe or not counts:
        return [], "no severe clashes; copied bound PSF/PDB as cleaned", counts
    top_key, top_count = counts.most_common(1)[0]
    top_category, top_segid, top_resid, top_resname = top_key
    dominant = top_count >= 3 and top_count / len(severe) >= 0.50
    if top_category in {"water", "ion"}:
        targets = sorted(key for key, _count in counts.items() if key[0] in {"water", "ion"})
        return targets, "deleted all severe-clashing water/ion residues", counts
    if dominant and top_category in {"lipid", "cholesterol"}:
        return [top_key], f"deleted dominant severe-clashing {top_category} molecule {top_resname} {top_segid}:{top_resid}", counts
    return [], "no dominant deletable severe clash source; copied bound PSF/PDB as cleaned", counts


def contact_rows(contacts: list[Contact]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for contact in contacts:
        ligand = contact.ligand_atom
        other = contact.other_atom
        rows.append(
            {
                "contact_type": contact.contact_type,
                "distance_A": f"{contact.distance:.3f}",
                "category": contact.other_category,
                "ligand_serial": ligand.serial,
                "ligand_atom_name": ligand.name,
                "ligand_resname": ligand.resname,
                "ligand_segid": ligand.segid,
                "ligand_resid": ligand.resid,
                "other_serial": other.serial,
                "other_atom_name": other.name,
                "other_resname": other.resname,
                "other_segid": other.segid,
                "other_resid": other.resid,
                "other_element": other.element,
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


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


def copy_bound_to_cleaned() -> None:
    shutil.copyfile(BOUND_PSF, CLEANED_PSF)
    shutil.copyfile(BOUND_PDB, CLEANED_PDB)


def make_selection(targets: list[tuple[str, str, str, str]]) -> str:
    terms = []
    for _category, segid, resid, resname in targets:
        terms.append(f'(segid "{segid}" and resid "{resid}" and resname "{resname}")')
    return " or ".join(terms)


def coords_have_nan(coords: list[tuple[float, float, float]]) -> bool:
    return any(not math.isfinite(v) for xyz in coords for v in xyz)


def force_rows(atoms: list[PdbAtom], forces: Any) -> tuple[list[dict[str, Any]], bool]:
    from openmm.unit import kilojoules_per_mole, nanometer

    vals = forces.value_in_unit(kilojoules_per_mole / nanometer)
    rows: list[dict[str, Any]] = []
    has_nan = False
    for i, (atom, force) in enumerate(zip(atoms, vals)):
        fx, fy, fz = float(force.x), float(force.y), float(force.z)
        mag = math.sqrt(fx * fx + fy * fy + fz * fz)
        finite = all(math.isfinite(x) for x in [fx, fy, fz, mag])
        has_nan = has_nan or not finite
        rows.append(
            {
                "index": i,
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
            }
        )
    rows.sort(key=lambda row: row["force_magnitude_kJ_mol_nm"], reverse=True)
    return rows, has_nan


def pbc_restraint_force(restraint_map: dict[int, float], positions: Any) -> Any:
    from openmm import CustomExternalForce
    from openmm.unit import nanometer

    force = CustomExternalForce("0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2")
    force.addPerParticleParameter("k")
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    for index, k_value in restraint_map.items():
        pos = positions[index].value_in_unit(nanometer)
        force.addParticle(index, [k_value, pos.x, pos.y, pos.z])
    return force


def r2_restraint_map(atoms: list[PdbAtom], protein_k: float = 500.0, ligand_k: float = 25.0, lipid_k: float = 5.0) -> dict[int, float]:
    restraints: dict[int, float] = {}
    for i, atom in enumerate(atoms):
        k = 0.0
        if is_protein(atom) and atom.name in PROTEIN_BACKBONE_NAMES:
            k = protein_k
        elif is_ligand(atom) and is_heavy(atom):
            k = ligand_k
        elif (is_lipid(atom) or is_cholesterol(atom)) and is_heavy(atom):
            k = lipid_k
        if k > 0:
            restraints[i] = k
    return restraints


def ligand_key_distances(atoms: list[PdbAtom], coords: list[tuple[float, float, float]]) -> dict[str, float]:
    box = read_gro_box_lengths_a()
    ligand_heavy = [i for i, atom in enumerate(atoms) if is_ligand(atom) and is_heavy(atom)]
    out: dict[str, float] = {}
    for resid, label in KEY_RESIDUES.items():
        key_indices = [
            i for i, atom in enumerate(atoms)
            if atom.segid in {"PROA", "PROB"} and atom.resid == resid and is_heavy(atom)
        ]
        out[label] = min(distance_a(coords[i], coords[j], box) for i in ligand_heavy for j in key_indices)
    return out


def retention_and_clashes(atoms: list[PdbAtom], ref_coords: list[tuple[float, float, float]], cur_coords: list[tuple[float, float, float]]) -> tuple[dict[str, Any], dict[str, Any]]:
    box = read_gro_box_lengths_a()
    ligand_heavy = [i for i, atom in enumerate(atoms) if is_ligand(atom) and is_heavy(atom)]
    key_indices = [
        i for i, atom in enumerate(atoms)
        if atom.segid in {"PROA", "PROB"} and atom.resid in set(KEY_RESIDUES) and is_heavy(atom)
    ]
    key_min = min(distance_a(cur_coords[i], cur_coords[j], box) for i in ligand_heavy for j in key_indices)
    chain_dists = {}
    for segid in sorted({atom.segid for atom in atoms if is_protein(atom)}):
        indices = [i for i, atom in enumerate(atoms) if atom.segid == segid and is_protein(atom) and is_heavy(atom)]
        chain_dists[segid] = min(distance_a(cur_coords[i], cur_coords[j], box) for i in ligand_heavy for j in indices)
    near_chains = sum(1 for dist in chain_dists.values() if dist <= RETENTION_CUTOFF_A)
    ligand_rmsd = math.sqrt(sum(distance_a(ref_coords[i], cur_coords[i]) ** 2 for i in ligand_heavy) / len(ligand_heavy))
    severe = 0
    close = 0
    other_heavy = [i for i, atom in enumerate(atoms) if not is_ligand(atom) and is_heavy(atom)]
    for i in ligand_heavy:
        for j in other_heavy:
            dist = distance_a(cur_coords[i], cur_coords[j], box)
            if dist < CLOSE_CUTOFF_A:
                close += 1
                if dist < SEVERE_CUTOFF_A:
                    severe += 1
    return (
        {
            "key_min_distance_A": key_min,
            "key_distances_A": ligand_key_distances(atoms, cur_coords),
            "near_chains": near_chains,
            "chain_distances_A": chain_dists,
            "ligand_rmsd_A": ligand_rmsd,
        },
        {"severe_lt_1A": severe, "close_lt_1p5A": close},
    )


def bool_text(value: bool) -> str:
    return "yes" if value else "no"
