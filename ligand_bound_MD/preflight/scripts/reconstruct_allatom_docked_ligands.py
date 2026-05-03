from __future__ import annotations

import csv
import math
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any


TRKB_ROOT = Path(r"C:\TRKB_WP2")
MD_ROOT = TRKB_ROOT / "ligand_bound_MD"
PREFLIGHT_ROOT = MD_ROOT / "preflight"

MOL2_DIR = MD_ROOT / "cgenff_parameterization" / "upload_mol2"
STR_DIR = MD_ROOT / "cgenff_parameterization" / "returned_str"
SELECTED_INPUTS_DIR = MD_ROOT / "selected_inputs"

LIGAND_POSE_ALLATOM_NAMED_DIR = PREFLIGHT_ROOT / "ligand_pose_allatom_named"
COMPLEX_PREVIEWS_ALLATOM_DIR = PREFLIGHT_ROOT / "complex_previews_allatom"
REPORTS_DIR = PREFLIGHT_ROOT / "reports"
SCRIPTS_DIR = PREFLIGHT_ROOT / "scripts"
TEMP_DIR = PREFLIGHT_ROOT / "temp_pose_reconstruction"

REPORT_CSV = REPORTS_DIR / "allatom_pose_reconstruction_report.csv"
SUMMARY_MD = REPORTS_DIR / "allatom_pose_reconstruction_summary.md"

APO_PDB_BY_SYSTEM = {
    "20chol": TRKB_ROOT / "TRKB_20CHOL" / "openmm_short_md_output" / "short_md_final.pdb",
    "40chol": TRKB_ROOT / "TRKB_40CHOL" / "openmm_short_md_output" / "short_md_final.pdb",
}

LIGANDS = [
    ("8.1", "L001"),
    ("14.2", "L002"),
    ("19.1", "L003"),
    ("9.2", "L004"),
    ("12.2", "L005"),
    ("2.3", "L006"),
    ("4.3", "L007"),
    ("6.2", "L008"),
    ("8.3", "L009"),
    ("17.2", "L010"),
]
SYSTEMS = ["20chol", "40chol"]

REPORT_FIELDS = [
    "ligand_id",
    "resname",
    "system",
    "method_used",
    "cgenff_mol2_atom_count",
    "cgenff_mol2_heavy_atom_count",
    "docked_pose_atom_count",
    "docked_pose_heavy_atom_count",
    "reconstructed_atom_count",
    "heavy_atom_count_match",
    "allatom_count_match",
    "element_sequence_check",
    "allatom_named_pose_generated",
    "complex_preview_generated",
    "output_allatom_named_pdb",
    "output_complex_preview_pdb",
    "status",
    "notes",
]


Coord = tuple[float, float, float]


@dataclass(frozen=True)
class Mol2Atom:
    index: int
    name: str
    element: str
    x: float
    y: float
    z: float
    atom_type: str


@dataclass(frozen=True)
class PoseAtom:
    serial: int
    name: str
    element: str
    x: float
    y: float
    z: float


@dataclass
class ReconstructionResult:
    method: str
    coords: list[Coord]
    element_sequence_ok: bool
    notes: list[str]


def detect_rdkit() -> tuple[bool, str]:
    try:
        import rdkit  # noqa: F401

        return True, "available"
    except Exception as exc:
        return False, repr(exc)


def detect_openbabel() -> tuple[bool, str]:
    exe = shutil.which("obabel")
    if exe is None:
        return False, "obabel executable not found on PATH"
    try:
        completed = subprocess.run(
            [exe, "-V"],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
    except Exception as exc:
        return False, repr(exc)
    output = (completed.stdout or completed.stderr or "").strip()
    if completed.returncode != 0:
        return False, output or f"obabel -V returned {completed.returncode}"
    return True, output


def normalize_element(raw: str, atom_name: str = "") -> str:
    token = raw.strip()
    if not token:
        token = atom_name.strip()
    token = token.split(".", 1)[0]
    token = token.split()[0] if token.split() else token
    token = "".join(ch for ch in token if ch.isalpha())
    if not token:
        token = "".join(ch for ch in atom_name if ch.isalpha())
    if not token:
        return ""

    upper = token.upper()
    autodock_map = {
        "A": "C",
        "C": "C",
        "N": "N",
        "NA": "N",
        "NS": "N",
        "OA": "O",
        "O": "O",
        "OS": "O",
        "S": "S",
        "SA": "S",
        "P": "P",
        "HD": "H",
        "H": "H",
        "HS": "H",
        "F": "F",
        "CL": "Cl",
        "BR": "Br",
        "I": "I",
    }
    if upper in autodock_map:
        return autodock_map[upper]
    if len(upper) >= 2 and upper[:2] in {"CL", "BR"}:
        return upper[:2].title()
    return upper[0]


def parse_mol2_atoms(path: Path) -> list[Mol2Atom]:
    atoms: list[Mol2Atom] = []
    if not path.exists():
        return atoms

    in_atoms = False
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        upper = line.upper()
        if upper.startswith("@<TRIPOS>"):
            in_atoms = upper == "@<TRIPOS>ATOM"
            continue
        if not in_atoms or not line:
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        atom_type = parts[5]
        atom_name = parts[1]
        atoms.append(
            Mol2Atom(
                index=int(parts[0]),
                name=atom_name,
                element=normalize_element(atom_type, atom_name),
                x=float(parts[2]),
                y=float(parts[3]),
                z=float(parts[4]),
                atom_type=atom_type,
            )
        )
    return atoms


def parse_pdb_atoms(path: Path) -> list[PoseAtom]:
    atoms: list[PoseAtom] = []
    if not path.exists():
        return atoms
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        try:
            serial = int(line[6:11])
            atom_name = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            parts = line.split()
            if len(parts) < 9:
                continue
            serial = int(parts[1])
            atom_name = parts[2]
            x, y, z = float(parts[5]), float(parts[6]), float(parts[7])
        element = normalize_element(line[76:78].strip(), atom_name)
        atoms.append(PoseAtom(serial=serial, name=atom_name, element=element, x=x, y=y, z=z))
    return atoms


def parse_pdbqt_atoms(path: Path) -> list[PoseAtom]:
    atoms: list[PoseAtom] = []
    if not path.exists():
        return atoms
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        try:
            serial = int(line[6:11])
            atom_name = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            parts = line.split()
            if len(parts) < 9:
                continue
            serial = int(parts[1])
            atom_name = parts[2]
            x, y, z = float(parts[5]), float(parts[6]), float(parts[7])
        atom_type = line[77:].strip()
        element = normalize_element(atom_type, atom_name)
        atoms.append(PoseAtom(serial=serial, name=atom_name, element=element, x=x, y=y, z=z))
    return atoms


def count_heavy_atoms(atoms: list[Mol2Atom] | list[PoseAtom]) -> int:
    return sum(1 for atom in atoms if atom.element != "H")


def heavy_atoms(atoms: list[Mol2Atom] | list[PoseAtom]) -> list[Mol2Atom] | list[PoseAtom]:
    return [atom for atom in atoms if atom.element != "H"]


def find_first(folder: Path, pattern: str) -> Path | None:
    matches = sorted(folder.glob(pattern)) if folder.exists() else []
    return matches[0] if matches else None


def find_pose_files(ligand_id: str, system: str) -> tuple[Path | None, Path | None, Path | None]:
    folder = SELECTED_INPUTS_DIR / ligand_id / system
    pose_pdb = find_first(folder, f"candidate_{ligand_id}_{system}_mode*.pdb")
    pose_pdbqt = find_first(folder, f"candidate_{ligand_id}_{system}_mode*.pdbqt")
    original_pdbqt = folder / f"{ligand_id}.pdbqt"
    return pose_pdb, pose_pdbqt, original_pdbqt if original_pdbqt.exists() else None


def rdkit_mol_from_mol2(mol2_path: Path) -> Any:
    from rdkit import Chem

    mol = Chem.MolFromMol2File(str(mol2_path), sanitize=True, removeHs=False)
    if mol is None:
        mol = Chem.MolFromMol2File(str(mol2_path), sanitize=False, removeHs=False)
        if mol is not None:
            Chem.SanitizeMol(mol)
    if mol is None:
        raise ValueError(f"RDKit could not read CGenFF MOL2: {mol2_path}")
    return mol


def rdkit_pose_mol_from_inputs(
    pose_pdb: Path | None,
    pose_pdbqt: Path | None,
    temp_dir: Path,
    openbabel_available: bool,
) -> Any:
    from rdkit import Chem

    if pose_pdb is not None and pose_pdb.exists():
        mol = Chem.MolFromPDBFile(str(pose_pdb), sanitize=True, removeHs=False)
        if mol is not None:
            return mol

    if pose_pdbqt is not None and pose_pdbqt.exists() and openbabel_available:
        temp_dir.mkdir(parents=True, exist_ok=True)
        converted = temp_dir / "pose_for_rdkit.pdb"
        subprocess.run(
            ["obabel", str(pose_pdbqt), "-O", str(converted)],
            check=True,
            capture_output=True,
            text=True,
            timeout=120,
        )
        mol = Chem.MolFromPDBFile(str(converted), sanitize=True, removeHs=False)
        if mol is not None:
            return mol

    raise ValueError("RDKit could not read a docked pose molecule from PDB/PDBQT")


def generic_heavy_copy(mol: Any) -> Any:
    from rdkit import Chem

    editable = Chem.RWMol()
    index_map: dict[int, int] = {}
    for atom in mol.GetAtoms():
        new_atom = Chem.Atom(atom.GetSymbol())
        new_atom.SetNoImplicit(True)
        new_idx = editable.AddAtom(new_atom)
        index_map[atom.GetIdx()] = new_idx
    for bond in mol.GetBonds():
        begin = index_map[bond.GetBeginAtomIdx()]
        end = index_map[bond.GetEndAtomIdx()]
        editable.AddBond(begin, end, Chem.BondType.SINGLE)
    return editable.GetMol()


def choose_rdkit_mapping(
    ref_mol: Any,
    pose_mol: Any,
    original_pdbqt: Path | None,
) -> tuple[dict[int, int], bool, str]:
    from rdkit import Chem

    ref_heavy_mol = Chem.RemoveHs(ref_mol, sanitize=False)
    pose_heavy_mol = Chem.RemoveHs(pose_mol, sanitize=False)
    if ref_heavy_mol.GetNumAtoms() != pose_heavy_mol.GetNumAtoms():
        raise ValueError(
            "RDKit heavy atom count mismatch "
            f"ref={ref_heavy_mol.GetNumAtoms()} pose={pose_heavy_mol.GetNumAtoms()}"
        )

    matches = list(ref_heavy_mol.GetSubstructMatches(pose_heavy_mol, uniquify=False))
    match_kind = "exact_bond_order"
    if not matches:
        ref_generic = generic_heavy_copy(ref_heavy_mol)
        pose_generic = generic_heavy_copy(pose_heavy_mol)
        matches = list(ref_generic.GetSubstructMatches(pose_generic, uniquify=False))
        match_kind = "generic_connectivity"
    if not matches:
        raise ValueError("RDKit could not map docked heavy atoms onto CGenFF heavy atoms")

    original_heavy = parse_pdbqt_atoms(original_pdbqt) if original_pdbqt is not None else []
    original_heavy = [atom for atom in original_heavy if atom.element != "H"]
    ref_conf = ref_heavy_mol.GetConformer()

    def score(match: tuple[int, ...]) -> float:
        if len(original_heavy) != len(match):
            return 0.0
        total = 0.0
        for pose_idx, ref_idx in enumerate(match):
            ref_pos = ref_conf.GetAtomPosition(ref_idx)
            pose_atom = original_heavy[pose_idx]
            total += (ref_pos.x - pose_atom.x) ** 2
            total += (ref_pos.y - pose_atom.y) ** 2
            total += (ref_pos.z - pose_atom.z) ** 2
        return math.sqrt(total / len(match))

    best_match = min(matches, key=score)
    element_ok = all(
        ref_heavy_mol.GetAtomWithIdx(ref_idx).GetSymbol()
        == pose_heavy_mol.GetAtomWithIdx(pose_idx).GetSymbol()
        for pose_idx, ref_idx in enumerate(best_match)
    )
    mapping = {ref_idx: pose_idx for pose_idx, ref_idx in enumerate(best_match)}
    score_note = ""
    if len(original_heavy) == len(best_match):
        score_note = f"; original-pose mapping RMSD={score(best_match):.4f} A"
    return mapping, element_ok, f"{match_kind}; substructure matches={len(matches)}{score_note}"


def position_tuple(point: Any) -> Coord:
    return (float(point.x), float(point.y), float(point.z))


def collect_local_anchors(ref_mol: Any, parent_idx: int, full_idx_to_target: dict[int, Coord]) -> list[int]:
    anchors: list[int] = []
    seen: set[int] = set()

    def add(idx: int) -> None:
        atom = ref_mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() != "H" and idx in full_idx_to_target and idx not in seen:
            seen.add(idx)
            anchors.append(idx)

    add(parent_idx)
    frontier = [parent_idx]
    for _depth in range(3):
        next_frontier: list[int] = []
        for idx in frontier:
            atom = ref_mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr.GetSymbol() != "H":
                    add(nbr_idx)
                    next_frontier.append(nbr_idx)
        frontier = next_frontier
        if len(anchors) >= 4:
            break
    return anchors


def transform_point_from_anchors(
    source_points: list[Coord],
    target_points: list[Coord],
    point: Coord,
    fallback_source: Coord,
    fallback_target: Coord,
) -> Coord:
    if len(source_points) < 3 or len(target_points) < 3:
        return (
            fallback_target[0] + point[0] - fallback_source[0],
            fallback_target[1] + point[1] - fallback_source[1],
            fallback_target[2] + point[2] - fallback_source[2],
        )

    try:
        import numpy as np

        source = np.asarray(source_points, dtype=float)
        target = np.asarray(target_points, dtype=float)
        moving = np.asarray(point, dtype=float)
        source_centroid = source.mean(axis=0)
        target_centroid = target.mean(axis=0)
        source_centered = source - source_centroid
        target_centered = target - target_centroid
        covariance = source_centered.T @ target_centered
        u, _s, vt = np.linalg.svd(covariance)
        rotation = vt.T @ u.T
        if np.linalg.det(rotation) < 0:
            vt[-1, :] *= -1
            rotation = vt.T @ u.T
        shifted = (moving - source_centroid) @ rotation + target_centroid
        return (float(shifted[0]), float(shifted[1]), float(shifted[2]))
    except Exception:
        return (
            fallback_target[0] + point[0] - fallback_source[0],
            fallback_target[1] + point[1] - fallback_source[1],
            fallback_target[2] + point[2] - fallback_source[2],
        )


def place_hydrogen(
    ref_mol: Any,
    atom_idx: int,
    full_idx_to_target: dict[int, Coord],
) -> Coord:
    conf = ref_mol.GetConformer()
    atom = ref_mol.GetAtomWithIdx(atom_idx)
    source = position_tuple(conf.GetAtomPosition(atom_idx))
    heavy_neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetSymbol() != "H"]
    if not heavy_neighbors:
        return source
    parent_idx = heavy_neighbors[0]
    parent_source = position_tuple(conf.GetAtomPosition(parent_idx))
    parent_target = full_idx_to_target[parent_idx]
    anchor_indices = collect_local_anchors(ref_mol, parent_idx, full_idx_to_target)
    source_points = [position_tuple(conf.GetAtomPosition(idx)) for idx in anchor_indices]
    target_points = [full_idx_to_target[idx] for idx in anchor_indices]
    return transform_point_from_anchors(source_points, target_points, source, parent_source, parent_target)


def reconstruct_with_rdkit(
    mol2_path: Path,
    pose_pdb: Path | None,
    pose_pdbqt: Path | None,
    original_pdbqt: Path | None,
    temp_dir: Path,
    openbabel_available: bool,
) -> ReconstructionResult:
    from rdkit import Chem

    ref_mol = rdkit_mol_from_mol2(mol2_path)
    pose_mol = rdkit_pose_mol_from_inputs(pose_pdb, pose_pdbqt, temp_dir, openbabel_available)
    ref_heavy_mol = Chem.RemoveHs(ref_mol, sanitize=False)
    pose_heavy_mol = Chem.RemoveHs(pose_mol, sanitize=False)
    ref_heavy_to_pose_heavy, element_ok, mapping_note = choose_rdkit_mapping(ref_mol, pose_mol, original_pdbqt)

    ref_heavy_to_full: dict[int, int] = {}
    heavy_counter = 0
    for atom in ref_mol.GetAtoms():
        if atom.GetSymbol() != "H":
            ref_heavy_to_full[heavy_counter] = atom.GetIdx()
            heavy_counter += 1

    pose_conf = pose_heavy_mol.GetConformer()
    full_idx_to_target: dict[int, Coord] = {}
    for ref_heavy_idx, pose_heavy_idx in ref_heavy_to_pose_heavy.items():
        full_idx = ref_heavy_to_full[ref_heavy_idx]
        full_idx_to_target[full_idx] = position_tuple(pose_conf.GetAtomPosition(pose_heavy_idx))

    coords: list[Coord | None] = [None] * ref_mol.GetNumAtoms()
    for full_idx, coord in full_idx_to_target.items():
        coords[full_idx] = coord

    for atom in ref_mol.GetAtoms():
        atom_idx = atom.GetIdx()
        if atom.GetSymbol() == "H":
            coords[atom_idx] = place_hydrogen(ref_mol, atom_idx, full_idx_to_target)

    if any(coord is None for coord in coords):
        missing = [str(idx + 1) for idx, coord in enumerate(coords) if coord is None]
        raise ValueError("Missing reconstructed coordinates for atom indices: " + ", ".join(missing))

    notes = [
        mapping_note,
        "heavy atom coordinates taken from docked pose",
        "hydrogen coordinates are preliminary local placements from CGenFF/reference geometry; later minimization/relaxation is still required",
    ]
    if not element_ok:
        notes.append("heavy atom element mismatch in RDKit mapping")
    return ReconstructionResult("rdkit", [coord for coord in coords if coord is not None], element_ok, notes)


def reconstruct_with_openbabel_fallback(
    cgenff_atoms: list[Mol2Atom],
    pose_pdbqt: Path | None,
    pose_pdb: Path | None,
    temp_dir: Path,
) -> ReconstructionResult:
    source = pose_pdbqt if pose_pdbqt is not None and pose_pdbqt.exists() else pose_pdb
    if source is None or not source.exists():
        raise ValueError("Open Babel fallback has no docked pose source file")
    temp_dir.mkdir(parents=True, exist_ok=True)
    temp_mol2 = temp_dir / "temp_docked_allH.mol2"
    subprocess.run(
        ["obabel", str(source), "-O", str(temp_mol2), "-h"],
        check=True,
        capture_output=True,
        text=True,
        timeout=120,
    )
    fallback_atoms = parse_mol2_atoms(temp_mol2)
    if len(fallback_atoms) != len(cgenff_atoms):
        raise ValueError(
            "Open Babel fallback atom count mismatch "
            f"cgenff={len(cgenff_atoms)} obabel={len(fallback_atoms)}"
        )
    cgenff_elements = [atom.element for atom in cgenff_atoms]
    fallback_elements = [atom.element for atom in fallback_atoms]
    element_ok = cgenff_elements == fallback_elements
    if not element_ok:
        raise ValueError("Open Babel fallback element sequence is not compatible with CGenFF MOL2 atom order")
    coords = [(atom.x, atom.y, atom.z) for atom in fallback_atoms]
    return ReconstructionResult(
        "openbabel_fallback",
        coords,
        True,
        [
            "Open Babel generated all-H docked pose; atom names/residue names rewritten from CGenFF MOL2",
            "hydrogen coordinates are preliminary; later minimization/relaxation is still required",
        ],
    )


def format_atom_name(name: str, element: str) -> str:
    clean = name[:4]
    if len(clean) == 4:
        return clean
    if len(element) == 1:
        return f" {clean:<3}"[:4]
    return f"{clean:<4}"[:4]


def build_named_pose_pdb_lines(atoms: list[Mol2Atom], coords: list[Coord], resname: str) -> list[str]:
    lines: list[str] = []
    for serial, (atom, coord) in enumerate(zip(atoms, coords), start=1):
        atom_name = format_atom_name(atom.name, atom.element)
        element = atom.element.rjust(2)[:2]
        lines.append(
            f"HETATM{serial:5d} {atom_name} {resname:>4s}Z{1:4d}    "
            f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
            f"{1.00:6.2f}{0.00:6.2f}          {element}"
        )
    lines.append("END")
    return lines


def write_named_pose(path: Path, atoms: list[Mol2Atom], coords: list[Coord], resname: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = build_named_pose_pdb_lines(atoms, coords, resname)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_complex_preview(apo_pdb: Path, named_pose: Path, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    apo_lines = [
        line
        for line in apo_pdb.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith("ATOM") or line.startswith("HETATM")
    ]
    ligand_lines = [
        line
        for line in named_pose.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith("ATOM") or line.startswith("HETATM")
    ]
    lines = [
        "REMARK visual-only complex preview; not a simulation-ready topology",
        "REMARK apo ATOM/HETATM records plus reconstructed all-atom ligand coordinates only",
        *apo_lines,
        *ligand_lines,
        "END",
    ]
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=REPORT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def process_ligand_system(
    ligand_id: str,
    resname: str,
    system: str,
    rdkit_available: bool,
    openbabel_available: bool,
) -> dict[str, Any]:
    mol2_path = MOL2_DIR / f"{resname}_{ligand_id}.mol2"
    str_path = STR_DIR / f"{resname}_{ligand_id}.str"
    pose_pdb, pose_pdbqt, original_pdbqt = find_pose_files(ligand_id, system)
    pose_atoms = parse_pdbqt_atoms(pose_pdbqt) if pose_pdbqt is not None else []
    if not pose_atoms and pose_pdb is not None:
        pose_atoms = parse_pdb_atoms(pose_pdb)
    cgenff_atoms = parse_mol2_atoms(mol2_path)

    cgenff_atom_count = len(cgenff_atoms)
    cgenff_heavy_count = count_heavy_atoms(cgenff_atoms)
    docked_atom_count = len(pose_atoms)
    docked_heavy_count = count_heavy_atoms(pose_atoms)
    heavy_match = cgenff_heavy_count > 0 and cgenff_heavy_count == docked_heavy_count

    notes: list[str] = []
    method = ""
    status = "failed"
    reconstructed_atom_count = 0
    allatom_match = False
    element_ok = False
    allatom_generated = False
    preview_generated = False
    allatom_path = (
        LIGAND_POSE_ALLATOM_NAMED_DIR
        / ligand_id
        / system
        / f"{resname}_{ligand_id}_{system}_allatom_named.pdb"
    )
    preview_path = (
        COMPLEX_PREVIEWS_ALLATOM_DIR
        / ligand_id
        / system
        / f"complex_preview_allatom_{ligand_id}_{system}.pdb"
    )

    if not mol2_path.exists():
        notes.append(f"missing CGenFF MOL2: {mol2_path}")
    if not str_path.exists():
        notes.append(f"missing CGenFF STR: {str_path}")
    if pose_pdb is None and pose_pdbqt is None:
        notes.append("missing docked pose PDB/PDBQT")
    if not heavy_match:
        notes.append(
            "heavy atom count mismatch; all-atom pose not generated "
            f"(CGenFF heavy={cgenff_heavy_count}, docked heavy={docked_heavy_count})"
        )

    if mol2_path.exists() and str_path.exists() and (pose_pdb is not None or pose_pdbqt is not None) and heavy_match:
        temp_dir = TEMP_DIR / ligand_id / system
        try:
            if rdkit_available:
                result = reconstruct_with_rdkit(
                    mol2_path,
                    pose_pdb,
                    pose_pdbqt,
                    original_pdbqt,
                    temp_dir,
                    openbabel_available,
                )
            elif openbabel_available:
                result = reconstruct_with_openbabel_fallback(cgenff_atoms, pose_pdbqt, pose_pdb, temp_dir)
            else:
                raise ValueError("neither RDKit nor Open Babel is available")
        except Exception as rdkit_exc:
            if rdkit_available and openbabel_available:
                notes.append(f"RDKit reconstruction failed: {rdkit_exc!r}")
                try:
                    result = reconstruct_with_openbabel_fallback(cgenff_atoms, pose_pdbqt, pose_pdb, temp_dir)
                except Exception as obabel_exc:
                    result = None
                    notes.append(f"Open Babel fallback failed: {obabel_exc!r}")
            else:
                result = None
                notes.append(f"reconstruction failed: {rdkit_exc!r}")

        if result is not None:
            method = result.method
            notes.extend(result.notes)
            reconstructed_atom_count = len(result.coords)
            allatom_match = reconstructed_atom_count == cgenff_atom_count
            element_ok = result.element_sequence_ok and allatom_match
            if allatom_match and element_ok:
                write_named_pose(allatom_path, cgenff_atoms, result.coords, resname)
                allatom_generated = allatom_path.exists()
                apo_pdb = APO_PDB_BY_SYSTEM[system]
                if apo_pdb.exists() and allatom_generated:
                    write_complex_preview(apo_pdb, allatom_path, preview_path)
                    preview_generated = preview_path.exists()
                elif not apo_pdb.exists():
                    notes.append(f"missing apo PDB for {system}: {apo_pdb}")
                status = "ok" if allatom_generated and preview_generated else "failed"
            else:
                notes.append(
                    "reconstructed atom count or element sequence check failed; named pose not written"
                )

    return {
        "ligand_id": ligand_id,
        "resname": resname,
        "system": system,
        "method_used": method,
        "cgenff_mol2_atom_count": cgenff_atom_count,
        "cgenff_mol2_heavy_atom_count": cgenff_heavy_count,
        "docked_pose_atom_count": docked_atom_count,
        "docked_pose_heavy_atom_count": docked_heavy_count,
        "reconstructed_atom_count": reconstructed_atom_count,
        "heavy_atom_count_match": int(heavy_match),
        "allatom_count_match": int(allatom_match),
        "element_sequence_check": int(element_ok),
        "allatom_named_pose_generated": int(allatom_generated),
        "complex_preview_generated": int(preview_generated),
        "output_allatom_named_pdb": str(allatom_path) if allatom_generated else "",
        "output_complex_preview_pdb": str(preview_path) if preview_generated else "",
        "status": status,
        "notes": "; ".join(notes),
    }


def write_summary(
    rows: list[dict[str, Any]],
    rdkit_available: bool,
    rdkit_message: str,
    openbabel_available: bool,
    openbabel_message: str,
) -> None:
    total = len(rows)
    pose_success = [row for row in rows if int(row["allatom_named_pose_generated"]) == 1]
    preview_success = [row for row in rows if int(row["complex_preview_generated"]) == 1]
    failures = [row for row in rows if row["status"] != "ok"]
    method_counts: dict[str, int] = {}
    for row in pose_success:
        method = str(row["method_used"] or "unknown")
        method_counts[method] = method_counts.get(method, 0) + 1

    generated_preview_lines = [
        f"- {row['ligand_id']} {row['system']}: {row['output_complex_preview_pdb']}"
        for row in preview_success
    ]
    failure_lines = [
        f"- {row['ligand_id']} {row['system']}: {row['notes'] or row['status']}"
        for row in failures
    ]
    if not failure_lines:
        failure_lines = ["- NONE"]

    method_text = ", ".join(f"{method}: {count}" for method, count in sorted(method_counts.items()))
    if not method_text:
        method_text = "none"

    lines = [
        "# All-Atom Docked Ligand Pose Reconstruction Summary",
        "",
        "## Why reconstruction was needed",
        "",
        "The docked poses are not CGenFF all-atom ligand coordinate sets. They are docking outputs with heavy atoms plus limited hydrogens and may use a docking/root-branch atom order. CGenFF MOL2/STR files describe all atoms with CGenFF atom names, so an atom-count mismatch blocks the next ligand-bound topology/read-test step.",
        "",
        "This reconstruction maps docked heavy atom coordinates back onto the CGenFF all-atom molecule and writes preliminary all-atom, CGenFF-named ligand PDB files. Heavy atom coordinates come from the docked pose. Hydrogen coordinates are preliminary local placements from the CGenFF/reference geometry and still require later minimization/relaxation before any production simulation workflow.",
        "",
        "## Tools used",
        "",
        f"- RDKit available: {'YES' if rdkit_available else 'NO'} ({rdkit_message})",
        f"- Open Babel available: {'YES' if openbabel_available else 'NO'} ({openbabel_message})",
        f"- Reconstruction methods used: {method_text}",
        "",
        "## Results",
        "",
        f"- Ligand-system combinations processed: {total}",
        f"- All-atom named poses generated: {len(pose_success)}",
        f"- Visual-only complex previews generated: {len(preview_success)}",
        "",
        "## Failures",
        "",
        *failure_lines,
        "",
        "## Generated complex previews",
        "",
        *(generated_preview_lines if generated_preview_lines else ["- NONE"]),
        "",
        "## Important limitation",
        "",
        "The complex preview PDB files are only for visualization. They are apo ATOM/HETATM records plus reconstructed ligand coordinates, not PSF/PDB/topology builds, not validated OpenMM systems, and not simulation-ready topologies.",
        "",
        "## Next step",
        "",
        "The next step is to attempt ligand-bound topology construction and an OpenMM read test using the CGenFF parameters and these preliminary all-atom named ligand coordinates.",
    ]
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    for directory in [
        LIGAND_POSE_ALLATOM_NAMED_DIR,
        COMPLEX_PREVIEWS_ALLATOM_DIR,
        REPORTS_DIR,
        SCRIPTS_DIR,
        TEMP_DIR,
    ]:
        directory.mkdir(parents=True, exist_ok=True)

    rdkit_available, rdkit_message = detect_rdkit()
    openbabel_available, openbabel_message = detect_openbabel()

    rows: list[dict[str, Any]] = []
    for ligand_id, resname in LIGANDS:
        for system in SYSTEMS:
            rows.append(process_ligand_system(ligand_id, resname, system, rdkit_available, openbabel_available))

    write_csv(REPORT_CSV, rows)
    write_summary(rows, rdkit_available, rdkit_message, openbabel_available, openbabel_message)

    total = len(rows)
    pose_success = sum(1 for row in rows if int(row["allatom_named_pose_generated"]) == 1)
    preview_success = sum(1 for row in rows if int(row["complex_preview_generated"]) == 1)
    failures = [row for row in rows if row["status"] != "ok"]

    print(f"RDKit available: {'YES' if rdkit_available else 'NO'} - {rdkit_message}")
    print(f"Open Babel available: {'YES' if openbabel_available else 'NO'} - {openbabel_message}")
    print(f"Top10 x 2 systems total: {total}")
    print(f"all-atom named pose success count: {pose_success}")
    print(f"complex preview success count: {preview_success}")
    print("failure list:")
    if failures:
        for row in failures:
            print(f"- {row['ligand_id']} {row['system']}: {row['notes'] or row['status']}")
    else:
        print("- NONE")
    print(f"report CSV path: {REPORT_CSV}")
    print(f"summary markdown path: {SUMMARY_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
