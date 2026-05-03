from __future__ import annotations

import csv
import json
import math
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Any


PILOT_ROOT = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol")
ROOT = PILOT_ROOT / "restraint_strategy_fix"
SCRIPTS_DIR = ROOT / "scripts"
LOGS_DIR = ROOT / "logs"
OUTPUTS_DIR = ROOT / "outputs"
REPORTS_DIR = ROOT / "reports"
BRANCH_B_CLEANED_PSF = PILOT_ROOT / "targeted_clash_cleanup" / "outputs" / "branch_B_remove_POPC94" / "branch_B_cleaned.psf"
STAGE4_PDB = PILOT_ROOT / "branchB_relaxation" / "outputs" / "branchB_stage4_minimized.pdb"
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
    "ALA","ARG","ASN","ASP","ASH","CYS","CYM","CYX","GLN","GLU","GLH","GLY","HIS","HSD","HSE","HSP",
    "ILE","LEU","LYS","LYN","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
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
    for path in [ROOT, SCRIPTS_DIR, LOGS_DIR, OUTPUTS_DIR, REPORTS_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def normalize_element(raw: str, atom_name: str = "") -> str:
    token = (raw or atom_name).strip()
    letters = "".join(ch for ch in token if ch.isalpha())
    if not letters:
        return ""
    upper = letters.upper()
    if upper.startswith("CL"): return "Cl"
    if upper.startswith("BR"): return "Br"
    if upper.startswith("NA"): return "Na"
    return upper[0]


def parse_pdb_atom_line(line: str, index: int) -> PdbAtom | None:
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
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
    atoms = []
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


def atom_category(atom: PdbAtom) -> str:
    if is_ligand(atom): return "ligand"
    if is_protein(atom): return "protein"
    if is_cholesterol(atom): return "cholesterol"
    if is_lipid(atom): return "lipid"
    if is_water_or_ion(atom): return "water/ion"
    return "other"


def is_memb49(atom: PdbAtom) -> bool:
    return atom.segid == "MEMB" and atom.resid == "49" and atom.resname == "POPC"


def is_proa437_438(atom: PdbAtom) -> bool:
    return atom.segid == "PROA" and atom.resid in {"437", "438"}


def read_gro_box_lengths_a(path: Path = GRO_BOX) -> tuple[float, float, float]:
    vals = [float(x) for x in path.read_text(encoding="utf-8", errors="replace").splitlines()[-1].split()]
    return vals[0] * 10.0, vals[1] * 10.0, vals[2] * 10.0


def read_gro_box_vectors(path: Path = GRO_BOX) -> Any:
    from openmm import Vec3
    from openmm.unit import nanometer
    vals = [float(x) for x in path.read_text(encoding="utf-8", errors="replace").splitlines()[-1].split()]
    if len(vals) == 3:
        a, b, c = vals
        return (Vec3(a,0,0)*nanometer, Vec3(0,b,0)*nanometer, Vec3(0,0,c)*nanometer)
    v1x,v2y,v3z,v1y,v1z,v2x,v2z,v3x,v3y = vals
    return (Vec3(v1x,v1y,v1z)*nanometer, Vec3(v2x,v2y,v2z)*nanometer, Vec3(v3x,v3y,v3z)*nanometer)


def set_psf_box_from_gro(psf: Any) -> None:
    from openmm.unit import nanometer
    vecs = read_gro_box_vectors()
    if hasattr(psf, "setBoxVectors"):
        psf.setBoxVectors(*vecs)
    else:
        psf.setBox(vecs[0].value_in_unit(nanometer)[0]*nanometer, vecs[1].value_in_unit(nanometer)[1]*nanometer, vecs[2].value_in_unit(nanometer)[2]*nanometer)


def choose_opencl_platform(platform_cls: Any) -> Any:
    names = [platform_cls.getPlatform(i).getName() for i in range(platform_cls.getNumPlatforms())]
    if "OpenCL" not in names:
        raise RuntimeError(f"OpenCL platform required; available: {names}")
    return platform_cls.getPlatformByName("OpenCL")


def build_system(psf_path: Path = BRANCH_B_CLEANED_PSF, pdb_path: Path = STAGE4_PDB, constraints: Any = None, rigid_water: bool = False) -> tuple[Any, Any, Any]:
    from openmm.app import CharmmParameterSet, CharmmPsfFile, PDBFile, PME
    from openmm.unit import nanometer
    psf = CharmmPsfFile(str(psf_path))
    pdb = PDBFile(str(pdb_path))
    set_psf_box_from_gro(psf)
    params = CharmmParameterSet(*[str(p) for p in PARAMETER_FILES])
    kwargs = dict(nonbondedMethod=PME, nonbondedCutoff=1.2*nanometer, switchDistance=1.0*nanometer, constraints=constraints)
    if rigid_water:
        kwargs["rigidWater"] = True
    return psf, pdb, psf.createSystem(params, **kwargs)


def positions_to_angstrom(positions: Any) -> list[tuple[float, float, float]]:
    from openmm.unit import angstrom
    return [(p.x, p.y, p.z) for p in positions.value_in_unit(angstrom)]


def distance_a(a: tuple[float,float,float], b: tuple[float,float,float], box: tuple[float,float,float] | None = None) -> float:
    total = 0.0
    for av, bv, L in zip(a, b, box or (0.0,0.0,0.0)):
        d = av - bv
        if box is not None and L > 0:
            d -= round(d / L) * L
        total += d*d
    return math.sqrt(total)


def coords_have_nan(coords: list[tuple[float,float,float]]) -> bool:
    return any(not math.isfinite(v) for xyz in coords for v in xyz)


def force_rows(atoms: list[PdbAtom], forces: Any) -> tuple[list[dict[str, Any]], bool]:
    from openmm.unit import kilojoules_per_mole, nanometer
    vals = forces.value_in_unit(kilojoules_per_mole / nanometer)
    rows = []
    has_nan = False
    for i, (atom, f) in enumerate(zip(atoms, vals)):
        fx, fy, fz = float(f.x), float(f.y), float(f.z)
        mag = math.sqrt(fx*fx + fy*fy + fz*fz)
        finite = all(math.isfinite(x) for x in [fx, fy, fz, mag])
        has_nan = has_nan or not finite
        rows.append({
            "index": i, "serial": atom.serial, "atom_name": atom.name, "resname": atom.resname,
            "segid": atom.segid, "resid": atom.resid, "element": atom.element, "category": atom_category(atom),
            "force_x_kJ_mol_nm": fx, "force_y_kJ_mol_nm": fy, "force_z_kJ_mol_nm": fz,
            "force_magnitude_kJ_mol_nm": mag,
            "is_MEMB49": int(is_memb49(atom)), "is_PROA437_438": int(is_proa437_438(atom)), "is_ligand": int(is_ligand(atom)),
        })
    rows.sort(key=lambda r: r["force_magnitude_kJ_mol_nm"], reverse=True)
    return rows, has_nan


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def write_json(path: Path, data: dict[str, Any]) -> None:
    def clean(v: Any) -> Any:
        if isinstance(v, dict): return {k: clean(x) for k, x in v.items()}
        if isinstance(v, list): return [clean(x) for x in v]
        if isinstance(v, tuple): return list(v)
        if isinstance(v, float): return v if math.isfinite(v) else None
        return v
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(clean(data), indent=2) + "\n", encoding="utf-8")


def direct_rmsd(indices: list[int], ref: list[tuple[float,float,float]], cur: list[tuple[float,float,float]]) -> float:
    if not indices: return float("nan")
    return math.sqrt(sum(distance_a(ref[i], cur[i])**2 for i in indices) / len(indices))


def retention_and_clashes(atoms: list[PdbAtom], ref_coords: list[tuple[float,float,float]], cur_coords: list[tuple[float,float,float]]) -> tuple[dict[str, Any], dict[str, Any]]:
    box = read_gro_box_lengths_a()
    lig_heavy = [i for i,a in enumerate(atoms) if is_ligand(a) and is_heavy(a)]
    # key residues by known resid in this processed structure
    key_indices = [i for i,a in enumerate(atoms) if a.segid in {"PROA","PROB"} and a.resid in {"434","438","441"} and is_heavy(a)]
    key_min = min(distance_a(cur_coords[i], cur_coords[j], box) for i in lig_heavy for j in key_indices)
    chain_dists = {}
    for seg in sorted({a.segid for a in atoms if is_protein(a)}):
        inds = [i for i,a in enumerate(atoms) if a.segid == seg and is_protein(a) and is_heavy(a)]
        chain_dists[seg] = min(distance_a(cur_coords[i], cur_coords[j], box) for i in lig_heavy for j in inds)
    near_chains = sum(1 for d in chain_dists.values() if d <= RETENTION_CUTOFF_A)
    rmsd = direct_rmsd(lig_heavy, ref_coords, cur_coords)
    lig_all = [i for i,a in enumerate(atoms) if is_ligand(a)]
    other = [i for i,a in enumerate(atoms) if not is_ligand(a)]
    severe = close = 0
    for i in lig_all:
        for j in other:
            d = distance_a(cur_coords[i], cur_coords[j], box)
            if d < CLOSE_CUTOFF_A:
                close += 1
                if d < SEVERE_CUTOFF_A:
                    severe += 1
    return (
        {"key_min_distance_A": key_min, "near_chains": near_chains, "ligand_rmsd_A": rmsd, "chain_distances": chain_dists},
        {"severe_lt_1A": severe, "close_lt_1p5A": close},
    )


def pbc_restraint_force(restraint_map: dict[int, float], positions: Any) -> Any:
    from openmm import CustomExternalForce
    from openmm.unit import nanometer
    force = CustomExternalForce("0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2")
    force.addPerParticleParameter("k")
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    for i, k in restraint_map.items():
        p = positions[i].value_in_unit(nanometer)
        force.addParticle(i, [k, p.x, p.y, p.z])
    return force


def write_pdb(template: Path, coords: list[tuple[float,float,float]], output: Path, remarks: list[str]) -> None:
    it = iter(coords)
    lines = remarks[:]
    for line in template.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM","HETATM")):
            x,y,z = next(it)
            lines.append(line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[54:])
    lines.append("END")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def fmt(v: Any, digits: int = 6) -> str:
    try: n = float(v)
    except Exception: return "nan"
    return f"{n:.{digits}g}" if math.isfinite(n) else "nan"
