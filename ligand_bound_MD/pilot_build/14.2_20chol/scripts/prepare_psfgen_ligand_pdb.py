from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


PILOT_ROOT = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol")
OUTPUTS_DIR = PILOT_ROOT / "outputs"

INPUT_LIGAND_PDB = Path(
    r"C:\TRKB_WP2\ligand_bound_MD\preflight\ligand_pose_allatom_named\14.2\20chol\L002_14.2_20chol_allatom_named.pdb"
)
LIGAND_STR = Path(r"C:\TRKB_WP2\ligand_bound_MD\cgenff_parameterization\returned_str\L002_14.2.str")
OUTPUT_LIGAND_PDB = OUTPUTS_DIR / "L002_14.2_20chol_psfgen_ready.pdb"

RESNAME = "L002"
CHAIN_ID = "Z"
RESID = 1
SEGID = "LIG"


@dataclass(frozen=True)
class AtomRecord:
    serial: int
    name: str
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
    return upper[0]


def parse_ligand_pdb(path: Path) -> list[AtomRecord]:
    atoms: list[AtomRecord] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        atoms.append(
            AtomRecord(
                serial=int(line[6:11]),
                name=line[12:16].strip(),
                element=normalize_element(line[76:78].strip(), line[12:16].strip()),
                x=float(line[30:38]),
                y=float(line[38:46]),
                z=float(line[46:54]),
            )
        )
    return atoms


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


def format_atom_name(name: str, element: str) -> str:
    clean = name[:4]
    if len(clean) == 4:
        return clean
    if len(element) == 1:
        return f" {clean:<3}"[:4]
    return f"{clean:<4}"[:4]


def format_pdb_line(serial: int, name: str, element: str, x: float, y: float, z: float) -> str:
    atom_name = format_atom_name(name, element)
    return (
        f"HETATM{serial:5d} {atom_name} {RESNAME:>4s}{CHAIN_ID}{RESID:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"{1.00:6.2f}{0.00:6.2f}      {SEGID:<4s}{element:>2s}"
    )


def prepare_psfgen_ready_pdb() -> Path:
    atoms = parse_ligand_pdb(INPUT_LIGAND_PDB)
    str_names = parse_str_atom_names(LIGAND_STR, RESNAME)
    if len(atoms) != len(str_names):
        raise ValueError(
            f"Atom count mismatch: ligand PDB has {len(atoms)} atoms; STR residue has {len(str_names)} atoms"
        )
    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    lines = [
        "REMARK psfgen-ready ligand coordinate PDB generated for pilot build only",
        "REMARK atom names were copied from L002_14.2.str atom order; coordinates are from the reconstructed docked all-atom ligand pose",
    ]
    for serial, (atom, str_name) in enumerate(zip(atoms, str_names), start=1):
        lines.append(format_pdb_line(serial, str_name, atom.element, atom.x, atom.y, atom.z))
    lines.append("END")
    OUTPUT_LIGAND_PDB.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return OUTPUT_LIGAND_PDB


def main() -> int:
    out = prepare_psfgen_ready_pdb()
    print(f"psfgen-ready ligand PDB: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
