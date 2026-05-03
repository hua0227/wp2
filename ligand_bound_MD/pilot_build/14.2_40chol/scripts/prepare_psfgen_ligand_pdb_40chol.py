from __future__ import annotations

from dataclasses import dataclass

from common_40chol import (
    CHAIN_ID,
    INPUT_LIGAND_PDB,
    LIGAND_PSFGEN_READY_PDB,
    LIGAND_RESID,
    LIGAND_STR,
    RESNAME,
    SEGID,
    atom_lines,
    ensure_dirs,
    format_pdb_line,
    normalize_element,
    parse_pdb_xyz,
    parse_str_atom_names,
)


@dataclass(frozen=True)
class LigandInputAtom:
    element: str
    x: float
    y: float
    z: float


def read_ligand_atoms():
    atoms = []
    for line in atom_lines(INPUT_LIGAND_PDB):
        x, y, z = parse_pdb_xyz(line)
        atoms.append(LigandInputAtom(normalize_element(line[76:78].strip(), line[12:16].strip()), x, y, z))
    return atoms


def prepare_ligand_pdb():
    ensure_dirs()
    atoms = read_ligand_atoms()
    str_names = parse_str_atom_names(LIGAND_STR, RESNAME)
    if len(atoms) != len(str_names):
        raise ValueError(f"Atom count mismatch: ligand PDB has {len(atoms)} atoms; STR residue has {len(str_names)} atoms")
    lines = [
        "REMARK psfgen-ready ligand coordinate PDB generated for 14.2 / L002 / 40chol pilot only",
        "REMARK atom names copied from L002_14.2.str atom order; coordinates from 40chol all-atom named ligand pose",
    ]
    for serial, (atom, str_name) in enumerate(zip(atoms, str_names), start=1):
        lines.append(
            format_pdb_line("HETATM", serial, str_name, RESNAME, CHAIN_ID, LIGAND_RESID, atom.x, atom.y, atom.z, SEGID, atom.element)
        )
    lines.append("END")
    LIGAND_PSFGEN_READY_PDB.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return LIGAND_PSFGEN_READY_PDB


def main() -> int:
    out = prepare_ligand_pdb()
    print(f"psfgen-ready ligand PDB: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
