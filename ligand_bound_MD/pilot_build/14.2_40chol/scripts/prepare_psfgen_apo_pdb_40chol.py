from __future__ import annotations

from common_40chol import (
    APO_PSF,
    APO_PSFGEN_READY_PDB,
    APO_SHORTMD_PDB,
    APO_TEMPLATE_PDB,
    atom_lines,
    ensure_dirs,
    parse_pdb_xyz,
    parse_psf_atoms,
)


def prepare_apo_pdb():
    ensure_dirs()
    psf_atoms = parse_psf_atoms(APO_PSF)
    template_atoms = atom_lines(APO_TEMPLATE_PDB)
    short_atoms = atom_lines(APO_SHORTMD_PDB)
    if len(psf_atoms) != len(template_atoms) or len(psf_atoms) != len(short_atoms):
        raise ValueError(
            "Atom count mismatch: "
            f"PSF has {len(psf_atoms)} atoms; step5_assembly.pdb has {len(template_atoms)} atoms; "
            f"short_md_final.pdb has {len(short_atoms)} atoms"
        )
    lines = [
        "REMARK psfgen-ready apo coordinate PDB generated for 14.2 / L002 / 40chol pilot only",
        "REMARK atom count validated against step5_assembly.psf",
        "REMARK atom identity/segid/resid columns from step5_assembly.pdb; coordinates from short_md_final.pdb",
    ]
    for template_line, coord_line in zip(template_atoms, short_atoms):
        x, y, z = parse_pdb_xyz(coord_line)
        lines.append(template_line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + template_line[54:])
    lines.append("END")
    APO_PSFGEN_READY_PDB.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return APO_PSFGEN_READY_PDB


def main() -> int:
    out = prepare_apo_pdb()
    print(f"psfgen-ready apo PDB: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
