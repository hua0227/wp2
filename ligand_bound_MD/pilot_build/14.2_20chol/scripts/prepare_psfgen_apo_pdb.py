from __future__ import annotations

from pathlib import Path


PILOT_ROOT = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol")
OUTPUTS_DIR = PILOT_ROOT / "outputs"

TEMPLATE_PDB = Path(r"C:\TRKB_WP2\TRKB_20CHOL\step5_assembly.pdb")
SHORTMD_PDB = Path(r"C:\TRKB_WP2\TRKB_20CHOL\openmm_short_md_output\short_md_final.pdb")
OUTPUT_PDB = OUTPUTS_DIR / "TRKB_20chol_shortmd_psfgen_ready.pdb"


def atom_lines(path: Path) -> list[str]:
    return [
        line.rstrip("\n")
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith("ATOM") or line.startswith("HETATM")
    ]


def merge_atom_line(template_line: str, coordinate_line: str) -> str:
    return template_line[:30] + coordinate_line[30:54] + template_line[54:]


def prepare_apo_pdb() -> Path:
    template_atoms = atom_lines(TEMPLATE_PDB)
    shortmd_atoms = atom_lines(SHORTMD_PDB)
    if len(template_atoms) != len(shortmd_atoms):
        raise ValueError(
            f"Atom count mismatch: template has {len(template_atoms)} atoms; short-MD PDB has {len(shortmd_atoms)} atoms"
        )

    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    lines = [
        "REMARK psfgen-ready apo coordinate PDB generated for pilot build only",
        "REMARK atom identity/segid/resid columns from step5_assembly.pdb; coordinates from short_md_final.pdb",
    ]
    lines.extend(merge_atom_line(template, coords) for template, coords in zip(template_atoms, shortmd_atoms))
    lines.append("END")
    OUTPUT_PDB.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return OUTPUT_PDB


def main() -> int:
    out = prepare_apo_pdb()
    print(f"psfgen-ready apo PDB: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
