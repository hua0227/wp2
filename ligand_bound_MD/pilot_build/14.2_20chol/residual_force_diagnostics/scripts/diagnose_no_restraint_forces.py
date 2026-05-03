from __future__ import annotations

import csv
import math

from residual_common import (
    CLEANED_PSF,
    REPORTS_DIR,
    STAGE4_PDB,
    build_charmm_system,
    choose_opencl_platform,
    ensure_dirs,
    force_rows,
    fmt,
    is_memb49,
    read_pdb_atoms,
    write_csv,
)


TOP50_CSV = REPORTS_DIR / "no_restraint_force_top50.csv"
SUMMARY_MD = REPORTS_DIR / "no_restraint_force_summary.md"
WITH_RESTRAINT_SUMMARY_CSV = REPORTS_DIR / "force_decomposition_summary.csv"


def read_with_restraint_memb49_max() -> dict[str, str] | None:
    if not WITH_RESTRAINT_SUMMARY_CSV.exists():
        return None
    rows = list(csv.DictReader(WITH_RESTRAINT_SUMMARY_CSV.open(encoding="utf-8-sig")))
    if not rows:
        return None
    return max(rows, key=lambda row: float(row["MEMB49_max_force_kJ_mol_nm"]))


def desc(row: dict[str, object]) -> str:
    return f"{row['atom_name']} {row['resname']} {row['segid']}:{row['resid']} index={row['index']}"


def main() -> int:
    import openmm
    from openmm import Platform, VerletIntegrator
    from openmm.app import HBonds
    from openmm.unit import kilojoules_per_mole, picoseconds

    ensure_dirs()
    atoms = read_pdb_atoms(STAGE4_PDB)
    psf, pdb, _params, system = build_charmm_system(CLEANED_PSF, STAGE4_PDB, constraints=HBonds, rigid_water=True)
    platform = choose_opencl_platform(Platform)
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = openmm.Context(system, integrator, platform)
    context.setPositions(pdb.positions)
    state = context.getState(getEnergy=True, getForces=True)
    energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    rows, has_nan = force_rows(atoms, state.getForces(asNumpy=False))
    top50 = [{"rank": rank, **row} for rank, row in enumerate(rows[:50], start=1)]
    write_csv(TOP50_CSV, top50)
    memb49_rows = [row for row in rows if is_memb49(atoms[int(row["index"])])]
    memb49_max = max(memb49_rows, key=lambda row: row["force_magnitude_kJ_mol_nm"])
    top = top50[0]
    memb49_in_top50 = any(row["is_MEMB49"] == 1 for row in top50)
    with_restraint = read_with_restraint_memb49_max()
    lines = [
        "# No-Restraint Force Diagnostics",
        "",
        "No CustomExternalForce restraints were added. No minimization or dynamics was run.",
        "",
        f"- OpenMM platform: {platform.getName()}",
        f"- Potential energy: {fmt(energy, 12)} kJ/mol",
        f"- Energy finite: {'YES' if math.isfinite(energy) else 'NO'}",
        f"- Force NaN/nonfinite present: {'YES' if has_nan else 'NO'}",
        f"- Top force atom: {desc(top)}",
        f"- Top force magnitude: {fmt(top['force_magnitude_kJ_mol_nm'], 12)} kJ/mol/nm",
        f"- MEMB:49 max force atom: {desc(memb49_max)}",
        f"- MEMB:49 max force magnitude: {fmt(memb49_max['force_magnitude_kJ_mol_nm'], 12)} kJ/mol/nm",
        f"- MEMB:49 appears in top 50: {'YES' if memb49_in_top50 else 'NO'}",
        "",
        "## With-Restraint Comparison",
        "",
    ]
    if with_restraint:
        lines.extend(
            [
                f"- With-restraint dominant MEMB:49 force group: {with_restraint['force_label']} ({with_restraint['force_class']})",
                f"- With-restraint MEMB:49 max force: {with_restraint['MEMB49_max_force_kJ_mol_nm']} kJ/mol/nm",
                f"- With-restraint MEMB:49 max atom: {with_restraint['MEMB49_max_force_atom']}",
            ]
        )
    else:
        lines.append("- With-restraint decomposition CSV was not available.")
    lines.extend(
        [
            "",
            "## Top 10 No-Restraint Force Atoms",
            "",
            "| Rank | Atom | Residue | Category | Force kJ/mol/nm | MEMB49 |",
            "|---:|---|---|---|---:|---:|",
        ]
    )
    for row in top50[:10]:
        lines.append(
            f"| {row['rank']} | {row['atom_name']} index {row['index']} | {row['resname']} {row['segid']}:{row['resid']} | "
            f"{row['category']} | {fmt(row['force_magnitude_kJ_mol_nm'], 8)} | {row['is_MEMB49']} |"
        )
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"no_restraint_energy_kJ_mol: {energy:.12g}")
    print(f"no_restraint_top_force_atom: {desc(top)} force={top['force_magnitude_kJ_mol_nm']:.12g}")
    print(f"MEMB49_max_force_without_restraints_kJ_mol_nm: {memb49_max['force_magnitude_kJ_mol_nm']:.12g}")
    print(f"MEMB49_max_force_without_restraints_atom: {desc(memb49_max)}")
    print(f"MEMB49_in_top50_without_restraints: {memb49_in_top50}")
    print(f"force_nan: {has_nan}")
    print(f"Top50 CSV: {TOP50_CSV}")
    print(f"Summary Markdown: {SUMMARY_MD}")
    del context
    del integrator
    return 0 if math.isfinite(energy) and not has_nan else 1


if __name__ == "__main__":
    raise SystemExit(main())
