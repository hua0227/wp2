from __future__ import annotations

import math

from branchB_common import (
    BRANCH_B_MINIMIZED_PDB,
    CLEANED_PSF,
    REPORTS_DIR,
    atom_category,
    build_charmm_system,
    choose_opencl_platform,
    ensure_dirs,
    force_rows,
    fmt,
    is_hydrogen,
    is_ligand,
    read_pdb_atoms,
    write_csv,
)


CSV_OUT = REPORTS_DIR / "branchB_force_diagnostics.csv"
MD_OUT = REPORTS_DIR / "branchB_force_diagnostics.md"


def main() -> int:
    import openmm
    from openmm import Platform, VerletIntegrator
    from openmm.app import HBonds
    from openmm.unit import kilojoules_per_mole, picoseconds

    ensure_dirs()
    atoms = read_pdb_atoms(BRANCH_B_MINIMIZED_PDB)
    psf, pdb, _params, system = build_charmm_system(CLEANED_PSF, BRANCH_B_MINIMIZED_PDB, constraints=HBonds, rigid_water=True)
    platform = choose_opencl_platform(Platform)
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = openmm.Context(system, integrator, platform)
    context.setPositions(pdb.positions)
    state = context.getState(getEnergy=True, getForces=True)
    potential_energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    rows, has_force_nan = force_rows(atoms, state.getForces(asNumpy=False))
    top50 = rows[:50]
    fieldnames = [
        "rank",
        "index",
        "serial",
        "atom_name",
        "resname",
        "segid",
        "resid",
        "element",
        "category",
        "force_magnitude_kJ_mol_nm",
        "force_x_kJ_mol_nm",
        "force_y_kJ_mol_nm",
        "force_z_kJ_mol_nm",
        "is_ligand_H",
    ]
    csv_rows = []
    for rank, row in enumerate(top50, start=1):
        out = {"rank": rank, **row}
        csv_rows.append(out)
    write_csv(CSV_OUT, csv_rows, fieldnames)

    ligand_h_rows = [row for row in top50 if row["is_ligand_H"]]
    top10_ligand_h_count = sum(1 for row in top50[:10] if row["is_ligand_H"])
    top20_ligand_h_count = sum(1 for row in top50[:20] if row["is_ligand_H"])
    highest = top50[0]
    category_counts: dict[str, int] = {}
    for row in top50:
        category_counts[row["category"]] = category_counts.get(row["category"], 0) + 1
    ligand_h_concentrated = top10_ligand_h_count >= 5 or (highest["is_ligand_H"] == 1 and top20_ligand_h_count >= 5)

    lines = [
        "# Branch B Force Diagnostics",
        "",
        "Input coordinates: Branch B minimized PDB. No velocities were set and no dynamics was run.",
        "",
        f"- OpenMM platform: {platform.getName()}",
        f"- Potential energy: {fmt(potential_energy, 12)} kJ/mol",
        f"- Potential energy finite: {'YES' if math.isfinite(potential_energy) else 'NO'}",
        f"- Force NaN/nonfinite present: {'YES' if has_force_nan else 'NO'}",
        f"- Highest-force atom: {highest['atom_name']} {highest['resname']} {highest['segid']}:{highest['resid']} index {highest['index']}",
        f"- Highest-force magnitude: {fmt(highest['force_magnitude_kJ_mol_nm'], 12)} kJ/mol/nm",
        f"- Top-50 category counts: {category_counts}",
        f"- Ligand H atoms in top 10 / top 20 / top 50: {top10_ligand_h_count} / {top20_ligand_h_count} / {len(ligand_h_rows)}",
        f"- Highest forces concentrated in ligand H atoms: {'YES' if ligand_h_concentrated else 'NO'}",
        "",
        "## Top 50 Force Atoms",
        "",
        "| Rank | Atom | Residue | Category | Force kJ/mol/nm | Ligand H |",
        "|---:|---|---|---|---:|---:|",
    ]
    for rank, row in enumerate(top50, start=1):
        lines.append(
            f"| {rank} | {row['atom_name']} index {row['index']} | {row['resname']} {row['segid']}:{row['resid']} | "
            f"{row['category']} | {fmt(row['force_magnitude_kJ_mol_nm'], 8)} | {row['is_ligand_H']} |"
        )
    if ligand_h_rows:
        lines.extend(["", "## Ligand H Atoms In Top 50", "", "| Rank | Ligand atom | Force kJ/mol/nm |", "|---:|---|---:|"])
        for row in ligand_h_rows:
            rank = top50.index(row) + 1
            lines.append(f"| {rank} | {row['atom_name']} | {fmt(row['force_magnitude_kJ_mol_nm'], 8)} |")
    MD_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"initial_branchB_energy_kJ_mol: {potential_energy:.12g}")
    print(
        "initial_max_force_atom: "
        f"{highest['atom_name']} {highest['resname']} {highest['segid']}:{highest['resid']} "
        f"index={highest['index']} force={highest['force_magnitude_kJ_mol_nm']:.12g}"
    )
    print(f"force_nan: {has_force_nan}")
    print(f"ligand_H_concentrated: {ligand_h_concentrated}")
    print(f"CSV: {CSV_OUT}")
    print(f"Markdown: {MD_OUT}")
    del context
    del integrator
    return 0 if math.isfinite(potential_energy) and not has_force_nan else 1


if __name__ == "__main__":
    raise SystemExit(main())
