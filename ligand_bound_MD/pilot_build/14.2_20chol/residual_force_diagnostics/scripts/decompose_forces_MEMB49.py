from __future__ import annotations

import math

from residual_common import (
    CLEANED_PSF,
    REPORTS_DIR,
    STAGE4_PDB,
    add_position_restraints,
    build_charmm_system,
    choose_opencl_platform,
    distance_a,
    ensure_dirs,
    force_rows,
    fmt,
    is_ligand,
    is_memb49,
    read_gro_box_lengths_a,
    read_pdb_atoms,
    stage4_restraint_map,
    write_csv,
    write_json,
)


SUMMARY_CSV = REPORTS_DIR / "force_decomposition_summary.csv"
MEMB49_CSV = REPORTS_DIR / "force_decomposition_MEMB49.csv"
TOP50_CSV = REPORTS_DIR / "force_decomposition_top50.csv"
SUMMARY_JSON = REPORTS_DIR / "force_decomposition_summary.json"


def force_label(force: object) -> str:
    name = force.__class__.__name__
    if name in {"NonbondedForce", "CustomNonbondedForce"}:
        return "NonbondedForce/CustomNonbondedForce"
    if name == "CustomExternalForce":
        return "CustomExternalForce/restraints"
    if name in {"HarmonicBondForce", "HarmonicAngleForce", "PeriodicTorsionForce", "CMAPTorsionForce"}:
        return name
    return "Other Force"


def top_description(row: dict[str, object]) -> str:
    return f"{row['atom_name']} {row['resname']} {row['segid']}:{row['resid']} index={row['index']}"


def main() -> int:
    import openmm
    from openmm import Platform, VerletIntegrator
    from openmm.app import HBonds
    from openmm.unit import kilojoules_per_mole, picoseconds

    ensure_dirs()
    atoms = read_pdb_atoms(STAGE4_PDB)
    psf, pdb, _params, system = build_charmm_system(CLEANED_PSF, STAGE4_PDB, constraints=HBonds, rigid_water=True)
    for force_index in range(system.getNumForces()):
        system.getForce(force_index).setForceGroup(force_index)
    restraint_group = system.getNumForces()
    restraint_count = add_position_restraints(system, pdb.positions, stage4_restraint_map(atoms), force_group=restraint_group)
    force_metadata: list[dict[str, object]] = []
    for force_index in range(system.getNumForces()):
        force = system.getForce(force_index)
        force_metadata.append(
            {
                "force_index": force_index,
                "force_group": force.getForceGroup(),
                "force_class": force.__class__.__name__,
                "force_label": force_label(force),
            }
        )

    platform = choose_opencl_platform(Platform)
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = openmm.Context(system, integrator, platform)
    context.setPositions(pdb.positions)

    total_state = context.getState(getEnergy=True, getForces=True)
    total_energy = total_state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    total_rows, total_has_nan = force_rows(atoms, total_state.getForces(asNumpy=False))
    top50_rows = []
    for rank, row in enumerate(total_rows[:50], start=1):
        top50_rows.append(
            {
                "rank": rank,
                **row,
            }
        )

    memb49_indices = [index for index, atom in enumerate(atoms) if is_memb49(atom)]
    ligand_indices = [index for index, atom in enumerate(atoms) if is_ligand(atom)]
    coords = [(atom.x, atom.y, atom.z) for atom in atoms]
    box = read_gro_box_lengths_a()
    ligand_memb49_min_distance = min(
        distance_a(coords[i], coords[j], box) for i in memb49_indices for j in ligand_indices
    )

    summary_rows: list[dict[str, object]] = []
    atom_rows: list[dict[str, object]] = []
    group_results: list[dict[str, object]] = []
    for meta in force_metadata:
        group = int(meta["force_group"])
        state = context.getState(getEnergy=True, getForces=True, groups=1 << group)
        energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        rows, has_nan = force_rows(atoms, state.getForces(asNumpy=False))
        max_row = rows[0]
        memb49_force_rows = [row for row in rows if row["index"] in memb49_indices]
        ligand_force_rows = [row for row in rows if row["index"] in ligand_indices]
        memb49_max = max(memb49_force_rows, key=lambda row: row["force_magnitude_kJ_mol_nm"])
        ligand_max = max(ligand_force_rows, key=lambda row: row["force_magnitude_kJ_mol_nm"])
        summary_rows.append(
            {
                "force_group": group,
                "force_index": meta["force_index"],
                "force_class": meta["force_class"],
                "force_label": meta["force_label"],
                "energy_kJ_mol": f"{energy:.12g}",
                "force_nan": int(has_nan),
                "global_max_force_atom": top_description(max_row),
                "global_max_force_kJ_mol_nm": f"{max_row['force_magnitude_kJ_mol_nm']:.12g}",
                "MEMB49_max_force_atom": top_description(memb49_max),
                "MEMB49_max_force_kJ_mol_nm": f"{memb49_max['force_magnitude_kJ_mol_nm']:.12g}",
                "ligand_max_force_atom": top_description(ligand_max),
                "ligand_max_force_kJ_mol_nm": f"{ligand_max['force_magnitude_kJ_mol_nm']:.12g}",
                "ligand_MEMB49_min_distance_A": f"{ligand_memb49_min_distance:.6g}",
                "restraint_atom_count": restraint_count if meta["force_label"] == "CustomExternalForce/restraints" else "",
            }
        )
        for row in memb49_force_rows:
            atom_rows.append(
                {
                    "scope": "MEMB49",
                    "force_group": group,
                    "force_class": meta["force_class"],
                    "force_label": meta["force_label"],
                    **row,
                }
            )
        for row in ligand_force_rows:
            atom_rows.append(
                {
                    "scope": "LIGAND_L002",
                    "force_group": group,
                    "force_class": meta["force_class"],
                    "force_label": meta["force_label"],
                    **row,
                }
            )
        group_results.append({"meta": meta, "energy": energy, "has_nan": has_nan, "memb49_max": memb49_max, "ligand_max": ligand_max})

    dominant = max(summary_rows, key=lambda row: float(row["MEMB49_max_force_kJ_mol_nm"]))
    write_csv(SUMMARY_CSV, summary_rows)
    write_csv(MEMB49_CSV, atom_rows)
    write_csv(TOP50_CSV, top50_rows)
    write_json(
        SUMMARY_JSON,
        {
            "platform": platform.getName(),
            "total_energy_kJ_mol": total_energy,
            "total_force_nan": total_has_nan,
            "total_max_force_atom": top50_rows[0],
            "dominant_MEMB49_force_group": dominant,
            "ligand_MEMB49_min_distance_A": ligand_memb49_min_distance,
            "restraint_atom_count": restraint_count,
            "force_metadata": force_metadata,
        },
    )
    print(f"MEMB49_max_force_with_restraints_kJ_mol_nm: {dominant['MEMB49_max_force_kJ_mol_nm']}")
    print(f"dominant_force_group: {dominant['force_label']} ({dominant['force_class']}, group {dominant['force_group']})")
    print(f"MEMB49_dominant_force_atom: {dominant['MEMB49_max_force_atom']}")
    print(f"total_max_force_atom: {top_description(top50_rows[0])} force={top50_rows[0]['force_magnitude_kJ_mol_nm']:.12g}")
    print(f"total_energy_with_stage4_restraints_kJ_mol: {total_energy:.12g}")
    print(f"force_nan: {total_has_nan}")
    print(f"ligand_MEMB49_min_distance_A: {ligand_memb49_min_distance:.6g}")
    print(f"Summary CSV: {SUMMARY_CSV}")
    print(f"MEMB49/Ligand CSV: {MEMB49_CSV}")
    print(f"Top50 CSV: {TOP50_CSV}")
    del context
    del integrator
    return 0 if math.isfinite(total_energy) and not total_has_nan else 1


if __name__ == "__main__":
    raise SystemExit(main())
