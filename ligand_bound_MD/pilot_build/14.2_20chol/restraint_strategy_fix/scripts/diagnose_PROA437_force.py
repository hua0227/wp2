from __future__ import annotations

import math

from strategy_common import (
    BRANCH_B_CLEANED_PSF,
    REPORTS_DIR,
    STAGE4_PDB,
    atom_category,
    build_system,
    choose_opencl_platform,
    distance_a,
    ensure_dirs,
    force_rows,
    fmt,
    is_ligand,
    is_proa437_438,
    read_gro_box_lengths_a,
    read_pdb_atoms,
    write_csv,
)


FORCE_CSV = REPORTS_DIR / "PROA437_force_decomposition.csv"
CONTACTS_CSV = REPORTS_DIR / "PROA437_contacts.csv"
MD_OUT = REPORTS_DIR / "PROA437_force_diagnostics.md"


def force_label(force: object) -> str:
    name = force.__class__.__name__
    if name in {"NonbondedForce", "CustomNonbondedForce"}:
        return "nonbonded"
    if name == "HarmonicBondForce":
        return "bonded_internal"
    if name == "HarmonicAngleForce":
        return "angle"
    if name in {"PeriodicTorsionForce", "CustomTorsionForce", "CMAPTorsionForce"}:
        return "torsion"
    return "other"


def desc(row: dict[str, object]) -> str:
    return f"{row['atom_name']} {row['resname']} {row['segid']}:{row['resid']} index={row['index']}"


def main() -> int:
    import openmm
    from openmm import Platform, VerletIntegrator
    from openmm.app import HBonds
    from openmm.unit import kilojoules_per_mole, picoseconds

    ensure_dirs()
    atoms = read_pdb_atoms(STAGE4_PDB)
    _psf, pdb, system = build_system(BRANCH_B_CLEANED_PSF, STAGE4_PDB, constraints=HBonds, rigid_water=True)
    for i in range(system.getNumForces()):
        system.getForce(i).setForceGroup(i)
    platform = choose_opencl_platform(Platform)
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = openmm.Context(system, integrator, platform)
    context.setPositions(pdb.positions)

    pro_indices = [i for i, a in enumerate(atoms) if is_proa437_438(a)]
    ligand_indices = [i for i, a in enumerate(atoms) if is_ligand(a)]
    target_o_index = next((i for i, a in enumerate(atoms) if a.segid == "PROA" and a.resid == "437" and a.name == "O"), None)
    all_rows = []
    group_rows = []
    for group in range(system.getNumForces()):
        force = system.getForce(group)
        state = context.getState(getEnergy=True, getForces=True, groups=1 << group)
        energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        rows, has_nan = force_rows(atoms, state.getForces(asNumpy=False))
        pro_rows = [r for r in rows if int(r["index"]) in pro_indices]
        lig_rows = [r for r in rows if int(r["index"]) in ligand_indices]
        o_row = next((r for r in pro_rows if int(r["index"]) == target_o_index), None)
        pro_max = max(pro_rows, key=lambda r: r["force_magnitude_kJ_mol_nm"])
        lig_max = max(lig_rows, key=lambda r: r["force_magnitude_kJ_mol_nm"])
        group_rows.append({
            "force_group": group,
            "force_class": force.__class__.__name__,
            "force_source_classification": force_label(force),
            "energy_kJ_mol": f"{energy:.12g}",
            "force_nan": int(has_nan),
            "PROA437_438_max_force_atom": desc(pro_max),
            "PROA437_438_max_force_kJ_mol_nm": f"{pro_max['force_magnitude_kJ_mol_nm']:.12g}",
            "PROA437_O_force_kJ_mol_nm": f"{o_row['force_magnitude_kJ_mol_nm']:.12g}" if o_row else "nan",
            "ligand_max_force_atom": desc(lig_max),
            "ligand_max_force_kJ_mol_nm": f"{lig_max['force_magnitude_kJ_mol_nm']:.12g}",
        })
        for r in pro_rows:
            all_rows.append({"scope": "PROA437_438", "force_group": group, "force_class": force.__class__.__name__, "force_source_classification": force_label(force), **r})
        for r in lig_rows:
            all_rows.append({"scope": "LIGAND_L002", "force_group": group, "force_class": force.__class__.__name__, "force_source_classification": force_label(force), **r})

    total_state = context.getState(getEnergy=True, getForces=True)
    total_energy = total_state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    total_rows, total_nan = force_rows(atoms, total_state.getForces(asNumpy=False))
    total_pro = [r for r in total_rows if int(r["index"]) in pro_indices]
    total_o = next(r for r in total_pro if int(r["index"]) == target_o_index)
    dominant_o = max(group_rows, key=lambda r: float(r["PROA437_O_force_kJ_mol_nm"]) if r["PROA437_O_force_kJ_mol_nm"] != "nan" else -1)

    coords = [(a.x, a.y, a.z) for a in atoms]
    box = read_gro_box_lengths_a()
    contact_rows = []
    for i in pro_indices:
        for j, other in enumerate(atoms):
            if j in pro_indices:
                continue
            d = distance_a(coords[i], coords[j], box)
            if d <= 5.0:
                contact_rows.append({
                    "distance_A": f"{d:.6g}",
                    "target_index": i,
                    "target_atom": atoms[i].name,
                    "target_resid": atoms[i].resid,
                    "other_index": j,
                    "other_atom": other.name,
                    "other_resname": other.resname,
                    "other_segid": other.segid,
                    "other_resid": other.resid,
                    "other_category": atom_category(other),
                    "other_is_ligand": int(is_ligand(other)),
                })
    contact_rows.sort(key=lambda r: float(r["distance_A"]))
    force_fieldnames = [
        "scope","force_group","force_class","force_source_classification","energy_kJ_mol","force_nan",
        "PROA437_438_max_force_atom","PROA437_438_max_force_kJ_mol_nm","PROA437_O_force_kJ_mol_nm",
        "ligand_max_force_atom","ligand_max_force_kJ_mol_nm",
        "index","serial","atom_name","resname","segid","resid","element","category",
        "force_x_kJ_mol_nm","force_y_kJ_mol_nm","force_z_kJ_mol_nm","force_magnitude_kJ_mol_nm",
        "is_MEMB49","is_PROA437_438","is_ligand",
    ]
    write_csv(FORCE_CSV, group_rows + all_rows, force_fieldnames)
    write_csv(CONTACTS_CSV, contact_rows)

    lines = [
        "# PROA437/438 Force Diagnostics",
        "",
        "No restraints, no minimization, and no MD were run.",
        "",
        f"- Total potential energy: {fmt(total_energy, 12)} kJ/mol",
        f"- Force NaN: {'YES' if total_nan else 'NO'}",
        f"- Total PROA:437 O force: {fmt(total_o['force_magnitude_kJ_mol_nm'], 12)} kJ/mol/nm",
        f"- PROA:437 O dominant force group: {dominant_o['force_source_classification']} ({dominant_o['force_class']}, group {dominant_o['force_group']})",
        f"- Dominant group PROA:437 O force: {dominant_o['PROA437_O_force_kJ_mol_nm']} kJ/mol/nm",
        f"- Ligand contacts within 5 A of PROA:437/438: {sum(1 for r in contact_rows if r['other_is_ligand'] == 1)}",
        "",
        "## Force Group Summary",
        "",
        "| Group | Force class | Source | Energy kJ/mol | PROA437 O force | PROA437/438 max atom | PROA437/438 max force |",
        "|---:|---|---|---:|---:|---|---:|",
    ]
    for r in group_rows:
        lines.append(f"| {r['force_group']} | {r['force_class']} | {r['force_source_classification']} | {r['energy_kJ_mol']} | {r['PROA437_O_force_kJ_mol_nm']} | {r['PROA437_438_max_force_atom']} | {r['PROA437_438_max_force_kJ_mol_nm']} |")
    lines.extend(["", "## Closest Contacts", "", "| Distance A | Target | Other | Category | Ligand |", "|---:|---|---|---|---:|"])
    for r in contact_rows[:30]:
        lines.append(f"| {r['distance_A']} | PROA:{r['target_resid']} {r['target_atom']} | {r['other_atom']} {r['other_resname']} {r['other_segid']}:{r['other_resid']} | {r['other_category']} | {r['other_is_ligand']} |")
    MD_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"PROA437_dominant_force_group: {dominant_o['force_source_classification']} ({dominant_o['force_class']}, group {dominant_o['force_group']})")
    print(f"PROA437_O_total_force_kJ_mol_nm: {total_o['force_magnitude_kJ_mol_nm']:.12g}")
    print(f"PROA437_O_dominant_group_force_kJ_mol_nm: {dominant_o['PROA437_O_force_kJ_mol_nm']}")
    print(f"ligand_contacts_within_5A: {sum(1 for r in contact_rows if r['other_is_ligand'] == 1)}")
    print(f"Force CSV: {FORCE_CSV}")
    print(f"Contacts CSV: {CONTACTS_CSV}")
    print(f"Markdown: {MD_OUT}")
    del context
    del integrator
    return 0 if not total_nan and math.isfinite(total_energy) else 1


if __name__ == "__main__":
    raise SystemExit(main())
