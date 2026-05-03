from __future__ import annotations

import math

from strategy_common import (
    BRANCH_B_CLEANED_PSF,
    OUTPUTS_DIR,
    REPORTS_DIR,
    STAGE4_PDB,
    build_system,
    choose_opencl_platform,
    coords_have_nan,
    ensure_dirs,
    force_rows,
    is_cholesterol,
    is_heavy,
    is_ligand,
    is_lipid,
    is_memb49,
    is_proa437_438,
    is_protein,
    pbc_restraint_force,
    positions_to_angstrom,
    read_pdb_atoms,
    retention_and_clashes,
    write_csv,
    write_json,
    write_pdb,
    PROTEIN_BACKBONE_NAMES,
)


COMPARISON_CSV = REPORTS_DIR / "pbcaware_minimization_comparison.csv"
SUMMARY_JSON = REPORTS_DIR / "pbcaware_minimization_comparison.json"
R1_PDB = OUTPUTS_DIR / "R1_pbcaware_minimized.pdb"
R2_PDB = OUTPUTS_DIR / "R2_pbcaware_minimized.pdb"
MAX_ITERATIONS = 20000


def restraint_map(atoms, branch: str) -> dict[int, float]:
    m: dict[int, float] = {}
    for i, atom in enumerate(atoms):
        k = 0.0
        if is_protein(atom) and atom.name in PROTEIN_BACKBONE_NAMES:
            k = 500.0
        elif is_ligand(atom) and is_heavy(atom):
            k = 25.0
        elif branch == "R2" and (is_lipid(atom) or is_cholesterol(atom)) and is_heavy(atom):
            k = 5.0
        if k > 0:
            m[i] = k
    return m


def desc(row: dict[str, object]) -> str:
    return f"{row['atom_name']} {row['resname']} {row['segid']}:{row['resid']} index={row['index']}"


def evaluate(atoms, coords, ref_coords, state, energy) -> dict[str, object]:
    rows, force_nan = force_rows(atoms, state.getForces(asNumpy=False))
    top = rows[0]
    memb49_rows = [r for r in rows if is_memb49(atoms[int(r["index"])])]
    pro_rows = [r for r in rows if is_proa437_438(atoms[int(r["index"])])]
    memb49 = max(memb49_rows, key=lambda r: r["force_magnitude_kJ_mol_nm"])
    pro = max(pro_rows, key=lambda r: r["force_magnitude_kJ_mol_nm"])
    retention, clashes = retention_and_clashes(atoms, ref_coords, coords)
    return {
        "final_energy_kJ_mol": energy,
        "energy_finite": math.isfinite(energy),
        "force_nan": force_nan,
        "coordinate_nan": coords_have_nan(coords),
        "max_force_kJ_mol_nm": top["force_magnitude_kJ_mol_nm"],
        "top_force_atom": desc(top),
        "MEMB49_max_force_kJ_mol_nm": memb49["force_magnitude_kJ_mol_nm"],
        "MEMB49_max_force_atom": desc(memb49),
        "PROA437_438_max_force_kJ_mol_nm": pro["force_magnitude_kJ_mol_nm"],
        "PROA437_438_max_force_atom": desc(pro),
        "ligand_key_min_distance_A": retention["key_min_distance_A"],
        "ligand_near_chains_count": retention["near_chains"],
        "ligand_rmsd_vs_BranchB_start_A": retention["ligand_rmsd_A"],
        "severe_lt_1A": clashes["severe_lt_1A"],
        "close_lt_1p5A": clashes["close_lt_1p5A"],
    }


def run_branch(branch: str, output_pdb):
    import openmm
    from openmm import LocalEnergyMinimizer, Platform, VerletIntegrator
    from openmm.app import HBonds
    from openmm.unit import kilojoules_per_mole, nanometer, picoseconds

    atoms = read_pdb_atoms(STAGE4_PDB)
    _psf, pdb, system = build_system(BRANCH_B_CLEANED_PSF, STAGE4_PDB, constraints=HBonds, rigid_water=True)
    rmap = restraint_map(atoms, branch)
    pbc_force = pbc_restraint_force(rmap, pdb.positions)
    pbc_force.setForceGroup(system.getNumForces())
    system.addForce(pbc_force)
    platform = choose_opencl_platform(Platform)
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = openmm.Context(system, integrator, platform)
    context.setPositions(pdb.positions)
    ref_coords = positions_to_angstrom(pdb.positions)
    initial = context.getState(getEnergy=True, getForces=True)
    initial_energy = initial.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    print(f"{branch} initial restrained energy: {initial_energy:.12g} kJ/mol; restraints={len(rmap)}")
    LocalEnergyMinimizer.minimize(context, tolerance=1.0 * kilojoules_per_mole / nanometer, maxIterations=MAX_ITERATIONS)
    state = context.getState(getEnergy=True, getForces=True, getPositions=True)
    energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    coords = positions_to_angstrom(state.getPositions(asNumpy=False))
    result = evaluate(atoms, coords, ref_coords, state, energy)
    result.update({
        "branch": branch,
        "restraint_atom_count": len(rmap),
        "initial_energy_kJ_mol": initial_energy,
        "output_pdb": str(output_pdb),
    })
    write_pdb(
        STAGE4_PDB,
        coords,
        output_pdb,
        [
            "REMARK PBC-aware restraint minimization",
            f"REMARK branch {branch}",
            "REMARK no MD was run",
            f"REMARK final energy {energy:.12g} kJ/mol",
        ],
    )
    print(
        f"{branch} final energy {energy:.12g}; max force {result['max_force_kJ_mol_nm']:.12g} at {result['top_force_atom']}; "
        f"MEMB49 {result['MEMB49_max_force_kJ_mol_nm']:.12g}; PROA437/438 {result['PROA437_438_max_force_kJ_mol_nm']:.12g}; "
        f"severe/close {result['severe_lt_1A']}/{result['close_lt_1p5A']}; chains {result['ligand_near_chains_count']}"
    )
    del context
    del integrator
    return result


def main() -> int:
    ensure_dirs()
    results = [run_branch("R1", R1_PDB), run_branch("R2", R2_PDB)]
    rows = []
    for r in results:
        rows.append({k: (f"{v:.12g}" if isinstance(v, float) else int(v) if isinstance(v, bool) else v) for k, v in r.items()})
    write_csv(COMPARISON_CSV, rows)
    write_json(SUMMARY_JSON, {"results": results, "max_iterations": MAX_ITERATIONS})
    print(f"Comparison CSV: {COMPARISON_CSV}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
