from __future__ import annotations

import math
from pathlib import Path
from typing import Any

from branchB_common import (
    BRANCH_B_MINIMIZED_PDB,
    CLEANED_PSF,
    OUTPUTS_DIR,
    REPORTS_DIR,
    add_position_restraints,
    build_charmm_system,
    choose_opencl_platform,
    clash_metrics,
    ensure_dirs,
    fmt,
    force_rows,
    image_ligand_coords,
    ligand_image_translation_a,
    positions_to_angstrom,
    read_gro_box_lengths_a,
    read_pdb_atoms,
    restraint_map_for_stage,
    retention_metrics,
    write_csv,
    write_json,
    write_pdb_with_coords,
)


SUMMARY_CSV = REPORTS_DIR / "staged_minimization_summary.csv"
SUMMARY_JSON = REPORTS_DIR / "staged_minimization_summary.json"
STAGE_OUTPUTS = {
    1: OUTPUTS_DIR / "branchB_stage1_minimized.pdb",
    2: OUTPUTS_DIR / "branchB_stage2_minimized.pdb",
    3: OUTPUTS_DIR / "branchB_stage3_minimized.pdb",
    4: OUTPUTS_DIR / "branchB_stage4_minimized.pdb",
}

MAX_ITERATIONS = 10000
TOLERANCE_KJ_MOL_NM = 1.0


def coords_to_positions(coords_a: list[tuple[float, float, float]]) -> Any:
    from openmm import Vec3
    from openmm.unit import nanometer

    return [Vec3(x / 10.0, y / 10.0, z / 10.0) for x, y, z in coords_a] * nanometer


def top_force_description(row: dict[str, Any]) -> str:
    return (
        f"{row['atom_name']} {row['resname']} {row['segid']}:{row['resid']} "
        f"index={row['index']}"
    )


def run() -> dict[str, Any]:
    import openmm
    from openmm import LocalEnergyMinimizer, Platform, VerletIntegrator
    from openmm.app import HBonds
    from openmm.unit import kilojoules_per_mole, nanometer, picoseconds

    ensure_dirs()
    atoms = read_pdb_atoms(BRANCH_B_MINIMIZED_PDB)
    box_lengths_a = read_gro_box_lengths_a()
    _psf_ref, pdb_ref, _params_ref, physical_system = build_charmm_system(
        CLEANED_PSF,
        BRANCH_B_MINIMIZED_PDB,
        constraints=HBonds,
        rigid_water=True,
    )
    platform = choose_opencl_platform(Platform)
    physical_integrator = VerletIntegrator(0.001 * picoseconds)
    physical_context = openmm.Context(physical_system, physical_integrator, platform)

    reference_positions = pdb_ref.positions
    reference_coords_a = positions_to_angstrom(reference_positions)
    current_positions = pdb_ref.positions
    rows: list[dict[str, Any]] = []
    all_results: dict[str, Any] = {
        "platform": platform.getName(),
        "max_iterations_per_stage": MAX_ITERATIONS,
        "tolerance_kJ_mol_nm": TOLERANCE_KJ_MOL_NM,
        "stages": [],
    }

    print("Branch B staged deeper minimization")
    print("No velocities are set. No dynamics is run.")
    print(f"OpenMM platform: {platform.getName()}")
    print(f"Input PSF: {CLEANED_PSF}")
    print(f"Input PDB: {BRANCH_B_MINIMIZED_PDB}")
    print(f"Max iterations per stage: {MAX_ITERATIONS}")

    for stage in [1, 2, 3, 4]:
        print("")
        print(f"Stage {stage} starting")
        _psf, _pdb, _params, system = build_charmm_system(
            CLEANED_PSF,
            BRANCH_B_MINIMIZED_PDB,
            constraints=HBonds,
            rigid_water=True,
        )
        restraint_map = restraint_map_for_stage(atoms, stage)
        restraint_count = add_position_restraints(system, atoms, reference_positions, restraint_map)
        integrator = VerletIntegrator(0.001 * picoseconds)
        context = openmm.Context(system, integrator, platform)
        context.setPositions(current_positions)

        initial_state = context.getState(getEnergy=True, getForces=True)
        initial_restrained_energy = initial_state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        initial_restrained_forces, initial_restrained_nan = force_rows(atoms, initial_state.getForces(asNumpy=False))
        print(
            f"Stage {stage} initial restrained energy {initial_restrained_energy:.12g} kJ/mol; "
            f"max force {initial_restrained_forces[0]['force_magnitude_kJ_mol_nm']:.12g} "
            f"at {top_force_description(initial_restrained_forces[0])}; restraint atoms {restraint_count}"
        )

        LocalEnergyMinimizer.minimize(
            context,
            tolerance=TOLERANCE_KJ_MOL_NM * kilojoules_per_mole / nanometer,
            maxIterations=MAX_ITERATIONS,
        )
        final_state = context.getState(getEnergy=True, getForces=True, getPositions=True)
        final_restrained_energy = final_state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        final_restrained_forces, final_restrained_nan = force_rows(atoms, final_state.getForces(asNumpy=False))
        raw_final_positions = final_state.getPositions(asNumpy=False)
        raw_final_coords_a = positions_to_angstrom(raw_final_positions)
        image_translation = ligand_image_translation_a(atoms, reference_coords_a, raw_final_coords_a, box_lengths_a)
        imaged_coords_a = image_ligand_coords(atoms, raw_final_coords_a, image_translation)
        current_positions = coords_to_positions(imaged_coords_a)

        physical_context.setPositions(current_positions)
        physical_state = physical_context.getState(getEnergy=True, getForces=True)
        physical_energy = physical_state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        physical_forces, physical_force_nan = force_rows(atoms, physical_state.getForces(asNumpy=False))
        physical_max = physical_forces[0]
        clashes = clash_metrics(atoms, imaged_coords_a, box_lengths_a)
        retention = retention_metrics(atoms, reference_coords_a, imaged_coords_a, box_lengths_a)
        energy_finite = math.isfinite(final_restrained_energy) and math.isfinite(physical_energy)
        force_nan = final_restrained_nan or physical_force_nan
        output_pdb = STAGE_OUTPUTS[stage]
        write_pdb_with_coords(
            BRANCH_B_MINIMIZED_PDB,
            imaged_coords_a,
            output_pdb,
            [
                "REMARK Branch B staged minimization output",
                f"REMARK stage {stage}",
                "REMARK no dynamics was run",
                f"REMARK physical energy {physical_energy:.12g} kJ/mol",
                f"REMARK ligand image translation A {image_translation[0]:.4f} {image_translation[1]:.4f} {image_translation[2]:.4f}",
            ],
        )
        row = {
            "stage": stage,
            "restraint_atom_count": restraint_count,
            "initial_restrained_energy_kJ_mol": f"{initial_restrained_energy:.12g}",
            "final_restrained_energy_kJ_mol": f"{final_restrained_energy:.12g}",
            "physical_energy_kJ_mol": f"{physical_energy:.12g}",
            "energy_finite": int(energy_finite),
            "restrained_max_force_atom": top_force_description(final_restrained_forces[0]),
            "restrained_max_force_kJ_mol_nm": f"{final_restrained_forces[0]['force_magnitude_kJ_mol_nm']:.12g}",
            "physical_max_force_atom": top_force_description(physical_max),
            "physical_max_force_kJ_mol_nm": f"{physical_max['force_magnitude_kJ_mol_nm']:.12g}",
            "force_nan": int(force_nan),
            "severe_lt_1A": clashes["severe_lt_1A"],
            "close_lt_1p5A": clashes["close_lt_1p5A"],
            "h9_h11_popc_severe_lt_1A": clashes["h9_h11_popc_severe_lt_1A"],
            "h9_h11_popc_close_lt_1p5A": clashes["h9_h11_popc_close_lt_1p5A"],
            "ligand_key_min_distance_A": f"{retention['key_min_distance_a']:.6g}",
            "ligand_near_chains_count": retention["near_chain_count"],
            "ligand_rmsd_vs_branchB_start_A": f"{retention['ligand_rmsd_a']:.6g}",
            "ligand_retained_near_key_residues": int(retention["key_retained"]),
            "ligand_retained_near_two_chains": int(retention["near_two_chains"]),
            "ligand_image_translation_A": f"{image_translation[0]:.4f},{image_translation[1]:.4f},{image_translation[2]:.4f}",
            "output_pdb": str(output_pdb),
        }
        rows.append(row)
        all_results["stages"].append(
            {
                **row,
                "retention": retention,
                "clashes": clashes,
                "physical_max_force": physical_max,
                "restrained_max_force": final_restrained_forces[0],
            }
        )
        print(
            f"Stage {stage} final physical energy {physical_energy:.12g} kJ/mol; "
            f"physical max force {physical_max['force_magnitude_kJ_mol_nm']:.12g} at {top_force_description(physical_max)}; "
            f"severe/close {clashes['severe_lt_1A']}/{clashes['close_lt_1p5A']}; "
            f"key_min {retention['key_min_distance_a']:.4f} A; chains {retention['near_chain_count']}; "
            f"ligand RMSD {retention['ligand_rmsd_a']:.4f} A"
        )
        del context
        del integrator

    fieldnames = list(rows[0].keys())
    write_csv(SUMMARY_CSV, rows, fieldnames)
    write_json(SUMMARY_JSON, all_results)
    print("")
    print(f"Summary CSV: {SUMMARY_CSV}")
    print(f"Summary JSON: {SUMMARY_JSON}")
    del physical_context
    del physical_integrator
    return all_results


def main() -> int:
    result = run()
    stage4 = result["stages"][-1]
    ok = (
        int(stage4["energy_finite"]) == 1
        and int(stage4["force_nan"]) == 0
        and int(stage4["severe_lt_1A"]) == 0
        and int(stage4["ligand_retained_near_key_residues"]) == 1
        and int(stage4["ligand_retained_near_two_chains"]) == 1
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
