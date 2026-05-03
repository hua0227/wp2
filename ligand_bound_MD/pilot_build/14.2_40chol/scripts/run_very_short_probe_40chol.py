from __future__ import annotations

import math

from openmm import LangevinMiddleIntegrator, Platform
from openmm.app import HBonds, Simulation
from openmm.unit import kelvin, kilojoules_per_mole, picoseconds

from common_40chol import (
    CLEANED_PSF,
    MINIMIZED_PDB,
    PROBE_FINAL_PDB,
    REPORTS_DIR,
    build_system,
    choose_opencl_platform,
    force_rows,
    pbc_restraint_force,
    positions_to_angstrom,
    r2_restraint_map,
    read_pdb_atoms,
    retention_and_clashes,
    write_csv,
    write_pdb,
)


CSV_OUT = REPORTS_DIR / "very_short_probe_40chol.csv"
TIMESTEP_PS = 0.00005
TOTAL_STEPS = 4000
REPORT_INTERVAL = 100
SAFETY_INTERVAL = 10


def estimate_temperature_k(kinetic_energy_kj_mol: float, system) -> float:
    dof = 3 * system.getNumParticles() - system.getNumConstraints()
    for i in range(system.getNumForces()):
        if system.getForce(i).__class__.__name__ == "CMMotionRemover":
            dof -= 3
            break
    return 2 * kinetic_energy_kj_mol / (dof * 0.00831446261815324)


def main() -> int:
    print("40chol very-short 0.2 ps 50 K probe")
    print("OpenCL only. This is not production MD.")
    atoms = read_pdb_atoms(MINIMIZED_PDB)
    psf, pdb, system = build_system(CLEANED_PSF, MINIMIZED_PDB, constraints=HBonds, rigid_water=True)
    restraints = r2_restraint_map(atoms, protein_k=500.0, ligand_k=25.0, lipid_k=5.0)
    system.addForce(pbc_restraint_force(restraints, pdb.positions))
    platform = choose_opencl_platform(Platform)
    integrator = LangevinMiddleIntegrator(50 * kelvin, 10.0 / picoseconds, TIMESTEP_PS * picoseconds)
    integrator.setRandomNumberSeed(142400)
    simulation = Simulation(psf.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    simulation.context.setVelocitiesToTemperature(50 * kelvin, 142400)
    ref_coords = positions_to_angstrom(pdb.positions)

    fields = [
        "step",
        "temperature_K",
        "potential_energy_kJ_mol",
        "kinetic_energy_kJ_mol",
        "max_force_kJ_mol_nm",
        "top_force_atom",
        "severe_lt_1A",
        "close_lt_1p5A",
        "ligand_key_min_distance_A",
        "ligand_near_chains_count",
        "ligand_rmsd_A",
        "status",
    ]
    rows = []
    stopped = False
    stop_reason = ""
    step = 0
    last_positions = pdb.positions
    while step <= TOTAL_STEPS:
        state = simulation.context.getState(getEnergy=True, getForces=True, getPositions=True)
        pe = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        ke = state.getKineticEnergy().value_in_unit(kilojoules_per_mole)
        temp = estimate_temperature_k(ke, system)
        coords = positions_to_angstrom(state.getPositions(asNumpy=False))
        force_table, force_nan = force_rows(atoms, state.getForces(asNumpy=False))
        retention, clashes = retention_and_clashes(atoms, ref_coords, coords)
        statuses = []
        if not all(math.isfinite(value) for value in [pe, ke, temp]) or force_nan:
            statuses.append("nonfinite_energy_or_force")
            stopped = True
        if math.isfinite(temp) and temp > 200:
            statuses.append("temperature_above_200K")
            stopped = True
        if retention["near_chains"] == 0:
            statuses.append("ligand_near_chains_0")
            stopped = True
        if not statuses:
            statuses.append("ok")
        top = force_table[0]
        row = {
            "step": step,
            "temperature_K": f"{temp:.6g}",
            "potential_energy_kJ_mol": f"{pe:.12g}",
            "kinetic_energy_kJ_mol": f"{ke:.12g}",
            "max_force_kJ_mol_nm": f"{top['force_magnitude_kJ_mol_nm']:.12g}",
            "top_force_atom": f"{top['atom_name']} {top['resname']} {top['segid']}:{top['resid']} index={top['index']}",
            "severe_lt_1A": clashes["severe_lt_1A"],
            "close_lt_1p5A": clashes["close_lt_1p5A"],
            "ligand_key_min_distance_A": f"{retention['key_min_distance_A']:.6g}",
            "ligand_near_chains_count": retention["near_chains"],
            "ligand_rmsd_A": f"{retention['ligand_rmsd_A']:.6g}",
            "status": ";".join(statuses),
        }
        if step == 0 or step % REPORT_INTERVAL == 0 or stopped or step == TOTAL_STEPS:
            rows.append(row)
            print(
                f"step {step} temp {temp:.3f} PE {pe:.6g} maxF {top['force_magnitude_kJ_mol_nm']:.6g} "
                f"chains {retention['near_chains']} status {row['status']}"
            )
        last_positions = state.getPositions(asNumpy=False)
        if stopped or step == TOTAL_STEPS:
            stop_reason = row["status"] if stopped else ""
            break
        chunk = min(SAFETY_INTERVAL, TOTAL_STEPS - step)
        simulation.step(chunk)
        step += chunk

    write_csv(CSV_OUT, rows, fields)
    write_pdb(
        MINIMIZED_PDB,
        positions_to_angstrom(last_positions),
        PROBE_FINAL_PDB,
        ["REMARK 40chol very-short PBC-aware probe final coordinates", "REMARK not production MD", f"REMARK stop reason {stop_reason or 'none'}"],
    )
    completed = (not stopped) and step == TOTAL_STEPS
    print("very_short_probe_executed: True")
    print(f"very_short_probe_passed: {completed}")
    print(f"last_step: {step}/{TOTAL_STEPS}")
    print(f"stop_reason: {stop_reason or 'none'}")
    print(f"probe CSV: {CSV_OUT}")
    return 0 if completed else 1


if __name__ == "__main__":
    raise SystemExit(main())
