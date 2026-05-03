from __future__ import annotations

import csv
import math

from strategy_common import (
    BRANCH_B_CLEANED_PSF,
    LOGS_DIR,
    OUTPUTS_DIR,
    REPORTS_DIR,
    build_system,
    choose_opencl_platform,
    force_rows,
    is_cholesterol,
    is_heavy,
    is_ligand,
    is_lipid,
    is_protein,
    pbc_restraint_force,
    positions_to_angstrom,
    read_pdb_atoms,
    retention_and_clashes,
    write_csv,
    PROTEIN_BACKBONE_NAMES,
)


INPUT_PDB = OUTPUTS_DIR / "R2_pbcaware_minimized.pdb"
CSV_OUT = REPORTS_DIR / "very_short_probe_pbcaware.csv"
FINAL_PDB = OUTPUTS_DIR / "R2_very_short_probe_final.pdb"
TIMESTEP_PS = 0.00005
TOTAL_STEPS = 4000  # 0.2 ps
REPORT_INTERVAL = 100
SAFETY_INTERVAL = 10


def restraint_map(atoms):
    m = {}
    for i, atom in enumerate(atoms):
        k = 0.0
        if is_protein(atom) and atom.name in PROTEIN_BACKBONE_NAMES:
            k = 1000.0
        elif is_ligand(atom) and is_heavy(atom):
            k = 100.0
        elif (is_lipid(atom) or is_cholesterol(atom)) and is_heavy(atom):
            k = 10.0
        if k > 0:
            m[i] = k
    return m


def estimate_temp(ke, system):
    dof = 3 * system.getNumParticles() - system.getNumConstraints()
    for i in range(system.getNumForces()):
        if system.getForce(i).__class__.__name__ == "CMMotionRemover":
            dof -= 3
            break
    return 2 * ke / (dof * 0.00831446261815324)


def main() -> int:
    import openmm
    from openmm import LangevinMiddleIntegrator, Platform
    from openmm.app import HBonds, Simulation
    from openmm.unit import kelvin, kilojoules_per_mole, picoseconds

    atoms = read_pdb_atoms(INPUT_PDB)
    psf, pdb, system = build_system(BRANCH_B_CLEANED_PSF, INPUT_PDB, constraints=HBonds, rigid_water=True)
    rmap = restraint_map(atoms)
    system.addForce(pbc_restraint_force(rmap, pdb.positions))
    platform = choose_opencl_platform(Platform)
    integrator = LangevinMiddleIntegrator(50 * kelvin, 10.0 / picoseconds, TIMESTEP_PS * picoseconds)
    integrator.setRandomNumberSeed(142050)
    sim = Simulation(psf.topology, system, integrator, platform)
    sim.context.setPositions(pdb.positions)
    sim.context.setVelocitiesToTemperature(50 * kelvin, 142050)
    ref = positions_to_angstrom(pdb.positions)
    fields = ["step","temperature_K","potential_energy_kJ_mol","kinetic_energy_kJ_mol","max_force_kJ_mol_nm","top_force_atom","severe_lt_1A","close_lt_1p5A","ligand_key_min_distance_A","ligand_near_chains_count","ligand_rmsd_A","status"]
    rows = []
    stopped = False
    stop_reason = ""
    last_positions = pdb.positions
    step = 0
    while step <= TOTAL_STEPS:
        state = sim.context.getState(getEnergy=True, getForces=True, getPositions=True)
        pe = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        ke = state.getKineticEnergy().value_in_unit(kilojoules_per_mole)
        temp = estimate_temp(ke, system)
        coords = positions_to_angstrom(state.getPositions(asNumpy=False))
        frows, fnan = force_rows(atoms, state.getForces(asNumpy=False))
        retention, clashes = retention_and_clashes(atoms, ref, coords)
        statuses = []
        if not all(math.isfinite(x) for x in [pe, ke, temp]) or fnan:
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
        row = {
            "step": step, "temperature_K": f"{temp:.6g}", "potential_energy_kJ_mol": f"{pe:.12g}", "kinetic_energy_kJ_mol": f"{ke:.12g}",
            "max_force_kJ_mol_nm": f"{frows[0]['force_magnitude_kJ_mol_nm']:.12g}",
            "top_force_atom": f"{frows[0]['atom_name']} {frows[0]['resname']} {frows[0]['segid']}:{frows[0]['resid']} index={frows[0]['index']}",
            "severe_lt_1A": clashes["severe_lt_1A"], "close_lt_1p5A": clashes["close_lt_1p5A"],
            "ligand_key_min_distance_A": f"{retention['key_min_distance_A']:.6g}", "ligand_near_chains_count": retention["near_chains"],
            "ligand_rmsd_A": f"{retention['ligand_rmsd_A']:.6g}", "status": ";".join(statuses)
        }
        if step == 0 or step % REPORT_INTERVAL == 0 or stopped or step == TOTAL_STEPS:
            rows.append(row)
            print(f"step {step} temp {temp:.3f} PE {pe:.6g} maxF {frows[0]['force_magnitude_kJ_mol_nm']:.6g} chains {retention['near_chains']} status {row['status']}")
        last_positions = state.getPositions(asNumpy=False)
        if stopped or step == TOTAL_STEPS:
            stop_reason = row["status"] if stopped else ""
            break
        chunk = min(SAFETY_INTERVAL, TOTAL_STEPS - step)
        sim.step(chunk)
        step += chunk
    write_csv(CSV_OUT, rows, fields)
    from strategy_common import write_pdb
    write_pdb(INPUT_PDB, positions_to_angstrom(last_positions), FINAL_PDB, ["REMARK very-short PBC-aware probe final safe coordinates", "REMARK not production MD", f"REMARK stop reason {stop_reason or 'none'}"])
    completed = (not stopped) and step == TOTAL_STEPS
    print(f"very_short_probe_executed: True")
    print(f"very_short_probe_passed: {completed}")
    print(f"last_step: {step}/{TOTAL_STEPS}")
    print(f"stop_reason: {stop_reason or 'none'}")
    return 0 if completed else 1


if __name__ == "__main__":
    raise SystemExit(main())
