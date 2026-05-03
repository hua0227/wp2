from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

from run_deeper_restrained_minimization import (
    GRO_BOX,
    PARAMETER_FILES,
    RETENTION_CUTOFF_A,
    add_position_restraints,
    choose_opencl_platform,
    clash_metrics,
    image_ligand_coords,
    json_safe,
    ligand_image_translation_a,
    positions_to_angstrom,
    read_gro_box_lengths_a,
    read_pdb_atoms,
    retention_metrics,
    set_psf_box_from_gro,
    write_pdb_with_coords,
)


PILOT_ROOT = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol")
CLEANUP_ROOT = PILOT_ROOT / "targeted_clash_cleanup"
BRANCH_B_DIR = CLEANUP_ROOT / "outputs" / "branch_B_remove_POPC94"
REPORTS_DIR = CLEANUP_ROOT / "reports"

INPUT_PSF = BRANCH_B_DIR / "branch_B_cleaned.psf"
INPUT_PDB = BRANCH_B_DIR / "branch_B_minimized.pdb"
FINAL_PDB = BRANCH_B_DIR / "branch_B_1ps50K_final.pdb"
SUMMARY_JSON = BRANCH_B_DIR / "branch_B_1ps50K_summary.json"
RETENTION_CSV = REPORTS_DIR / "branch_B_1ps50K_retention.csv"

TIMESTEP_PS = 0.0001  # 0.1 fs
TOTAL_STEPS = 10000  # 1 ps
TEMPERATURE_K = 50.0
TEMPERATURE_ABORT_K = 200.0
SAFETY_INTERVAL_STEPS = 10
REPORT_INTERVAL_STEPS = 100

FIELDS = [
    "step",
    "temperature_K",
    "potential_energy_kJ_mol",
    "kinetic_energy_kJ_mol",
    "total_energy_kJ_mol",
    "ligand_key_min_distance_A",
    "ligand_near_chains_count",
    "ligand_rmsd_A",
    "severe_lt_1A",
    "close_lt_1p5A",
    "h9_h11_popc_severe_lt_1A",
    "h9_h11_popc_close_lt_1p5A",
    "status",
]


def estimate_temperature_k(kinetic_energy_kj_mol: float, system: Any) -> float:
    dof = 3 * system.getNumParticles() - system.getNumConstraints()
    for force_index in range(system.getNumForces()):
        if system.getForce(force_index).__class__.__name__ == "CMMotionRemover":
            dof -= 3
            break
    if dof <= 0:
        return float("nan")
    gas_constant_kj_mol_k = 0.00831446261815324
    return 2.0 * kinetic_energy_kj_mol / (dof * gas_constant_kj_mol_k)


def finite_values(*values: float) -> bool:
    return all(math.isfinite(value) for value in values)


def sample_row(
    step: int,
    potential: float,
    kinetic: float,
    total: float,
    temperature: float,
    retention: dict[str, Any],
    clashes: dict[str, Any],
    status: str,
) -> dict[str, Any]:
    return {
        "step": step,
        "temperature_K": f"{temperature:.6g}" if math.isfinite(temperature) else "nan",
        "potential_energy_kJ_mol": f"{potential:.12g}" if math.isfinite(potential) else "nan",
        "kinetic_energy_kJ_mol": f"{kinetic:.12g}" if math.isfinite(kinetic) else "nan",
        "total_energy_kJ_mol": f"{total:.12g}" if math.isfinite(total) else "nan",
        "ligand_key_min_distance_A": f"{retention['key_min_distance_a']:.6g}" if math.isfinite(retention["key_min_distance_a"]) else "nan",
        "ligand_near_chains_count": retention["near_chain_count"],
        "ligand_rmsd_A": f"{retention['ligand_rmsd_a']:.6g}" if math.isfinite(retention["ligand_rmsd_a"]) else "nan",
        "severe_lt_1A": clashes["severe_lt_1A"],
        "close_lt_1p5A": clashes["close_lt_1p5A"],
        "h9_h11_popc_severe_lt_1A": clashes["h9_h11_popc_severe_lt_1A"],
        "h9_h11_popc_close_lt_1p5A": clashes["h9_h11_popc_close_lt_1p5A"],
        "status": status,
    }


def run_check() -> dict[str, Any]:
    import openmm
    from openmm import LangevinMiddleIntegrator, Platform
    from openmm.app import CharmmParameterSet, CharmmPsfFile, HBonds, PDBFile, PME, Simulation
    from openmm.unit import kelvin, kilojoules_per_mole, nanometer, picoseconds

    print("Ultra-conservative 1 ps stability check for Branch B only")
    print("This is not production MD.")
    print(f"Input PSF: {INPUT_PSF}")
    print(f"Input PDB: {INPUT_PDB}")
    print("OpenMM platform required: OpenCL")
    print(f"Timestep: {TIMESTEP_PS * 1000.0:.3f} fs")
    print(f"Temperature: {TEMPERATURE_K:.1f} K")
    print(f"Temperature safety stop: > {TEMPERATURE_ABORT_K:.1f} K")
    print("No DCD trajectory will be written.")

    atoms = read_pdb_atoms(INPUT_PDB)
    psf = CharmmPsfFile(str(INPUT_PSF))
    pdb = PDBFile(str(INPUT_PDB))
    if len(atoms) != len(pdb.positions):
        raise ValueError(f"Atom count mismatch: parsed {len(atoms)} atoms, OpenMM positions {len(pdb.positions)}")
    set_psf_box_from_gro(psf)
    params = CharmmParameterSet(*[str(path) for path in PARAMETER_FILES])
    system = psf.createSystem(
        params,
        nonbondedMethod=PME,
        nonbondedCutoff=1.2 * nanometer,
        switchDistance=1.0 * nanometer,
        constraints=HBonds,
        rigidWater=True,
    )
    restraint_counts = add_position_restraints(system, atoms, pdb.positions)
    for key, value in restraint_counts.items():
        print(f"{key}: {value}")

    print("Available OpenMM platforms:")
    for index in range(Platform.getNumPlatforms()):
        print(f"- {Platform.getPlatform(index).getName()}")
    platform = choose_opencl_platform(Platform)
    print(f"Selected OpenMM platform: {platform.getName()}")
    print(f"OpenMM version: {openmm.version.version}")

    integrator = LangevinMiddleIntegrator(TEMPERATURE_K * kelvin, 10.0 / picoseconds, TIMESTEP_PS * picoseconds)
    integrator.setRandomNumberSeed(140250)
    simulation = Simulation(psf.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    simulation.context.applyConstraints(1e-8)
    simulation.context.setVelocitiesToTemperature(TEMPERATURE_K * kelvin, 140250)

    reference_coords_a = positions_to_angstrom(pdb.positions)
    box_lengths_a = read_gro_box_lengths_a(GRO_BOX)
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    BRANCH_B_DIR.mkdir(parents=True, exist_ok=True)

    last_finite_positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=False)
    last_finite_step = 0
    last_row: dict[str, Any] | None = None
    max_temperature = float("-inf")
    stopped_early = False
    stop_reason = ""
    nonfinite_detected = False

    with RETENTION_CSV.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        step = 0
        while step <= TOTAL_STEPS:
            state = simulation.context.getState(getEnergy=True, getPositions=True)
            potential = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
            kinetic = state.getKineticEnergy().value_in_unit(kilojoules_per_mole)
            total = potential + kinetic
            temperature = estimate_temperature_k(kinetic, system)
            if math.isfinite(temperature):
                max_temperature = max(max_temperature, temperature)
            positions = state.getPositions(asNumpy=False)
            coords_a = positions_to_angstrom(positions)
            image_translation = ligand_image_translation_a(atoms, reference_coords_a, coords_a, box_lengths_a)
            imaged_coords_a = image_ligand_coords(atoms, coords_a, image_translation)
            retention = retention_metrics(atoms, reference_coords_a, imaged_coords_a, box_lengths_a)
            clashes = clash_metrics(atoms, imaged_coords_a, box_lengths_a)
            status_parts: list[str] = []
            if not finite_values(potential, kinetic, total, temperature):
                status_parts.append("nonfinite_energy_or_temperature")
                nonfinite_detected = True
                stopped_early = True
            if math.isfinite(temperature) and temperature > TEMPERATURE_ABORT_K:
                status_parts.append("temperature_above_200K")
                stopped_early = True
            if clashes["severe_lt_1A"] > 0:
                status_parts.append("severe_clash_present")
            if retention["near_chain_count"] < 2:
                status_parts.append("ligand_near_chains_lt_2")
            if not status_parts:
                status_parts.append("ok")
            status = ";".join(status_parts)
            row = sample_row(step, potential, kinetic, total, temperature, retention, clashes, status)
            if step == 0 or step % REPORT_INTERVAL_STEPS == 0 or stopped_early or step == TOTAL_STEPS:
                writer.writerow(row)
                f.flush()
                print(
                    "step {step} PE {pe:.6g} KE {ke:.6g} Temp {temp:.3f} K key_min {key:.3f} A "
                    "chains {chains} RMSD {rmsd:.3f} severe {severe} close {close} status {status}".format(
                        step=step,
                        pe=potential,
                        ke=kinetic,
                        temp=temperature,
                        key=retention["key_min_distance_a"],
                        chains=retention["near_chain_count"],
                        rmsd=retention["ligand_rmsd_a"],
                        severe=clashes["severe_lt_1A"],
                        close=clashes["close_lt_1p5A"],
                        status=status,
                    )
                )
            if not nonfinite_detected:
                last_finite_positions = positions
                last_finite_step = step
                last_row = row
            if stopped_early or step == TOTAL_STEPS:
                stop_reason = status if stopped_early else ""
                break
            chunk = min(SAFETY_INTERVAL_STEPS, TOTAL_STEPS - step)
            simulation.step(chunk)
            step += chunk

    final_coords_a = positions_to_angstrom(last_finite_positions)
    final_translation = ligand_image_translation_a(atoms, reference_coords_a, final_coords_a, box_lengths_a)
    final_imaged_coords_a = image_ligand_coords(atoms, final_coords_a, final_translation)
    final_retention = retention_metrics(atoms, reference_coords_a, final_imaged_coords_a, box_lengths_a)
    final_clashes = clash_metrics(atoms, final_imaged_coords_a, box_lengths_a)
    write_pdb_with_coords(
        INPUT_PDB,
        final_imaged_coords_a,
        FINAL_PDB,
        [
            "REMARK ultra-conservative 1 ps 50 K stability check",
            "REMARK not production MD",
            f"REMARK last finite step {last_finite_step}",
            f"REMARK safety stop reason {stop_reason or 'none'}",
        ],
    )
    completed = (not stopped_early) and last_finite_step == TOTAL_STEPS
    result = {
        "completed_1ps": completed,
        "last_finite_step": last_finite_step,
        "total_steps": TOTAL_STEPS,
        "timestep_fs": TIMESTEP_PS * 1000.0,
        "temperature_k": TEMPERATURE_K,
        "temperature_abort_k": TEMPERATURE_ABORT_K,
        "max_temperature_k": max_temperature if math.isfinite(max_temperature) else None,
        "nonfinite_detected": nonfinite_detected,
        "stop_reason": stop_reason,
        "last_row": last_row,
        "final_retention": final_retention,
        "final_clashes": final_clashes,
        "final_pdb": str(FINAL_PDB),
        "retention_csv": str(RETENTION_CSV),
    }
    SUMMARY_JSON.write_text(json.dumps(json_safe(result), indent=2) + "\n", encoding="utf-8")
    print("1 ps check completed:", completed)
    print("last finite step:", last_finite_step)
    print("nonfinite detected:", nonfinite_detected)
    print("stop reason:", stop_reason or "none")
    print("final severe clashes:", final_clashes["severe_lt_1A"])
    print("final close contacts:", final_clashes["close_lt_1p5A"])
    print("final ligand near chains:", final_retention["near_chain_count"])
    print("final ligand key min distance (A):", f"{final_retention['key_min_distance_a']:.4f}")
    print("summary JSON:", SUMMARY_JSON)
    del simulation
    del integrator
    return result


def main() -> int:
    result = run_check()
    return 0 if result["completed_1ps"] and not result["nonfinite_detected"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
