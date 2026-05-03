from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any


COMMON_DIR = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\batch_2.3_6.2_17.2\scripts")
if str(COMMON_DIR) not in sys.path:
    sys.path.insert(0, str(COMMON_DIR))

from generic_pilot_common import (  # noqa: E402
    LIGAND_ID,
    RESNAME,
    build_system,
    choose_opencl_platform,
    clash_metrics,
    coords_have_nan,
    direct_rmsd,
    ensure_dirs as ensure_pilot_dirs,
    estimate_temperature_k,
    is_heavy,
    is_ligand,
    key_distance_metrics,
    make_config,
    near_chain_count,
    pathset,
    pbc_restraint_force,
    positions_to_angstrom,
    r2_restraint_assignment,
    read_gro_box_lengths_a,
    read_pdb_atoms,
    write_csv,
    write_json,
    write_pdb,
)


TOTAL_STEPS_100PS = 200000
TIMESTEP_100PS = 0.0005
REPORT_INTERVAL_STEPS = 2000  # 1 ps
DCD_INTERVAL_STEPS = 4000  # 2 ps
SAFETY_INTERVAL_STEPS = 100
TEMPERATURE_K = 310.0

MONITOR_FIELDS = [
    "step",
    "time_ps",
    "temperature_K",
    "potential_energy_kJ_mol",
    "kinetic_energy_kJ_mol",
    "total_energy_kJ_mol",
    "ligand_distance_TYR13_A",
    "ligand_distance_VAL17_A",
    "ligand_distance_SER20_A",
    "ligand_min_distance_to_key_residues_A",
    "ligand_near_chains_count",
    "ligand_heavy_rmsd_vs_100ps_input_A",
    "ligand_heavy_rmsd_vs_original_R2_minimized_input_A",
    "severe_lt_1A",
    "close_lt_1p5A",
    "status",
]


def build_continuation_note(continued_from_checkpoint: bool, fallback_used: bool, fallback_reason: str) -> str:
    lines = [f"Continued from checkpoint: {continued_from_checkpoint}", f"Fallback used: {fallback_used}"]
    if fallback_reason:
        lines.append(f"Fallback reason: {fallback_reason}")
    return "\n".join(lines)


def md50ps_input_paths(config) -> dict[str, Path]:
    root = config.pilot_root / "md50ps_R2"
    return {
        "root": root,
        "final_pdb": root / "outputs" / f"TRKB_{config.system_name}_{RESNAME}_{LIGAND_ID}_50ps_R2_final.pdb",
        "state_xml": root / "outputs" / f"TRKB_{config.system_name}_{RESNAME}_{LIGAND_ID}_50ps_R2_state.xml",
        "checkpoint": root / "outputs" / f"TRKB_{config.system_name}_{RESNAME}_{LIGAND_ID}_50ps_R2.chk",
        "result_json": root / "reports" / f"md50ps_{config.system_name}_R2_result.json",
    }


def md100ps_paths(config) -> dict[str, Path]:
    root = config.pilot_root / "md100ps_R2"
    return {
        "root": root,
        "scripts": root / "scripts",
        "logs": root / "logs",
        "outputs": root / "outputs",
        "reports": root / "reports",
        "log": root / "logs" / f"md100ps_{config.system_name}_R2.log",
        "monitor_csv": root / "reports" / f"md100ps_{config.system_name}_R2_monitor.csv",
        "summary_md": root / "reports" / f"md100ps_{config.system_name}_R2_summary.md",
        "final_pdb": root / "outputs" / f"TRKB_{config.system_name}_{RESNAME}_{LIGAND_ID}_100ps_R2_final.pdb",
        "dcd": root / "outputs" / f"TRKB_{config.system_name}_{RESNAME}_{LIGAND_ID}_100ps_R2.dcd",
        "state_xml": root / "outputs" / f"TRKB_{config.system_name}_{RESNAME}_{LIGAND_ID}_100ps_R2_state.xml",
        "checkpoint": root / "outputs" / f"TRKB_{config.system_name}_{RESNAME}_{LIGAND_ID}_100ps_R2.chk",
        "result_json": root / "reports" / f"md100ps_{config.system_name}_R2_result.json",
    }


def ensure_dirs(config) -> None:
    ensure_pilot_dirs(config)
    for path in md100ps_paths(config).values():
        if isinstance(path, Path) and path.suffix == "":
            path.mkdir(parents=True, exist_ok=True)


def summarize_numeric(rows: list[dict[str, Any]], field: str) -> dict[str, float]:
    values = []
    for row in rows:
        try:
            value = float(row.get(field))
        except Exception:
            continue
        if math.isfinite(value):
            values.append(value)
    if not values:
        return {"initial": float("nan"), "final": float("nan"), "min": float("nan"), "max": float("nan"), "mean": float("nan")}
    return {
        "initial": values[0],
        "final": values[-1],
        "min": min(values),
        "max": max(values),
        "mean": sum(values) / len(values),
    }


def initialize_state(config, simulation, pdb) -> tuple[bool, bool, str]:
    input_paths = md50ps_input_paths(config)
    checkpoint_in = input_paths["checkpoint"]
    state_xml_in = input_paths["state_xml"]
    if checkpoint_in.exists():
        try:
            simulation.loadCheckpoint(str(checkpoint_in))
            return True, False, ""
        except Exception as exc:
            reason = f"checkpoint load failed: {exc!r}"
    elif state_xml_in.exists():
        reason = "checkpoint missing; state XML exists but checkpoint is preferred"
    else:
        reason = "checkpoint and state XML missing"
    from openmm.unit import kelvin

    simulation.context.setPositions(pdb.positions)
    simulation.context.setVelocitiesToTemperature(TEMPERATURE_K * kelvin, config.seed + 100)
    return False, True, reason


def write_summary(config, result: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    paths = md100ps_paths(config)
    input_paths = md50ps_input_paths(config)
    temp = summarize_numeric(rows, "temperature_K")
    pe = summarize_numeric(rows, "potential_energy_kJ_mol")
    tyr = summarize_numeric(rows, "ligand_distance_TYR13_A")
    val = summarize_numeric(rows, "ligand_distance_VAL17_A")
    ser = summarize_numeric(rows, "ligand_distance_SER20_A")
    key = summarize_numeric(rows, "ligand_min_distance_to_key_residues_A")
    rmsd_input = summarize_numeric(rows, "ligand_heavy_rmsd_vs_100ps_input_A")
    rmsd_min = summarize_numeric(rows, "ligand_heavy_rmsd_vs_original_R2_minimized_input_A")
    near_chains_all_2 = bool(rows) and all(int(row.get("ligand_near_chains_count", 0)) == 2 for row in rows)
    severe_seen = any(int(row.get("severe_lt_1A", 0)) > 0 for row in rows)
    close_seen = any(int(row.get("close_lt_1p5A", 0)) > 0 for row in rows)
    lines = [
        f"# {config.system_name} 100 ps R2 Restrained MD Pilot Summary",
        "",
        "## Purpose",
        "",
        f"Run a 100 ps restrained MD pilot continuation for {LIGAND_ID} / {RESNAME} / {config.system_name} after the successful 50 ps R2 restrained MD pilot.",
        "",
        "## Continuation",
        "",
        f"- Continued from checkpoint: {result['continued_from_checkpoint']}",
        f"- 50 ps checkpoint input: `{input_paths['checkpoint']}`",
        f"- 50 ps state XML input: `{input_paths['state_xml']}`",
        f"- 50 ps final PDB input: `{input_paths['final_pdb']}`",
        f"- Fallback used: {result['fallback_used']}",
        f"- Fallback reason: {result['fallback_reason'] or 'none'}",
        "",
        "## Method",
        "",
        f"- OpenCL used: {result['platform'] == 'OpenCL'}",
        f"- Platform: {result['platform']}",
        f"- Timestep: {TIMESTEP_100PS} ps (0.5 fs)",
        "- Integrator: LangevinMiddleIntegrator",
        "- Nonbonded method: PME",
        "- Constraints: HBonds",
        "- rigidWater: True",
        "- Temperature: 310 K",
        "- R2 restraints: protein backbone medium, ligand heavy atoms weak, lipid/cholesterol heavy atoms extremely weak, water/ions unrestrained.",
        "- PBC-aware restraint expression: `0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2`",
        "- Monitoring cadence: 1 ps",
        "- DCD cadence: 2 ps",
        "- This is restrained MD pilot work, not production MD.",
        "",
        "## Outcome",
        "",
        f"- Completed 100 ps: {result['completed_100ps']}",
        f"- Last step: {result['last_step']}/{TOTAL_STEPS_100PS}",
        f"- Stop reason: {result['stop_reason']}",
        f"- NaN detected: {result['has_nan']}",
        f"- Temperature K initial/final/min/max/mean: {temp['initial']:.6g}/{temp['final']:.6g}/{temp['min']:.6g}/{temp['max']:.6g}/{temp['mean']:.6g}",
        f"- Potential energy kJ/mol initial/final/mean: {pe['initial']:.12g}/{pe['final']:.12g}/{pe['mean']:.12g}",
        "",
        "## Ligand Retention",
        "",
        f"- TYR13 distance A initial/final/min/max/mean: {tyr['initial']:.6g}/{tyr['final']:.6g}/{tyr['min']:.6g}/{tyr['max']:.6g}/{tyr['mean']:.6g}",
        f"- VAL17 distance A initial/final/min/max/mean: {val['initial']:.6g}/{val['final']:.6g}/{val['min']:.6g}/{val['max']:.6g}/{val['mean']:.6g}",
        f"- SER20 distance A initial/final/min/max/mean: {ser['initial']:.6g}/{ser['final']:.6g}/{ser['min']:.6g}/{ser['max']:.6g}/{ser['mean']:.6g}",
        f"- Key-residue min distance A initial/final/min/max/mean: {key['initial']:.6g}/{key['final']:.6g}/{key['min']:.6g}/{key['max']:.6g}/{key['mean']:.6g}",
        f"- Ligand near chains always 2: {near_chains_all_2}",
        f"- Ligand RMSD vs 100 ps input A initial/final/max: {rmsd_input['initial']:.6g}/{rmsd_input['final']:.6g}/{rmsd_input['max']:.6g}",
        f"- Ligand RMSD vs original R2 minimized input A initial/final/max: {rmsd_min['initial']:.6g}/{rmsd_min['final']:.6g}/{rmsd_min['max']:.6g}",
        "",
        "## Clash Monitoring",
        "",
        f"- Severe clashes appeared: {severe_seen}",
        f"- Close contacts appeared: {close_seen}",
        f"- Final severe/close contacts: {result['final_row'].get('severe_lt_1A', 'nan')}/{result['final_row'].get('close_lt_1p5A', 'nan')}",
        "",
        "## Recommendation",
        "",
        f"- Recommend 20chol vs 40chol comparison analysis: {result['recommend_comparison_analysis']}",
        "- This stage is still restrained MD pilot work, not production MD.",
    ]
    paths["summary_md"].write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_100ps(config) -> dict[str, Any]:
    from openmm import LangevinMiddleIntegrator, XmlSerializer
    from openmm.app import DCDReporter, HBonds, Simulation
    from openmm.unit import kelvin, kilojoules_per_mole, picoseconds

    ensure_dirs(config)
    paths = md100ps_paths(config)
    pilot_paths = pathset(config)
    input_paths = md50ps_input_paths(config)
    input_psf = pilot_paths["cleaned_psf"]
    input_pdb = input_paths["final_pdb"]
    minimized_input_pdb = pilot_paths["minimized_pdb"]

    atoms = read_pdb_atoms(input_pdb)
    psf, pdb, system = build_system(config, input_psf, input_pdb, constraints=HBonds, rigid_water=True)
    restraints = r2_restraint_assignment(atoms)
    system.addForce(pbc_restraint_force(restraints, pdb.positions))
    platform = choose_opencl_platform()
    integrator = LangevinMiddleIntegrator(TEMPERATURE_K * kelvin, 10.0 / picoseconds, TIMESTEP_100PS * picoseconds)
    integrator.setRandomNumberSeed(config.seed + 100)
    simulation = Simulation(psf.topology, system, integrator, platform)
    continued_from_checkpoint, fallback_used, fallback_reason = initialize_state(config, simulation, pdb)
    simulation.reporters.append(DCDReporter(str(paths["dcd"]), DCD_INTERVAL_STEPS))

    reference_input_atoms = read_pdb_atoms(input_pdb)
    reference_input_coords = [(atom.x, atom.y, atom.z) for atom in reference_input_atoms]
    reference_minimized_atoms = read_pdb_atoms(minimized_input_pdb)
    reference_minimized_coords = [(atom.x, atom.y, atom.z) for atom in reference_minimized_atoms]
    ligand_heavy = [index for index, atom in enumerate(atoms) if is_ligand(atom) and is_heavy(atom)]
    box = read_gro_box_lengths_a(config.gro_box)

    rows: list[dict[str, Any]] = []
    last_positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=False)
    last_row: dict[str, Any] | None = None
    stopped = False
    stop_reason = "none"
    global_step = 0
    severe_streak = 0
    far_streak = 0
    has_nan = False
    while global_step < TOTAL_STEPS_100PS:
        chunk = min(SAFETY_INTERVAL_STEPS, TOTAL_STEPS_100PS - global_step)
        simulation.step(chunk)
        global_step += chunk
        if global_step % REPORT_INTERVAL_STEPS != 0 and global_step != TOTAL_STEPS_100PS:
            continue
        state = simulation.context.getState(getEnergy=True, getPositions=True, getVelocities=True)
        pe = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        ke = state.getKineticEnergy().value_in_unit(kilojoules_per_mole)
        temp = estimate_temperature_k(ke, system)
        coords = positions_to_angstrom(state.getPositions(asNumpy=False))
        key_dist = key_distance_metrics(atoms, coords, box)
        near_chains = near_chain_count(atoms, coords, box)
        clashes = clash_metrics(atoms, coords, box)
        rmsd_input = direct_rmsd(ligand_heavy, reference_input_coords, coords)
        rmsd_minimized = direct_rmsd(ligand_heavy, reference_minimized_coords, coords)
        statuses = []
        if not all(math.isfinite(value) for value in [pe, ke, temp]) or coords_have_nan(coords):
            statuses.append("nonfinite_state")
            has_nan = True
        if math.isfinite(temp) and temp > 500.0:
            statuses.append("temperature_above_500K")
        if near_chains == 0:
            statuses.append("ligand_near_chains_0")
        if key_dist["key_min"] > 10.0:
            far_streak += 1
            statuses.append("warning_all_key_distances_gt_10A")
            if far_streak >= 2:
                statuses.append("persistent_all_key_distances_gt_10A")
        else:
            far_streak = 0
        if clashes["severe_lt_1A"] > 0:
            severe_streak += 1
            statuses.append("warning_severe_clash_gt_0")
            if severe_streak >= 2:
                statuses.append("persistent_severe_clash_gt_0")
        else:
            severe_streak = 0
        if not statuses:
            statuses.append("ok")
        row = {
            "step": global_step,
            "time_ps": f"{global_step * TIMESTEP_100PS:.6g}",
            "temperature_K": f"{temp:.6g}",
            "potential_energy_kJ_mol": f"{pe:.12g}",
            "kinetic_energy_kJ_mol": f"{ke:.12g}",
            "total_energy_kJ_mol": f"{pe + ke:.12g}",
            "ligand_distance_TYR13_A": f"{key_dist['TYR13']:.6g}",
            "ligand_distance_VAL17_A": f"{key_dist['VAL17']:.6g}",
            "ligand_distance_SER20_A": f"{key_dist['SER20']:.6g}",
            "ligand_min_distance_to_key_residues_A": f"{key_dist['key_min']:.6g}",
            "ligand_near_chains_count": near_chains,
            "ligand_heavy_rmsd_vs_100ps_input_A": f"{rmsd_input:.6g}",
            "ligand_heavy_rmsd_vs_original_R2_minimized_input_A": f"{rmsd_minimized:.6g}",
            "severe_lt_1A": clashes["severe_lt_1A"],
            "close_lt_1p5A": clashes["close_lt_1p5A"],
            "status": ";".join(statuses),
        }
        rows.append(row)
        last_row = row
        last_positions = state.getPositions(asNumpy=False)
        if any(
            token in row["status"]
            for token in [
                "nonfinite_state",
                "temperature_above_500K",
                "ligand_near_chains_0",
                "persistent_all_key_distances_gt_10A",
                "persistent_severe_clash_gt_0",
            ]
        ):
            stopped = True
            stop_reason = row["status"]
            break
    final_state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)
    paths["state_xml"].write_text(XmlSerializer.serialize(final_state), encoding="utf-8")
    simulation.saveCheckpoint(str(paths["checkpoint"]))
    write_csv(paths["monitor_csv"], rows, MONITOR_FIELDS)
    write_pdb(
        input_pdb,
        positions_to_angstrom(last_positions),
        paths["final_pdb"],
        [
            f"REMARK {LIGAND_ID} {config.system_name} 100 ps R2 restrained MD final coordinates",
            "REMARK not production MD",
            f"REMARK completed_100ps {((not stopped) and global_step == TOTAL_STEPS_100PS)}",
            f"REMARK continued_from_checkpoint {continued_from_checkpoint}",
            f"REMARK fallback_used {fallback_used}",
            f"REMARK stop_reason {stop_reason}",
        ],
    )
    completed = (not stopped) and global_step == TOTAL_STEPS_100PS
    final_row = last_row or {}
    recommend_comparison_analysis = (
        completed
        and not has_nan
        and int(final_row.get("ligand_near_chains_count", 0)) == 2
        and int(final_row.get("severe_lt_1A", 1)) == 0
        and int(final_row.get("close_lt_1p5A", 1)) == 0
    )
    rmsd_stats_input = summarize_numeric(rows, "ligand_heavy_rmsd_vs_100ps_input_A")
    rmsd_stats_min = summarize_numeric(rows, "ligand_heavy_rmsd_vs_original_R2_minimized_input_A")
    result = {
        "completed_100ps": completed,
        "continued_from_checkpoint": continued_from_checkpoint,
        "fallback_used": fallback_used,
        "fallback_reason": fallback_reason,
        "last_step": global_step,
        "stop_reason": stop_reason,
        "has_nan": has_nan,
        "platform": platform.getName(),
        "final_row": final_row,
        "ligand_rmsd_final_A": float(final_row["ligand_heavy_rmsd_vs_100ps_input_A"]) if final_row else float("nan"),
        "ligand_rmsd_max_A": rmsd_stats_input["max"],
        "ligand_rmsd_vs_original_minimized_final_A": float(final_row["ligand_heavy_rmsd_vs_original_R2_minimized_input_A"]) if final_row else float("nan"),
        "ligand_rmsd_vs_original_minimized_max_A": rmsd_stats_min["max"],
        "recommend_comparison_analysis": recommend_comparison_analysis,
    }
    write_summary(config, result, rows)
    write_json(paths["result_json"], result)
    return result
