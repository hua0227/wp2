from __future__ import annotations

import math
import sys
import traceback
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from typing import Any


PILOT_ROOT = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_40chol")
COMMON_SCRIPT_DIR = PILOT_ROOT / "scripts"
if str(COMMON_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(COMMON_SCRIPT_DIR))

from common_40chol import (  # noqa: E402
    CLEANED_PSF,
    GRO_BOX,
    LIGAND_STR,
    MINIMIZED_PDB,
    PROTEIN_BACKBONE_NAMES,
    TOPPAR_DIR,
    PdbAtom,
    build_system,
    choose_opencl_platform,
    distance_a,
    is_cholesterol,
    is_heavy,
    is_ligand,
    is_lipid,
    is_protein,
    positions_to_angstrom,
    read_gro_box_lengths_a,
    read_pdb_atoms,
    write_csv,
    write_pdb,
)


ROOT = PILOT_ROOT / "md100ps_R2"
SCRIPTS_DIR = ROOT / "scripts"
LOGS_DIR = ROOT / "logs"
OUTPUTS_DIR = ROOT / "outputs"
REPORTS_DIR = ROOT / "reports"

INPUT_PSF = CLEANED_PSF
INPUT_PDB = PILOT_ROOT / "md50ps_R2" / "outputs" / "TRKB_40chol_L002_14.2_50ps_R2_final.pdb"
CHECKPOINT_IN = PILOT_ROOT / "md50ps_R2" / "outputs" / "TRKB_40chol_L002_14.2_50ps_R2.chk"
STATE_XML_IN = PILOT_ROOT / "md50ps_R2" / "outputs" / "TRKB_40chol_L002_14.2_50ps_R2_state.xml"
R2_MINIMIZED_INPUT_PDB = MINIMIZED_PDB

LOG_PATH = LOGS_DIR / "md100ps_40chol_R2.log"
MONITOR_CSV = REPORTS_DIR / "md100ps_40chol_R2_monitor.csv"
SUMMARY_MD = REPORTS_DIR / "md100ps_40chol_R2_summary.md"
FINAL_PDB = OUTPUTS_DIR / "TRKB_40chol_L002_14.2_100ps_R2_final.pdb"
DCD_PATH = OUTPUTS_DIR / "TRKB_40chol_L002_14.2_100ps_R2.dcd"
STATE_XML = OUTPUTS_DIR / "TRKB_40chol_L002_14.2_100ps_R2_state.xml"
CHECKPOINT_PATH = OUTPUTS_DIR / "TRKB_40chol_L002_14.2_100ps_R2.chk"

PBC_RESTRAINT_EXPRESSION = "0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2"
PROTEIN_BACKBONE_K = 500.0
LIGAND_HEAVY_K = 25.0
LIPID_CHOL_HEAVY_K = 5.0

TIMESTEP_PS = 0.0005
TOTAL_STEPS = 200000
REPORT_INTERVAL_STEPS = 2000  # 1 ps
SAFETY_INTERVAL_STEPS = 100
DCD_INTERVAL_STEPS = 4000  # 2 ps
TEMPERATURE_K = 310.0
RETENTION_CUTOFF_A = 6.0
KEY_DISTANCE_ABORT_A = 10.0
VAL17_WARNING_A = 10.0
SER20_WARNING_A = 10.0
TEMPERATURE_ABORT_K = 500.0
MAX_RMSD_RECOMMEND_A = 6.0

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
    "ligand_heavy_rmsd_vs_R2_minimized_A",
    "severe_lt_1A",
    "close_lt_1p5A",
    "status",
]


class Tee:
    def __init__(self, *streams: Any) -> None:
        self.streams = streams

    def write(self, data: str) -> int:
        for stream in self.streams:
            stream.write(data)
            stream.flush()
        return len(data)

    def flush(self) -> None:
        for stream in self.streams:
            stream.flush()


def ensure_dirs() -> None:
    for path in [SCRIPTS_DIR, LOGS_DIR, OUTPUTS_DIR, REPORTS_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def restraint_map(atoms: list[PdbAtom]) -> dict[int, float]:
    mapping: dict[int, float] = {}
    for index, atom in enumerate(atoms):
        k_value = 0.0
        if is_protein(atom) and atom.name in PROTEIN_BACKBONE_NAMES:
            k_value = PROTEIN_BACKBONE_K
        elif is_ligand(atom) and is_heavy(atom):
            k_value = LIGAND_HEAVY_K
        elif (is_lipid(atom) or is_cholesterol(atom)) and is_heavy(atom):
            k_value = LIPID_CHOL_HEAVY_K
        if k_value > 0.0:
            mapping[index] = k_value
    return mapping


def add_r2_pbc_restraints(system: Any, atoms: list[PdbAtom], positions: Any) -> dict[str, int]:
    from openmm import CustomExternalForce
    from openmm.unit import nanometer

    force = CustomExternalForce(PBC_RESTRAINT_EXPRESSION)
    force.addPerParticleParameter("k")
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    counts = {"protein_backbone": 0, "ligand_heavy": 0, "lipid_cholesterol_heavy": 0, "total": 0}
    mapping = restraint_map(atoms)
    for index, k_value in mapping.items():
        pos = positions[index].value_in_unit(nanometer)
        force.addParticle(index, [k_value, pos.x, pos.y, pos.z])
        atom = atoms[index]
        if is_protein(atom):
            counts["protein_backbone"] += 1
        elif is_ligand(atom):
            counts["ligand_heavy"] += 1
        else:
            counts["lipid_cholesterol_heavy"] += 1
        counts["total"] += 1
    system.addForce(force)
    return counts


def estimate_temperature_k(kinetic_energy_kj_mol: float, system: Any) -> float:
    dof = 3 * system.getNumParticles() - system.getNumConstraints()
    for index in range(system.getNumForces()):
        if system.getForce(index).__class__.__name__ == "CMMotionRemover":
            dof -= 3
            break
    return 2 * kinetic_energy_kj_mol / (dof * 0.00831446261815324)


def key_distance_metrics(atoms: list[PdbAtom], coords: list[tuple[float, float, float]], box: tuple[float, float, float]) -> dict[str, float]:
    ligand_heavy = [index for index, atom in enumerate(atoms) if is_ligand(atom) and is_heavy(atom)]
    residue_map = {"TYR13": "434", "VAL17": "438", "SER20": "441"}
    out: dict[str, float] = {}
    for label, resid in residue_map.items():
        key_indices = [
            index
            for index, atom in enumerate(atoms)
            if atom.segid in {"PROA", "PROB"} and atom.resid == resid and is_heavy(atom)
        ]
        out[label] = min(distance_a(coords[i], coords[j], box) for i in ligand_heavy for j in key_indices)
    out["key_min"] = min(out.values())
    return out


def ligand_rmsd(atoms: list[PdbAtom], reference: list[tuple[float, float, float]], coords: list[tuple[float, float, float]]) -> float:
    ligand_heavy = [index for index, atom in enumerate(atoms) if is_ligand(atom) and is_heavy(atom)]
    return math.sqrt(sum(distance_a(reference[index], coords[index]) ** 2 for index in ligand_heavy) / len(ligand_heavy))


def near_chain_count(atoms: list[PdbAtom], coords: list[tuple[float, float, float]], box: tuple[float, float, float]) -> int:
    ligand_heavy = [index for index, atom in enumerate(atoms) if is_ligand(atom) and is_heavy(atom)]
    count = 0
    for segid in sorted({atom.segid for atom in atoms if is_protein(atom)}):
        indices = [index for index, atom in enumerate(atoms) if atom.segid == segid and is_protein(atom) and is_heavy(atom)]
        min_dist = min(distance_a(coords[i], coords[j], box) for i in ligand_heavy for j in indices)
        if min_dist <= RETENTION_CUTOFF_A:
            count += 1
    return count


def clash_metrics(atoms: list[PdbAtom], coords: list[tuple[float, float, float]], box: tuple[float, float, float]) -> dict[str, int]:
    ligand_heavy = [index for index, atom in enumerate(atoms) if is_ligand(atom) and is_heavy(atom)]
    other_heavy = [index for index, atom in enumerate(atoms) if not is_ligand(atom) and is_heavy(atom)]
    severe = 0
    close = 0
    for i in ligand_heavy:
        for j in other_heavy:
            dist = distance_a(coords[i], coords[j], box)
            if dist < 1.5:
                close += 1
                if dist < 1.0:
                    severe += 1
    return {"severe_lt_1A": severe, "close_lt_1p5A": close}


def coords_finite(coords: list[tuple[float, float, float]]) -> bool:
    return all(math.isfinite(value) for xyz in coords for value in xyz)


def read_reference_coords(path: Path) -> list[tuple[float, float, float]]:
    atoms = read_pdb_atoms(path)
    return [(atom.x, atom.y, atom.z) for atom in atoms]


def monitor_row(
    step: int,
    temperature: float,
    potential: float,
    kinetic: float,
    key_distances: dict[str, float],
    chains: int,
    rmsd_vs_input: float,
    rmsd_vs_minimized: float,
    clashes: dict[str, int],
    status: str,
) -> dict[str, Any]:
    return {
        "step": step,
        "time_ps": f"{step * TIMESTEP_PS:.6g}",
        "temperature_K": f"{temperature:.6g}",
        "potential_energy_kJ_mol": f"{potential:.12g}",
        "kinetic_energy_kJ_mol": f"{kinetic:.12g}",
        "total_energy_kJ_mol": f"{potential + kinetic:.12g}",
        "ligand_distance_TYR13_A": f"{key_distances['TYR13']:.6g}",
        "ligand_distance_VAL17_A": f"{key_distances['VAL17']:.6g}",
        "ligand_distance_SER20_A": f"{key_distances['SER20']:.6g}",
        "ligand_min_distance_to_key_residues_A": f"{key_distances['key_min']:.6g}",
        "ligand_near_chains_count": chains,
        "ligand_heavy_rmsd_vs_100ps_input_A": f"{rmsd_vs_input:.6g}",
        "ligand_heavy_rmsd_vs_R2_minimized_A": f"{rmsd_vs_minimized:.6g}",
        "severe_lt_1A": clashes["severe_lt_1A"],
        "close_lt_1p5A": clashes["close_lt_1p5A"],
        "status": status,
    }


def to_float(value: Any, default: float = float("nan")) -> float:
    try:
        return float(value)
    except Exception:
        return default


def to_int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except Exception:
        return default


def summarize_numeric(rows: list[dict[str, Any]], field: str) -> dict[str, float]:
    values = [to_float(row.get(field)) for row in rows]
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return {"initial": float("nan"), "final": float("nan"), "min": float("nan"), "max": float("nan"), "mean": float("nan")}
    return {
        "initial": finite[0],
        "final": finite[-1],
        "min": min(finite),
        "max": max(finite),
        "mean": sum(finite) / len(finite),
    }


def fmt(value: Any, digits: int = 6) -> str:
    number = to_float(value)
    return f"{number:.{digits}g}" if math.isfinite(number) else "nan"


def write_summary(result: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    temp = summarize_numeric(rows, "temperature_K")
    pe = summarize_numeric(rows, "potential_energy_kJ_mol")
    tyr = summarize_numeric(rows, "ligand_distance_TYR13_A")
    val = summarize_numeric(rows, "ligand_distance_VAL17_A")
    ser = summarize_numeric(rows, "ligand_distance_SER20_A")
    key = summarize_numeric(rows, "ligand_min_distance_to_key_residues_A")
    rmsd = summarize_numeric(rows, "ligand_heavy_rmsd_vs_100ps_input_A")
    near_chains_all_2 = bool(rows) and all(to_int(row.get("ligand_near_chains_count"), -1) == 2 for row in rows)
    severe_seen = any(to_int(row.get("severe_lt_1A"), 0) > 0 for row in rows)
    close_seen = any(to_int(row.get("close_lt_1p5A"), 0) > 0 for row in rows)
    potential_finite = bool(rows) and all(math.isfinite(to_float(row.get("potential_energy_kJ_mol"))) for row in rows)
    lines = [
        "# 40chol 100 ps R2 Restrained MD Pilot Summary",
        "",
        "## Purpose",
        "",
        "Run a 100 ps restrained MD pilot continuation for 14.2 / L002 / 40chol after the successful 50 ps R2 restrained MD pilot.",
        "",
        "## Continuation",
        "",
        f"- Continued from checkpoint: {result.get('continued_from_checkpoint')}",
        f"- Checkpoint input: `{CHECKPOINT_IN}`",
        f"- State XML input: `{STATE_XML_IN}`",
        f"- Fallback used: {result.get('fallback_used')}",
        f"- Fallback reason: {result.get('fallback_reason') or 'none'}",
        "",
        "## Method",
        "",
        f"- OpenCL used: {result.get('platform') == 'OpenCL'}",
        f"- Platform: {result.get('platform', 'not reached')}",
        f"- Timestep: {TIMESTEP_PS} ps (0.5 fs)",
        "- Integrator: LangevinMiddleIntegrator",
        "- Nonbonded method: PME",
        "- Constraints: HBonds",
        "- rigidWater: True",
        "- Temperature: 310 K",
        "- R2 restraints: protein backbone medium, ligand heavy atoms weak, lipid/cholesterol heavy atoms extremely weak, water/ions unrestrained.",
        f"- PBC-aware restraint expression: `{PBC_RESTRAINT_EXPRESSION}`",
        "- This is restrained MD pilot work, not production MD.",
        "",
        "## Outcome",
        "",
        f"- Completed 100 ps: {result.get('completed_100ps')}",
        f"- Last step: {result.get('last_step')}/{TOTAL_STEPS}",
        f"- Stop reason: {result.get('stop_reason') or 'none'}",
        f"- NaN detected: {result.get('has_nan')}",
        f"- Potential energy finite across monitor rows: {potential_finite}",
        f"- Temperature K initial/final/min/max/mean: {fmt(temp['initial'])}/{fmt(temp['final'])}/{fmt(temp['min'])}/{fmt(temp['max'])}/{fmt(temp['mean'])}",
        f"- Potential energy kJ/mol initial/final/mean: {fmt(pe['initial'])}/{fmt(pe['final'])}/{fmt(pe['mean'])}",
        "",
        "## Ligand Retention",
        "",
        f"- TYR13 distance A initial/final/min/max/mean: {fmt(tyr['initial'])}/{fmt(tyr['final'])}/{fmt(tyr['min'])}/{fmt(tyr['max'])}/{fmt(tyr['mean'])}",
        f"- VAL17 distance A initial/final/min/max/mean: {fmt(val['initial'])}/{fmt(val['final'])}/{fmt(val['min'])}/{fmt(val['max'])}/{fmt(val['mean'])}",
        f"- SER20 distance A initial/final/min/max/mean: {fmt(ser['initial'])}/{fmt(ser['final'])}/{fmt(ser['min'])}/{fmt(ser['max'])}/{fmt(ser['mean'])}",
        f"- Key-residue min distance A initial/final/min/max/mean: {fmt(key['initial'])}/{fmt(key['final'])}/{fmt(key['min'])}/{fmt(key['max'])}/{fmt(key['mean'])}",
        f"- Ligand near chains always 2: {near_chains_all_2}",
        f"- Ligand RMSD vs 100 ps input A initial/final/max: {fmt(rmsd['initial'])}/{fmt(rmsd['final'])}/{fmt(rmsd['max'])}",
        "",
        "## Clash Monitoring",
        "",
        f"- Severe clashes appeared: {severe_seen}",
        f"- Close contacts appeared: {close_seen}",
        f"- Final severe/close contacts: {result.get('final_row', {}).get('severe_lt_1A', 'not available')}/{result.get('final_row', {}).get('close_lt_1p5A', 'not available')}",
        "",
        "## Recommendation",
        "",
        f"- Recommend 14.2 20chol vs 40chol comparison analysis: {result.get('recommend_comparison')}",
        "- This stage is still restrained MD pilot work, not production MD.",
    ]
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def initialize_state(simulation: Any, pdb: Any) -> tuple[bool, bool, str]:
    from openmm.unit import kelvin

    if CHECKPOINT_IN.exists():
        try:
            simulation.loadCheckpoint(str(CHECKPOINT_IN))
            return True, False, ""
        except Exception as exc:
            reason = f"checkpoint load failed: {exc!r}"
            print(f"WARNING: {reason}")
    elif STATE_XML_IN.exists():
        reason = "checkpoint missing"
        print("WARNING: checkpoint missing; state XML is available but checkpoint is preferred. Using PDB fallback.")
    else:
        reason = "checkpoint and state XML missing"
        print(f"WARNING: {reason}")
    simulation.context.setPositions(pdb.positions)
    simulation.context.setVelocitiesToTemperature(TEMPERATURE_K * kelvin, 142500)
    return False, True, reason


def run_md100ps() -> dict[str, Any]:
    import openmm
    from openmm import LangevinMiddleIntegrator, Platform
    from openmm.app import DCDReporter, HBonds, Simulation
    from openmm.unit import kelvin, kilojoules_per_mole, picoseconds

    ensure_dirs()
    atoms = read_pdb_atoms(INPUT_PDB)
    psf, pdb, system = build_system(INPUT_PSF, INPUT_PDB, constraints=HBonds, rigid_water=True)
    if len(atoms) != len(pdb.positions):
        raise ValueError(f"Atom count mismatch: parsed {len(atoms)}, OpenMM positions {len(pdb.positions)}")
    restraint_counts = add_r2_pbc_restraints(system, atoms, pdb.positions)
    box = read_gro_box_lengths_a()
    reference_100ps = positions_to_angstrom(pdb.positions)
    reference_minimized = read_reference_coords(R2_MINIMIZED_INPUT_PDB)

    print("100 ps R2 restrained MD pilot continuation for 14.2 / L002 / 40chol")
    print("This is not production MD. No CUDA will be used.")
    print(f"Input PSF: {INPUT_PSF}")
    print(f"Input PDB: {INPUT_PDB}")
    print(f"Checkpoint input: {CHECKPOINT_IN}")
    print(f"State XML input: {STATE_XML_IN}")
    print(f"GRO box: {GRO_BOX}")
    print(f"Toppar dir: {TOPPAR_DIR}")
    print(f"Ligand str: {LIGAND_STR}")
    print(f"OpenMM version: {openmm.version.version}")
    print("Available OpenMM platforms:")
    for index in range(Platform.getNumPlatforms()):
        print(f"- {Platform.getPlatform(index).getName()}")
    platform = choose_opencl_platform(Platform)
    print(f"Selected OpenMM platform: {platform.getName()}")
    print(f"R2 PBC-aware restraint expression: {PBC_RESTRAINT_EXPRESSION}")
    print("R2 restraint counts:")
    for key, value in restraint_counts.items():
        print(f"- {key}: {value}")

    integrator = LangevinMiddleIntegrator(TEMPERATURE_K * kelvin, 10.0 / picoseconds, TIMESTEP_PS * picoseconds)
    integrator.setRandomNumberSeed(142500)
    simulation = Simulation(psf.topology, system, integrator, platform)
    continued_from_checkpoint, fallback_used, fallback_reason = initialize_state(simulation, pdb)
    print(f"continued_from_checkpoint: {continued_from_checkpoint}")
    print(f"fallback_used: {fallback_used}")
    if fallback_reason:
        print(f"fallback_reason: {fallback_reason}")
    simulation.reporters.append(DCDReporter(str(DCD_PATH), DCD_INTERVAL_STEPS))

    rows: list[dict[str, Any]] = []
    last_positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=False)
    last_row: dict[str, Any] = {}
    has_nan = False
    stopped = False
    stop_reason = ""
    severe_warning_streak = 0
    global_step = 0

    def sample(force_write: bool) -> bool:
        nonlocal last_positions, last_row, has_nan, stop_reason, severe_warning_streak
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        potential = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        kinetic = state.getKineticEnergy().value_in_unit(kilojoules_per_mole)
        temperature = estimate_temperature_k(kinetic, system)
        positions = state.getPositions(asNumpy=False)
        coords = positions_to_angstrom(positions)
        finite = all(math.isfinite(value) for value in [potential, kinetic, temperature]) and coords_finite(coords)
        statuses: list[str] = []
        abort = False
        if finite:
            key_distances = key_distance_metrics(atoms, coords, box)
            chains = near_chain_count(atoms, coords, box)
            rmsd_vs_input = ligand_rmsd(atoms, reference_100ps, coords)
            rmsd_vs_minimized = ligand_rmsd(atoms, reference_minimized, coords)
            clashes = clash_metrics(atoms, coords, box)
            last_positions = positions
        else:
            key_distances = {"TYR13": float("nan"), "VAL17": float("nan"), "SER20": float("nan"), "key_min": float("nan")}
            chains = 0
            rmsd_vs_input = float("nan")
            rmsd_vs_minimized = float("nan")
            clashes = {"severe_lt_1A": -1, "close_lt_1p5A": -1}
            has_nan = True
            abort = True
            statuses.append("nonfinite_energy_temperature_or_coordinates")
        if math.isfinite(temperature) and temperature > TEMPERATURE_ABORT_K:
            statuses.append("temperature_above_500K")
            abort = True
        if chains == 0:
            statuses.append("ligand_near_chains_0")
            abort = True
        if all(math.isfinite(key_distances[label]) and key_distances[label] > KEY_DISTANCE_ABORT_A for label in ["TYR13", "VAL17", "SER20"]):
            statuses.append("all_key_distances_gt_10A")
            abort = True
        if math.isfinite(key_distances["VAL17"]) and key_distances["VAL17"] > VAL17_WARNING_A:
            statuses.append("warning_VAL17_distance_gt_10A")
        if math.isfinite(key_distances["SER20"]) and key_distances["SER20"] > SER20_WARNING_A:
            statuses.append("warning_SER20_distance_gt_10A")
        if clashes["severe_lt_1A"] > 0:
            statuses.append("warning_severe_clash_gt_0")
            severe_warning_streak += 1
            if severe_warning_streak >= 3:
                statuses.append("persistent_severe_clashes")
                abort = True
        else:
            severe_warning_streak = 0
        if not statuses:
            statuses.append("ok")
        row = monitor_row(
            global_step,
            temperature,
            potential,
            kinetic,
            key_distances,
            chains,
            rmsd_vs_input,
            rmsd_vs_minimized,
            clashes,
            ";".join(statuses),
        )
        last_row = row
        if force_write or abort:
            rows.append(row)
            print(
                "step {step} time {time} ps temp {temp:.3f} K PE {pe:.6g} "
                "TYR/VAL/SER {tyr:.3f}/{val:.3f}/{ser:.3f} A chains {chains} "
                "RMSD100 {rmsd100:.3f} A RMSDmin {rmsdmin:.3f} A severe/close {severe}/{close} status {status}".format(
                    step=global_step,
                    time=row["time_ps"],
                    temp=temperature,
                    pe=potential,
                    tyr=key_distances["TYR13"],
                    val=key_distances["VAL17"],
                    ser=key_distances["SER20"],
                    chains=chains,
                    rmsd100=rmsd_vs_input,
                    rmsdmin=rmsd_vs_minimized,
                    severe=clashes["severe_lt_1A"],
                    close=clashes["close_lt_1p5A"],
                    status=row["status"],
                )
            )
        if abort:
            stop_reason = row["status"]
        return abort

    sample(True)
    while global_step < TOTAL_STEPS:
        chunk = min(SAFETY_INTERVAL_STEPS, TOTAL_STEPS - global_step)
        simulation.step(chunk)
        global_step += chunk
        force_write = global_step % REPORT_INTERVAL_STEPS == 0 or global_step == TOTAL_STEPS
        if sample(force_write):
            stopped = True
            print(f"Safety stop triggered at step {global_step}: {stop_reason}")
            break

    completed = (not stopped) and global_step == TOTAL_STEPS
    final_coords = positions_to_angstrom(last_positions)
    write_pdb(
        INPUT_PDB,
        final_coords,
        FINAL_PDB,
        [
            "REMARK 40chol 100 ps R2 restrained MD pilot final safe coordinates",
            "REMARK not production MD",
            f"REMARK completed_100ps {completed}",
            f"REMARK continued_from_checkpoint {continued_from_checkpoint}",
            f"REMARK last_step {global_step}",
            f"REMARK stop_reason {stop_reason or 'none'}",
        ],
    )
    simulation.saveState(str(STATE_XML))
    simulation.saveCheckpoint(str(CHECKPOINT_PATH))
    write_csv(MONITOR_CSV, rows, MONITOR_FIELDS)
    rmsd_stats = summarize_numeric(rows, "ligand_heavy_rmsd_vs_100ps_input_A")
    recommend_comparison = (
        completed
        and continued_from_checkpoint
        and not has_nan
        and to_int(last_row.get("ligand_near_chains_count"), 0) == 2
        and to_int(last_row.get("severe_lt_1A"), 999) == 0
        and to_int(last_row.get("close_lt_1p5A"), 999) == 0
        and to_float(last_row.get("ligand_min_distance_to_key_residues_A")) <= RETENTION_CUTOFF_A
        and rmsd_stats["max"] <= MAX_RMSD_RECOMMEND_A
    )
    result = {
        "completed_100ps": completed,
        "continued_from_checkpoint": continued_from_checkpoint,
        "fallback_used": fallback_used,
        "fallback_reason": fallback_reason,
        "last_step": global_step,
        "stop_reason": stop_reason,
        "has_nan": has_nan,
        "platform": platform.getName(),
        "final_row": last_row,
        "recommend_comparison": recommend_comparison,
    }
    write_summary(result, rows)
    print("")
    print("Terminal summary")
    print(f"completed 100 ps: {'YES' if completed else 'NO'}")
    print(f"continued from checkpoint: {'YES' if continued_from_checkpoint else 'NO'}")
    print(f"final temperature: {last_row.get('temperature_K', 'nan')}")
    print(f"final potential energy: {last_row.get('potential_energy_kJ_mol', 'nan')}")
    print(f"NaN: {'YES' if has_nan else 'NO'}")
    print(
        "final distance to TYR13 / VAL17 / SER20: "
        f"{last_row.get('ligand_distance_TYR13_A', 'nan')} / "
        f"{last_row.get('ligand_distance_VAL17_A', 'nan')} / "
        f"{last_row.get('ligand_distance_SER20_A', 'nan')}"
    )
    print(f"ligand final near chains count: {last_row.get('ligand_near_chains_count', 'nan')}")
    print(f"ligand final RMSD: {last_row.get('ligand_heavy_rmsd_vs_100ps_input_A', 'nan')}")
    print(f"max ligand RMSD: {fmt(rmsd_stats['max'])}")
    print(f"final severe/close contacts: {last_row.get('severe_lt_1A', 'nan')}/{last_row.get('close_lt_1p5A', 'nan')}")
    print(f"recommend 20chol vs 40chol comparison analysis: {'YES' if recommend_comparison else 'NO'}")
    print(f"summary report path: {SUMMARY_MD}")
    return result


def failure_result(error: str) -> dict[str, Any]:
    if not MONITOR_CSV.exists():
        write_csv(MONITOR_CSV, [], MONITOR_FIELDS)
    result = {
        "completed_100ps": False,
        "continued_from_checkpoint": False,
        "fallback_used": False,
        "fallback_reason": error,
        "last_step": 0,
        "stop_reason": error,
        "has_nan": False,
        "platform": "not reached",
        "final_row": {},
        "recommend_comparison": False,
    }
    write_summary(result, [])
    return result


def main() -> int:
    ensure_dirs()
    with LOG_PATH.open("w", encoding="utf-8") as log_file:
        tee_out = Tee(sys.__stdout__, log_file)
        tee_err = Tee(sys.__stderr__, log_file)
        with redirect_stdout(tee_out), redirect_stderr(tee_err):
            try:
                result = run_md100ps()
            except Exception as exc:
                error = repr(exc)
                print("FAILED: 40chol 100 ps R2 restrained MD failed:", error)
                traceback.print_exc()
                result = failure_result(error)
    return 0 if result.get("completed_100ps") and not result.get("has_nan") else 1


if __name__ == "__main__":
    raise SystemExit(main())
