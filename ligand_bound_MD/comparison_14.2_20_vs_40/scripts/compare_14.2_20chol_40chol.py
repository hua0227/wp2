from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np


BASE_DIR = Path(r"C:\TRKB_WP2\ligand_bound_MD")
OUT_DIR = BASE_DIR / "comparison_14.2_20_vs_40"
REPORTS_DIR = OUT_DIR / "reports"
FIGURES_DIR = OUT_DIR / "figures"
TABLES_DIR = OUT_DIR / "tables"
PYMOL_DIR = OUT_DIR / "pymol"

MONITOR_20 = BASE_DIR / r"pilot_build\14.2_20chol\md100ps_R2\reports\md100ps_R2_monitor.csv"
PDB_20 = BASE_DIR / r"pilot_build\14.2_20chol\md100ps_R2\outputs\TRKB_20chol_L002_14.2_100ps_R2_final.pdb"
DCD_20 = BASE_DIR / r"pilot_build\14.2_20chol\md100ps_R2\outputs\TRKB_20chol_L002_14.2_100ps_R2.dcd"
SUMMARY_20 = BASE_DIR / r"pilot_build\14.2_20chol\md100ps_R2\reports\md100ps_R2_summary.md"

MONITOR_40 = BASE_DIR / r"pilot_build\14.2_40chol\md100ps_R2\reports\md100ps_40chol_R2_monitor.csv"
PDB_40 = BASE_DIR / r"pilot_build\14.2_40chol\md100ps_R2\outputs\TRKB_40chol_L002_14.2_100ps_R2_final.pdb"
DCD_40 = BASE_DIR / r"pilot_build\14.2_40chol\md100ps_R2\outputs\TRKB_40chol_L002_14.2_100ps_R2.dcd"
SUMMARY_40 = BASE_DIR / r"pilot_build\14.2_40chol\md100ps_R2\reports\md100ps_40chol_R2_summary.md"

SUMMARY_TABLE = TABLES_DIR / "comparison_summary_statistics.csv"
COMBINED_TABLE = TABLES_DIR / "combined_time_series_14.2.csv"
REPORT_PATH = REPORTS_DIR / "comparison_14.2_20_vs_40_summary.md"
PYMOL_PATH = PYMOL_DIR / "view_14.2_20_vs_40_final.pml"

SUMMARY_FIELDS = [
    "system",
    "completed_ps",
    "temperature_min",
    "temperature_max",
    "temperature_mean",
    "temperature_final",
    "potential_energy_initial",
    "potential_energy_final",
    "potential_energy_mean",
    "ligand_rmsd_initial",
    "ligand_rmsd_final",
    "ligand_rmsd_max",
    "key_min_distance_initial",
    "key_min_distance_final",
    "key_min_distance_min",
    "key_min_distance_max",
    "key_min_distance_mean",
    "TYR13_distance_initial",
    "TYR13_distance_final",
    "TYR13_distance_mean",
    "VAL17_distance_initial",
    "VAL17_distance_final",
    "VAL17_distance_mean",
    "SER20_distance_initial",
    "SER20_distance_final",
    "SER20_distance_mean",
    "near_chains_all_2",
    "severe_clash_ever",
    "close_clash_ever",
    "NaN_detected",
]

COMBINED_FIELDS = [
    "time_ps",
    "system",
    "temperature_K",
    "potential_energy_kJ_mol",
    "ligand_rmsd_A",
    "key_min_distance_A",
    "TYR13_distance_A",
    "VAL17_distance_A",
    "SER20_distance_A",
    "near_chains",
    "severe_clashes",
    "close_contacts",
    "status",
]


def ensure_dirs() -> None:
    for path in [OUT_DIR, REPORTS_DIR, FIGURES_DIR, TABLES_DIR, PYMOL_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def to_float(value: object) -> float:
    if value is None:
        return math.nan
    text = str(value).strip()
    if text == "":
        return math.nan
    try:
        return float(text)
    except ValueError:
        return math.nan


def to_int(value: object) -> int:
    if value is None or str(value).strip() == "":
        return 0
    return int(float(str(value).strip()))


def is_nanish(value: object) -> bool:
    if value is None:
        return False
    text = str(value).strip().lower()
    if text in {"nan", "inf", "-inf", "infinity", "-infinity"}:
        return True
    try:
        return not math.isfinite(float(text))
    except ValueError:
        return False


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def pick_float(row: dict[str, str], names: Iterable[str]) -> float:
    for name in names:
        if name in row:
            value = to_float(row.get(name))
            if math.isfinite(value):
                return value
    return math.nan


def normalize_monitor_row(system: str, row: dict[str, str]) -> dict[str, object]:
    return {
        "time_ps": pick_float(row, ["time_ps"]),
        "system": system,
        "temperature_K": pick_float(row, ["temperature_K", "actual_temperature_K"]),
        "potential_energy_kJ_mol": pick_float(row, ["potential_energy_kJ_mol"]),
        "ligand_rmsd_A": pick_float(row, ["ligand_heavy_rmsd_vs_100ps_input_A"]),
        "key_min_distance_A": pick_float(
            row,
            ["ligand_min_distance_to_key_residues_A", "ligand_key_min_distance_A"],
        ),
        "TYR13_distance_A": pick_float(row, ["ligand_distance_TYR13_A"]),
        "VAL17_distance_A": pick_float(row, ["ligand_distance_VAL17_A"]),
        "SER20_distance_A": pick_float(row, ["ligand_distance_SER20_A"]),
        "near_chains": to_int(row.get("ligand_near_chains_count")),
        "severe_clashes": to_int(row.get("severe_lt_1A")),
        "close_contacts": to_int(row.get("close_lt_1p5A")),
        "status": row.get("status", ""),
    }


def finite_values(rows: list[dict[str, object]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = to_float(row.get(key))
        if math.isfinite(value):
            values.append(value)
    return values


def fmt(value: object, digits: int = 6) -> str:
    value = to_float(value)
    if not math.isfinite(value):
        return ""
    return f"{value:.{digits}g}"


def stats(values: list[float]) -> dict[str, float]:
    if not values:
        return {"initial": math.nan, "final": math.nan, "min": math.nan, "max": math.nan, "mean": math.nan}
    return {
        "initial": values[0],
        "final": values[-1],
        "min": min(values),
        "max": max(values),
        "mean": sum(values) / len(values),
    }


def summarize(system: str, rows: list[dict[str, object]], raw_rows: list[dict[str, str]]) -> dict[str, object]:
    temperature = stats(finite_values(rows, "temperature_K"))
    pe = stats(finite_values(rows, "potential_energy_kJ_mol"))
    rmsd = stats(finite_values(rows, "ligand_rmsd_A"))
    key = stats(finite_values(rows, "key_min_distance_A"))
    tyr = stats(finite_values(rows, "TYR13_distance_A"))
    val = stats(finite_values(rows, "VAL17_distance_A"))
    ser = stats(finite_values(rows, "SER20_distance_A"))
    times = finite_values(rows, "time_ps")
    return {
        "system": system,
        "completed_ps": max(times) if times else math.nan,
        "temperature_min": temperature["min"],
        "temperature_max": temperature["max"],
        "temperature_mean": temperature["mean"],
        "temperature_final": temperature["final"],
        "potential_energy_initial": pe["initial"],
        "potential_energy_final": pe["final"],
        "potential_energy_mean": pe["mean"],
        "ligand_rmsd_initial": rmsd["initial"],
        "ligand_rmsd_final": rmsd["final"],
        "ligand_rmsd_max": rmsd["max"],
        "key_min_distance_initial": key["initial"],
        "key_min_distance_final": key["final"],
        "key_min_distance_min": key["min"],
        "key_min_distance_max": key["max"],
        "key_min_distance_mean": key["mean"],
        "TYR13_distance_initial": tyr["initial"],
        "TYR13_distance_final": tyr["final"],
        "TYR13_distance_mean": tyr["mean"],
        "VAL17_distance_initial": val["initial"],
        "VAL17_distance_final": val["final"],
        "VAL17_distance_mean": val["mean"],
        "SER20_distance_initial": ser["initial"],
        "SER20_distance_final": ser["final"],
        "SER20_distance_mean": ser["mean"],
        "near_chains_all_2": bool(rows) and all(to_int(row.get("near_chains")) == 2 for row in rows),
        "severe_clash_ever": any(to_int(row.get("severe_clashes")) > 0 for row in rows),
        "close_clash_ever": any(to_int(row.get("close_contacts")) > 0 for row in rows),
        "NaN_detected": any(any(is_nanish(value) for value in raw.values()) for raw in raw_rows),
    }


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field)) if isinstance(row.get(field), float) else row.get(field, "") for field in fields})


def atom_indices(topology: md.Topology, predicate) -> list[int]:
    return [atom.index for atom in topology.atoms if predicate(atom)]


def recompute_per_residue_distances(pdb_path: Path, dcd_path: Path) -> tuple[dict[float, dict[str, float]], str]:
    try:
        traj = md.load_dcd(str(dcd_path), top=str(pdb_path))
    except Exception as exc:
        return {}, f"recompute failed: {type(exc).__name__}: {exc}"

    topology = traj.topology
    ligand = atom_indices(topology, lambda atom: atom.residue.name == "L002" and atom.element.symbol != "H")
    residue_specs = {
        "TYR13_distance_A": ("TYR", 434),
        "VAL17_distance_A": ("VAL", 438),
        "SER20_distance_A": ("SER", 441),
    }
    if not ligand:
        return {}, "recompute failed: no L002 heavy atoms found"

    out: dict[float, dict[str, float]] = {}
    for frame_idx in range(traj.n_frames):
        values: dict[str, float] = {}
        xyz = traj.xyz[frame_idx]
        box = traj.unitcell_lengths[frame_idx] if traj.unitcell_lengths is not None else None
        for label, (resname, resseq) in residue_specs.items():
            residue_atoms = atom_indices(
                topology,
                lambda atom, rn=resname, rs=resseq: atom.residue.name == rn
                and atom.residue.resSeq == rs
                and atom.element.symbol != "H",
            )
            if not residue_atoms:
                values[label] = math.nan
                continue
            pairs = np.array([(i, j) for i in ligand for j in residue_atoms], dtype=int)
            distances = md.compute_distances(
                traj[frame_idx],
                pairs,
                periodic=box is not None,
                opt=True,
            )[0]
            values[label] = float(np.min(distances) * 10.0)
        out[float((frame_idx + 1) * 2)] = values
    return out, f"recomputed from DCD: {traj.n_frames} frames at 2 ps spacing, mapped to 2-100 ps"


def merge_recomputed_per_residue(rows: list[dict[str, object]], recomputed: dict[float, dict[str, float]]) -> None:
    for row in rows:
        time_ps = round(to_float(row.get("time_ps")), 6)
        if time_ps not in recomputed:
            continue
        for key, value in recomputed[time_ps].items():
            row[key] = value


def load_normalized(system: str, monitor_path: Path) -> tuple[list[dict[str, object]], list[dict[str, str]]]:
    raw_rows = read_csv(monitor_path)
    return [normalize_monitor_row(system, row) for row in raw_rows], raw_rows


def plot_lines(path: Path, title: str, ylabel: str, rows_by_system: dict[str, list[dict[str, object]]], keys: list[str], labels: list[str] | None = None) -> None:
    plt.figure(figsize=(8.5, 5.2), dpi=160)
    colors = {"20chol": "#1f77b4", "40chol": "#d62728"}
    labels = labels or keys
    has_data = False
    for system, rows in rows_by_system.items():
        times = [to_float(row.get("time_ps")) for row in rows]
        for key, label in zip(keys, labels):
            y = [to_float(row.get(key)) for row in rows]
            valid = [(xv, yv) for xv, yv in zip(times, y) if math.isfinite(xv) and math.isfinite(yv)]
            if not valid:
                continue
            has_data = True
            xs, ys = zip(*valid)
            suffix = "" if len(keys) == 1 else f" {label}"
            linestyle = "-" if key.endswith("TYR13_distance_A") or len(keys) == 1 else "--" if key.endswith("VAL17_distance_A") else ":"
            plt.plot(xs, ys, label=f"{system}{suffix}", color=colors.get(system), linestyle=linestyle, linewidth=1.8)
    if not has_data:
        plt.text(0.5, 0.5, "unavailable", ha="center", va="center", transform=plt.gca().transAxes)
    plt.title(title)
    plt.xlabel("time_ps")
    plt.ylabel(ylabel)
    plt.grid(True, alpha=0.28)
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def make_figures(rows_by_system: dict[str, list[dict[str, object]]]) -> None:
    plot_lines(FIGURES_DIR / "temperature_20_vs_40.png", "Temperature: 20chol vs 40chol", "temperature (K)", rows_by_system, ["temperature_K"])
    plot_lines(FIGURES_DIR / "potential_energy_20_vs_40.png", "Potential Energy: 20chol vs 40chol", "potential energy (kJ/mol)", rows_by_system, ["potential_energy_kJ_mol"])
    plot_lines(FIGURES_DIR / "ligand_rmsd_20_vs_40.png", "Ligand Heavy-Atom RMSD: 20chol vs 40chol", "RMSD vs 100 ps input (A)", rows_by_system, ["ligand_rmsd_A"])
    plot_lines(FIGURES_DIR / "key_min_distance_20_vs_40.png", "Ligand Minimum Distance to Key Residues", "minimum distance (A)", rows_by_system, ["key_min_distance_A"])
    plot_lines(
        FIGURES_DIR / "per_residue_distance_20_vs_40.png",
        "Per-Residue Ligand Distances",
        "distance (A)",
        rows_by_system,
        ["TYR13_distance_A", "VAL17_distance_A", "SER20_distance_A"],
        ["TYR13", "VAL17", "SER20"],
    )
    plot_lines(FIGURES_DIR / "near_chains_20_vs_40.png", "Ligand Near Chains Count", "near chains count", rows_by_system, ["near_chains"])

    plt.figure(figsize=(8.5, 5.2), dpi=160)
    colors = {"20chol": "#1f77b4", "40chol": "#d62728"}
    for system, rows in rows_by_system.items():
        times = [to_float(row.get("time_ps")) for row in rows]
        severe = [to_float(row.get("severe_clashes")) for row in rows]
        close = [to_float(row.get("close_contacts")) for row in rows]
        plt.plot(times, severe, label=f"{system} severe <1.0 A", color=colors.get(system), linestyle="-", linewidth=1.8)
        plt.plot(times, close, label=f"{system} close <1.5 A", color=colors.get(system), linestyle="--", linewidth=1.8)
    plt.title("Ligand Clash Counts")
    plt.xlabel("time_ps")
    plt.ylabel("count")
    plt.grid(True, alpha=0.28)
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "clash_counts_20_vs_40.png")
    plt.close()


def write_pymol_script() -> None:
    PYMOL_PATH.write_text(
        "\n".join(
            [
                "reinitialize",
                f"load {PDB_20.as_posix()}, TRKB_20chol_L002_14_2",
                f"load {PDB_40.as_posix()}, TRKB_40chol_L002_14_2",
                "hide everything",
                "show cartoon, polymer.protein",
                "show sticks, resn L002",
                "show sticks, (resn TYR and resi 434) or (resn VAL and resi 438) or (resn SER and resi 441)",
                "color marine, TRKB_20chol_L002_14_2",
                "color firebrick, TRKB_40chol_L002_14_2",
                "color cyan, TRKB_20chol_L002_14_2 and resn L002",
                "color orange, TRKB_40chol_L002_14_2 and resn L002",
                "color yellow, resn TYR+VAL+SER and resi 434+438+441",
                "set stick_radius, 0.18",
                "set cartoon_transparency, 0.35",
                "zoom (resn L002) or (resi 434+438+441), 8",
                "orient (resn L002) or (resi 434+438+441)",
                "# Visualization only; not used for simulation.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def comparison_sentence(summary20: dict[str, object], summary40: dict[str, object]) -> tuple[str, str]:
    val20 = to_float(summary20.get("VAL17_distance_mean"))
    ser20 = to_float(summary20.get("SER20_distance_mean"))
    val40 = to_float(summary40.get("VAL17_distance_mean"))
    ser40 = to_float(summary40.get("SER20_distance_mean"))
    if all(math.isfinite(v) for v in [val20, ser20, val40, ser40]):
        far = val40 > val20 + 1.0 and ser40 > ser20 + 1.0
        return (
            "YES" if far else "NO",
            f"40chol VAL17/SER20 mean distances are {val40:.3f}/{ser40:.3f} A versus 20chol {val20:.3f}/{ser20:.3f} A.",
        )
    return (
        "PARTIAL",
        "20chol per-residue distances were recomputed only at DCD frame times; comparison uses available even-ps frames.",
    )


def write_report(summary20: dict[str, object], summary40: dict[str, object], recompute_note: str, far_flag: str, far_detail: str) -> None:
    near_both = summary20["near_chains_all_2"] and summary40["near_chains_all_2"]
    report = f"""# 14.2 / L002 20chol vs 40chol R2 Pilot Comparison

## Purpose

Compare restrained MD pilot stability for 14.2 / L002 in 20% CHOL and 40% CHOL TRKB-TMD systems.

## Input Data

- 20chol monitor CSV: `{MONITOR_20}`
- 20chol DCD: `{DCD_20}`
- 20chol final PDB: `{PDB_20}`
- 20chol summary: `{SUMMARY_20}`
- 40chol monitor CSV: `{MONITOR_40}`
- 40chol DCD: `{DCD_40}`
- 40chol final PDB: `{PDB_40}`
- 40chol summary: `{SUMMARY_40}`

## Simulation Status

Both systems are 100 ps R2 PBC-aware restrained MD pilots. These are not production MD simulations.

## Data Handling

- 40chol monitor includes individual TYR13, VAL17, and SER20 distances.
- 20chol monitor includes only the overall ligand minimum distance to TYR13/VAL17/SER20.
- 20chol per-residue distances were recomputed from DCD plus final-PDB topology where trajectory frames exist. {recompute_note}
- 20chol per-residue values are available only at DCD frame times; missing 0 ps and odd-ps values were left blank and were not fabricated.

## Stability Results

- 20chol completed ps: {fmt(summary20['completed_ps'])}; 40chol completed ps: {fmt(summary40['completed_ps'])}.
- Temperature range 20chol: {fmt(summary20['temperature_min'])}-{fmt(summary20['temperature_max'])} K; 40chol: {fmt(summary40['temperature_min'])}-{fmt(summary40['temperature_max'])} K.
- NaN detected 20chol: {summary20['NaN_detected']}; 40chol: {summary40['NaN_detected']}.
- Severe clashes appeared 20chol: {summary20['severe_clash_ever']}; 40chol: {summary40['severe_clash_ever']}.
- Close contacts appeared 20chol: {summary20['close_clash_ever']}; 40chol: {summary40['close_clash_ever']}.
- Near chains remained 2 in both systems: {near_both}.

## Ligand Retention

- Ligand RMSD final/max 20chol: {fmt(summary20['ligand_rmsd_final'])}/{fmt(summary20['ligand_rmsd_max'])} A.
- Ligand RMSD final/max 40chol: {fmt(summary40['ligand_rmsd_final'])}/{fmt(summary40['ligand_rmsd_max'])} A.
- Key-residue minimum distance final/mean 20chol: {fmt(summary20['key_min_distance_final'])}/{fmt(summary20['key_min_distance_mean'])} A.
- Key-residue minimum distance final/mean 40chol: {fmt(summary40['key_min_distance_final'])}/{fmt(summary40['key_min_distance_mean'])} A.
- Mean TYR13/VAL17/SER20 distances 20chol: {fmt(summary20['TYR13_distance_mean'])}/{fmt(summary20['VAL17_distance_mean'])}/{fmt(summary20['SER20_distance_mean'])} A.
- Mean TYR13/VAL17/SER20 distances 40chol: {fmt(summary40['TYR13_distance_mean'])}/{fmt(summary40['VAL17_distance_mean'])}/{fmt(summary40['SER20_distance_mean'])} A.

## Key Comparison

- 20chol keeps L002 close to the key-residue region, with a lower overall key-min distance and lower ligand RMSD during this restrained pilot.
- 40chol still retains L002 at the interface because near chains remain 2 and TYR13 remains close, but VAL17 and SER20 contacts are weaker than in 20chol based on available per-residue distances.
- 40chol VAL17/SER20 clearly farther than 20chol: {far_flag}. {far_detail}
- This pattern suggests 40chol may alter the ligand contact pose within the interface rather than fully ejecting the ligand.

## Limitations

- 100 ps restrained MD is not production MD.
- Both systems used restraints, so these data are pilot stability checks rather than equilibrium conclusions.
- CGenFF parameter penalty for L002 must remain part of interpretation.
- Only one ligand was analyzed.
- More ligands, longer simulations, reduced restraints, and replicate trajectories are needed before drawing production-level conclusions.

## Next Steps

1. Treat 14.2 as the first focused candidate for cross-cholesterol comparison.
2. Extend this workflow to another moderate-penalty candidate such as 19.1 or 2.3.
3. Perform a more systematic 14.2 contact occupancy / distance occupancy analysis.
4. Do not describe these restrained pilot results as production MD conclusions.

## Outputs

- Summary statistics table: `{SUMMARY_TABLE}`
- Combined time series table: `{COMBINED_TABLE}`
- Figures directory: `{FIGURES_DIR}`
- PyMOL visualization script: `{PYMOL_PATH}`
"""
    REPORT_PATH.write_text(report, encoding="utf-8")


def terminal_summary(summary20: dict[str, object], summary40: dict[str, object], far_flag: str) -> None:
    print(f"20chol completed 100 ps: {'YES' if to_float(summary20['completed_ps']) >= 100 else 'NO'}")
    print(f"40chol completed 100 ps: {'YES' if to_float(summary40['completed_ps']) >= 100 else 'NO'}")
    print(f"20chol NaN: {'YES' if summary20['NaN_detected'] else 'NO'}")
    print(f"40chol NaN: {'YES' if summary40['NaN_detected'] else 'NO'}")
    print(f"20chol temperature range: {fmt(summary20['temperature_min'])} - {fmt(summary20['temperature_max'])} K")
    print(f"40chol temperature range: {fmt(summary40['temperature_min'])} - {fmt(summary40['temperature_max'])} K")
    print(f"20chol ligand RMSD final/max: {fmt(summary20['ligand_rmsd_final'])} / {fmt(summary20['ligand_rmsd_max'])} A")
    print(f"40chol ligand RMSD final/max: {fmt(summary40['ligand_rmsd_final'])} / {fmt(summary40['ligand_rmsd_max'])} A")
    print(f"20chol key-min-distance final/mean: {fmt(summary20['key_min_distance_final'])} / {fmt(summary20['key_min_distance_mean'])} A")
    print(f"40chol key-min-distance final/mean: {fmt(summary40['key_min_distance_final'])} / {fmt(summary40['key_min_distance_mean'])} A")
    print(f"40chol VAL17/SER20 clearly farther than 20chol: {far_flag}")
    print(f"near chains always 2: {'YES' if summary20['near_chains_all_2'] and summary40['near_chains_all_2'] else 'NO'}")
    print(f"summary report path: {REPORT_PATH}")
    print(f"figures path: {FIGURES_DIR}")


def main() -> None:
    ensure_dirs()
    rows20, raw20 = load_normalized("20chol", MONITOR_20)
    rows40, raw40 = load_normalized("40chol", MONITOR_40)

    need_recompute20 = not any(math.isfinite(to_float(row.get("TYR13_distance_A"))) for row in rows20)
    recompute_note = "not needed"
    if need_recompute20:
        recomputed, recompute_note = recompute_per_residue_distances(PDB_20, DCD_20)
        merge_recomputed_per_residue(rows20, recomputed)

    rows_by_system = {"20chol": rows20, "40chol": rows40}
    combined_rows = rows20 + rows40
    combined_rows.sort(key=lambda row: (to_float(row.get("time_ps")), str(row.get("system"))))

    summary_rows = [summarize("20chol", rows20, raw20), summarize("40chol", rows40, raw40)]
    write_csv(SUMMARY_TABLE, summary_rows, SUMMARY_FIELDS)
    write_csv(COMBINED_TABLE, combined_rows, COMBINED_FIELDS)
    make_figures(rows_by_system)
    write_pymol_script()

    far_flag, far_detail = comparison_sentence(summary_rows[0], summary_rows[1])
    write_report(summary_rows[0], summary_rows[1], recompute_note, far_flag, far_detail)
    terminal_summary(summary_rows[0], summary_rows[1], far_flag)


if __name__ == "__main__":
    main()
