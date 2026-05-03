from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


BASE_DIR = Path(r"C:\TRKB_WP2\ligand_bound_MD")
PILOT_BUILD_DIR = BASE_DIR / "pilot_build"
BATCH_DIR = PILOT_BUILD_DIR / "batch_2.3_6.2_17.2"
BATCH_REPORTS_DIR = BATCH_DIR / "reports"
BATCH_STATUS_CSV = BATCH_REPORTS_DIR / "batch_gate_status.csv"
BATCH_SUMMARY_CSV = BATCH_REPORTS_DIR / "batch_comparison_occupancy_summary.csv"
CANDIDATES = [("2.3", "L006"), ("6.2", "L008"), ("17.2", "L010")]
SYSTEMS = ["20chol", "40chol"]
RESIDUES = ["TYR13", "VAL17", "SER20"]
RESIDUE_THRESHOLDS = [4.0, 5.0, 6.0, 8.0]
KEY_THRESHOLDS = [4.0, 5.0, 6.0]
COLORS = {"20chol": "#1f77b4", "40chol": "#d62728"}

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
    "TYR13_distance_min",
    "TYR13_distance_max",
    "TYR13_distance_mean",
    "VAL17_distance_initial",
    "VAL17_distance_final",
    "VAL17_distance_min",
    "VAL17_distance_max",
    "VAL17_distance_mean",
    "SER20_distance_initial",
    "SER20_distance_final",
    "SER20_distance_min",
    "SER20_distance_max",
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

OCCUPANCY_FIELDS = [
    "system",
    "per_residue_frame_count",
    "TYR13_occ_lt4A",
    "TYR13_occ_lt5A",
    "TYR13_occ_lt6A",
    "TYR13_occ_lt8A",
    "VAL17_occ_lt4A",
    "VAL17_occ_lt5A",
    "VAL17_occ_lt6A",
    "VAL17_occ_lt8A",
    "SER20_occ_lt4A",
    "SER20_occ_lt5A",
    "SER20_occ_lt6A",
    "SER20_occ_lt8A",
    "key_min_occ_lt4A",
    "key_min_occ_lt5A",
    "key_min_occ_lt6A",
    "near_chains_eq2_fraction",
    "ligand_rmsd_mean",
    "ligand_rmsd_median",
    "ligand_rmsd_max",
    "ligand_rmsd_p95",
]

BATCH_FIELDS = [
    "ligand_id",
    "resname",
    "system",
    "completed_100ps",
    "NaN_detected",
    "near_chains_all_2",
    "temperature_mean",
    "ligand_rmsd_mean",
    "ligand_rmsd_final",
    "ligand_rmsd_max",
    "TYR13_occ_lt5A",
    "VAL17_occ_lt5A",
    "SER20_occ_lt5A",
    "TYR13_mean_distance",
    "VAL17_mean_distance",
    "SER20_mean_distance",
    "key_min_distance_mean",
    "key_min_occ_lt5A",
    "severe_clash_ever",
    "close_clash_ever",
    "contact_mode_label",
    "cholesterol_response_label",
    "candidate_interpretation",
]


def comparison_dir(ligand_id: str) -> Path:
    return BASE_DIR / f"comparison_{ligand_id}_20_vs_40"


def paths_for_ligand(ligand_id: str) -> dict[str, Path]:
    out_dir = comparison_dir(ligand_id)
    return {
        "out_dir": out_dir,
        "scripts": out_dir / "scripts",
        "reports": out_dir / "reports",
        "tables": out_dir / "tables",
        "figures": out_dir / "figures",
        "pymol": out_dir / "pymol",
        "summary_table": out_dir / "tables" / "comparison_summary_statistics.csv",
        "combined_table": out_dir / "tables" / f"combined_time_series_{ligand_id}.csv",
        "occupancy_table": out_dir / "tables" / "contact_occupancy_summary.csv",
        "report_path": out_dir / "reports" / f"comparison_{ligand_id}_20_vs_40_summary.md",
        "pymol_path": out_dir / "pymol" / f"view_{ligand_id}_20_vs_40_final.pml",
    }


def ensure_dir_tree(ligand_id: str) -> dict[str, Path]:
    paths = paths_for_ligand(ligand_id)
    for key in ["out_dir", "scripts", "reports", "tables", "figures", "pymol"]:
        paths[key].mkdir(parents=True, exist_ok=True)
    return paths


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


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


def finite(values: Iterable[float]) -> list[float]:
    return [value for value in values if math.isfinite(value)]


def mean(values: Iterable[float]) -> float:
    valid = finite(values)
    if not valid:
        return math.nan
    return float(np.mean(np.asarray(valid, dtype=float)))


def median(values: Iterable[float]) -> float:
    valid = finite(values)
    if not valid:
        return math.nan
    return float(np.median(np.asarray(valid, dtype=float)))


def percentile(values: Iterable[float], q: float) -> float:
    valid = finite(values)
    if not valid:
        return math.nan
    return float(np.percentile(np.asarray(valid, dtype=float), q))


def max_value(values: Iterable[float]) -> float:
    valid = finite(values)
    return max(valid) if valid else math.nan


def stats(values: list[float]) -> dict[str, float]:
    valid = finite(values)
    if not valid:
        return {"initial": math.nan, "final": math.nan, "min": math.nan, "max": math.nan, "mean": math.nan}
    return {
        "initial": valid[0],
        "final": valid[-1],
        "min": min(valid),
        "max": max(valid),
        "mean": sum(valid) / len(valid),
    }


def occupancy(values: list[float], threshold: float) -> dict[str, float]:
    valid = finite(values)
    hits = sum(1 for value in valid if value < threshold)
    total = len(valid)
    return {"hits": hits, "total": total, "fraction": hits / total if total else math.nan}


def fmt_text(value: object, digits: int = 3) -> str:
    number = to_float(value)
    if math.isfinite(number):
        return f"{number:.{digits}f}"
    if isinstance(value, bool):
        return "YES" if value else "NO"
    return str(value)


def pct_text(value: object) -> str:
    number = to_float(value)
    if math.isfinite(number):
        return f"{100.0 * number:.1f}%"
    return "NA"


def bool_text(value: object) -> str:
    return "YES" if bool(value) else "NO"


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
        "total_energy_kJ_mol": pick_float(row, ["total_energy_kJ_mol"]),
        "ligand_rmsd_A": pick_float(row, ["ligand_heavy_rmsd_vs_100ps_input_A"]),
        "key_min_distance_A": pick_float(row, ["ligand_min_distance_to_key_residues_A", "ligand_key_min_distance_A"]),
        "TYR13_distance_A": pick_float(row, ["ligand_distance_TYR13_A"]),
        "VAL17_distance_A": pick_float(row, ["ligand_distance_VAL17_A"]),
        "SER20_distance_A": pick_float(row, ["ligand_distance_SER20_A"]),
        "near_chains": to_int(row.get("ligand_near_chains_count")),
        "severe_clashes": to_int(row.get("severe_lt_1A")),
        "close_contacts": to_int(row.get("close_lt_1p5A")),
        "status": row.get("status", ""),
    }


def finite_values(rows: list[dict[str, object]], key: str) -> list[float]:
    return finite(to_float(row.get(key)) for row in rows)


def summarize_comparison(system: str, rows: list[dict[str, object]], raw_rows: list[dict[str, str]]) -> dict[str, object]:
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
        "TYR13_distance_min": tyr["min"],
        "TYR13_distance_max": tyr["max"],
        "TYR13_distance_mean": tyr["mean"],
        "VAL17_distance_initial": val["initial"],
        "VAL17_distance_final": val["final"],
        "VAL17_distance_min": val["min"],
        "VAL17_distance_max": val["max"],
        "VAL17_distance_mean": val["mean"],
        "SER20_distance_initial": ser["initial"],
        "SER20_distance_final": ser["final"],
        "SER20_distance_min": ser["min"],
        "SER20_distance_max": ser["max"],
        "SER20_distance_mean": ser["mean"],
        "near_chains_all_2": all(to_int(row.get("near_chains")) == 2 for row in rows),
        "severe_clash_ever": any(to_int(row.get("severe_clashes")) > 0 for row in rows),
        "close_clash_ever": any(to_int(row.get("close_contacts")) > 0 for row in rows),
        "NaN_detected": any(any(is_nanish(value) for value in row.values()) for row in raw_rows),
    }


def summarize_occupancy_system(system: str, rows: list[dict[str, object]]) -> dict[str, object]:
    summary: dict[str, object] = {"system": system}
    for residue in RESIDUES:
        values = finite_values(rows, f"{residue}_distance_A")
        if residue == RESIDUES[0]:
            summary["per_residue_frame_count"] = len(values)
        for threshold in RESIDUE_THRESHOLDS:
            result = occupancy(values, threshold)
            summary[f"{residue}_occ_lt{int(threshold)}A"] = result["fraction"]

    key_values = finite_values(rows, "key_min_distance_A")
    for threshold in KEY_THRESHOLDS:
        result = occupancy(key_values, threshold)
        summary[f"key_min_occ_lt{int(threshold)}A"] = result["fraction"]

    near_values = [to_int(row.get("near_chains")) for row in rows]
    summary["near_chains_eq2_fraction"] = (
        sum(1 for value in near_values if value == 2) / len(near_values) if near_values else math.nan
    )

    rmsd_values = finite_values(rows, "ligand_rmsd_A")
    summary["ligand_rmsd_mean"] = mean(rmsd_values)
    summary["ligand_rmsd_median"] = median(rmsd_values)
    summary["ligand_rmsd_max"] = max_value(rmsd_values)
    summary["ligand_rmsd_p95"] = percentile(rmsd_values, 95)
    return summary


def make_line_plot(
    rows_by_system: dict[str, list[dict[str, object]]],
    figure_path: Path,
    title: str,
    ylabel: str,
    ykey: str,
) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 5.2), dpi=160)
    for system in SYSTEMS:
        rows = rows_by_system[system]
        ax.plot(
            [to_float(row.get("time_ps")) for row in rows],
            [to_float(row.get(ykey)) for row in rows],
            label=system,
            color=COLORS[system],
            linewidth=1.8,
        )
    ax.set_xlabel("time (ps)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.28)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figure_path)
    plt.close(fig)


def make_clash_plot(rows_by_system: dict[str, list[dict[str, object]]], figure_path: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8.8, 6.6), dpi=160, sharex=True)
    for system in SYSTEMS:
        rows = rows_by_system[system]
        times = [to_float(row.get("time_ps")) for row in rows]
        axes[0].plot(times, [to_int(row.get("severe_clashes")) for row in rows], label=system, color=COLORS[system], linewidth=1.8)
        axes[1].plot(times, [to_int(row.get("close_contacts")) for row in rows], label=system, color=COLORS[system], linewidth=1.8)
    axes[0].set_ylabel("severe clashes\n(<1.0 A)")
    axes[1].set_ylabel("close contacts\n(<1.5 A)")
    axes[1].set_xlabel("time (ps)")
    axes[0].set_title("Clash/contact counts vs time")
    for ax in axes:
        ax.grid(True, alpha=0.28)
        ax.legend()
    fig.tight_layout()
    fig.savefig(figure_path)
    plt.close(fig)


def make_per_residue_distance_plot(rows_by_system: dict[str, list[dict[str, object]]], figure_path: Path) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(8.8, 9.2), dpi=160, sharex=True)
    for index, residue in enumerate(RESIDUES):
        ax = axes[index]
        for system in SYSTEMS:
            rows = rows_by_system[system]
            ax.plot(
                [to_float(row.get("time_ps")) for row in rows],
                [to_float(row.get(f"{residue}_distance_A")) for row in rows],
                label=system,
                color=COLORS[system],
                linewidth=1.8,
            )
        ax.set_ylabel(f"{residue} (A)")
        ax.grid(True, alpha=0.28)
        ax.legend()
    axes[0].set_title("Per-residue ligand distances vs time")
    axes[-1].set_xlabel("time (ps)")
    fig.tight_layout()
    fig.savefig(figure_path)
    plt.close(fig)


def make_occupancy_barplot(occ_rows: dict[str, dict[str, object]], figure_path: Path) -> None:
    labels = [
        "TYR13<5",
        "VAL17<5",
        "SER20<5",
        "key-min<5",
        "nearChains=2",
    ]
    series = {
        "20chol": [
            to_float(occ_rows["20chol"]["TYR13_occ_lt5A"]),
            to_float(occ_rows["20chol"]["VAL17_occ_lt5A"]),
            to_float(occ_rows["20chol"]["SER20_occ_lt5A"]),
            to_float(occ_rows["20chol"]["key_min_occ_lt5A"]),
            to_float(occ_rows["20chol"]["near_chains_eq2_fraction"]),
        ],
        "40chol": [
            to_float(occ_rows["40chol"]["TYR13_occ_lt5A"]),
            to_float(occ_rows["40chol"]["VAL17_occ_lt5A"]),
            to_float(occ_rows["40chol"]["SER20_occ_lt5A"]),
            to_float(occ_rows["40chol"]["key_min_occ_lt5A"]),
            to_float(occ_rows["40chol"]["near_chains_eq2_fraction"]),
        ],
    }
    x = np.arange(len(labels))
    width = 0.35
    fig, ax = plt.subplots(figsize=(9.3, 5.2), dpi=160)
    ax.bar(x - width / 2, series["20chol"], width=width, label="20chol", color=COLORS["20chol"])
    ax.bar(x + width / 2, series["40chol"], width=width, label="40chol", color=COLORS["40chol"])
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("occupancy fraction")
    ax.set_title("Contact occupancy comparison")
    ax.grid(True, axis="y", alpha=0.28)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figure_path)
    plt.close(fig)


def make_per_residue_boxplot(rows_by_system: dict[str, list[dict[str, object]]], figure_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 5.4), dpi=160)
    positions = []
    data = []
    labels = []
    colors = []
    for index, residue in enumerate(RESIDUES):
        base = 1 + index * 3
        for offset, system in enumerate(SYSTEMS):
            values = finite(to_float(row.get(f"{residue}_distance_A")) for row in rows_by_system[system])
            positions.append(base + offset * 0.9)
            data.append(values)
            labels.append(f"{residue}\n{system}")
            colors.append(COLORS[system])
    bp = ax.boxplot(data, positions=positions, patch_artist=True, widths=0.7, showfliers=False)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.65)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    ax.set_ylabel("distance (A)")
    ax.set_title("Per-residue distance distributions")
    ax.grid(True, axis="y", alpha=0.28)
    fig.tight_layout()
    fig.savefig(figure_path)
    plt.close(fig)


def make_rmsd_distribution_plot(rows_by_system: dict[str, list[dict[str, object]]], figure_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 5.2), dpi=160)
    upper = max(
        max_value(to_float(row.get("ligand_rmsd_A")) for row in rows_by_system[system])
        for system in SYSTEMS
    ) + 0.25
    bins = np.linspace(0.0, upper, 14)
    for system in SYSTEMS:
        values = finite(to_float(row.get("ligand_rmsd_A")) for row in rows_by_system[system])
        ax.hist(values, bins=bins, alpha=0.45, label=system, color=COLORS[system], density=True, edgecolor="black", linewidth=0.5)
        ax.axvline(mean(values), color=COLORS[system], linestyle="--", linewidth=1.5)
    ax.set_xlabel("ligand RMSD (A)")
    ax.set_ylabel("density")
    ax.set_title("Ligand RMSD distribution")
    ax.grid(True, axis="y", alpha=0.28)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figure_path)
    plt.close(fig)


def classify_contact_mode(occ_row: dict[str, object]) -> str:
    tyr = to_float(occ_row.get("TYR13_occ_lt5A"))
    val = to_float(occ_row.get("VAL17_occ_lt5A"))
    ser = to_float(occ_row.get("SER20_occ_lt5A"))
    if all(math.isfinite(v) and v >= 0.8 for v in [tyr, val, ser]):
        return "balanced_three_residue_contact"
    if math.isfinite(val) and math.isfinite(ser) and math.isfinite(tyr) and val >= 0.8 and ser >= 0.8 and tyr < 0.5:
        return "VAL17_SER20_dominant"
    if math.isfinite(tyr) and tyr >= 0.8 and ((math.isfinite(val) and val <= 0.2) or (math.isfinite(ser) and ser <= 0.2)):
        return "TYR13_shifted"
    return "weak_or_mixed_contact"


def classify_cholesterol_response(rows_by_system: dict[str, dict[str, object]]) -> str:
    row20 = rows_by_system["20chol"]
    row40 = rows_by_system["40chol"]
    mode20 = str(row20.get("contact_mode_label", ""))
    mode40 = str(row40.get("contact_mode_label", ""))
    if mode20 == "balanced_three_residue_contact" and mode40 == "balanced_three_residue_contact":
        return "cholesterol_robust_balanced"
    if bool(row20.get("near_chains_all_2")) and bool(row40.get("near_chains_all_2")):
        deltas = []
        for key in ["TYR13_occ_lt5A", "VAL17_occ_lt5A", "SER20_occ_lt5A"]:
            v20 = to_float(row20.get(key))
            v40 = to_float(row40.get(key))
            if math.isfinite(v20) and math.isfinite(v40):
                deltas.append(v20 - v40)
        if deltas and max(deltas) >= 0.25:
            return "interface_retained_but_weaker_contact"
    if mode20 != mode40 and mode20 and mode40:
        return "cholesterol_sensitive_contact_shift"
    return "weak_or_inconclusive"


def candidate_interpretation(ligand_id: str, resname: str, rows_by_system: dict[str, dict[str, object]], response_label: str) -> str:
    mode20 = rows_by_system["20chol"]["contact_mode_label"]
    mode40 = rows_by_system["40chol"]["contact_mode_label"]
    if response_label == "cholesterol_robust_balanced":
        return f"{ligand_id} / {resname} keeps strong three-residue contact in both cholesterol conditions and looks like the most balanced pilot candidate."
    if response_label == "cholesterol_sensitive_contact_shift":
        return f"{ligand_id} / {resname} retains the interface but shifts contact mode between 20chol and 40chol, suggesting cholesterol-sensitive pose remodeling."
    if response_label == "interface_retained_but_weaker_contact":
        return f"{ligand_id} / {resname} keeps interface retention in both systems, but one or more key-residue occupancies soften in 40chol."
    return f"{ligand_id} / {resname} shows mixed or weaker restrained-pilot contact behavior ({mode20} vs {mode40}); treat it as inconclusive without longer or less restrained sampling."


def infer_contact_style(contact_mode_label: str, near_chains_all_2: bool) -> str:
    if contact_mode_label == "balanced_three_residue_contact":
        return "balanced contact"
    if contact_mode_label in {"VAL17_SER20_dominant", "TYR13_shifted"}:
        return "contact shift"
    if near_chains_all_2:
        return "interface retention only"
    return "weak contact"


def discover_inputs_for_ligand(ligand_id: str, resname: str) -> dict[str, dict[str, Path]]:
    discovered: dict[str, dict[str, Path]] = {}
    for system in SYSTEMS:
        md_dir = PILOT_BUILD_DIR / f"{ligand_id}_{system}" / "md100ps_R2"
        direct = {
            "md_dir": md_dir,
            "monitor_csv": md_dir / "reports" / f"md100ps_{system}_R2_monitor.csv",
            "result_json": md_dir / "reports" / f"md100ps_{system}_R2_result.json",
            "summary_md": md_dir / "reports" / f"md100ps_{system}_R2_summary.md",
            "final_pdb": md_dir / "outputs" / f"TRKB_{system}_{resname}_{ligand_id}_100ps_R2_final.pdb",
            "dcd": md_dir / "outputs" / f"TRKB_{system}_{resname}_{ligand_id}_100ps_R2.dcd",
        }
        if all(path.exists() for key, path in direct.items() if key != "md_dir"):
            discovered[system] = direct
            continue
        discovered[system] = direct
    return discovered


def make_pymol_script(ligand_id: str, resname: str, input_paths: dict[str, dict[str, Path]], pymol_path: Path) -> None:
    pdb20 = input_paths["20chol"]["final_pdb"].as_posix()
    pdb40 = input_paths["40chol"]["final_pdb"].as_posix()
    obj20 = f"TRKB_20chol_{resname}_{ligand_id.replace('.', '_')}"
    obj40 = f"TRKB_40chol_{resname}_{ligand_id.replace('.', '_')}"
    pymol_path.write_text(
        "\n".join(
            [
                "reinitialize",
                f"load {pdb20}, {obj20}",
                f"load {pdb40}, {obj40}",
                "hide everything",
                "show cartoon, polymer.protein",
                f"show sticks, resn {resname}",
                "show sticks, (resn TYR and resi 434) or (resn VAL and resi 438) or (resn SER and resi 441)",
                f"color marine, {obj20}",
                f"color firebrick, {obj40}",
                f"color cyan, {obj20} and resn {resname}",
                f"color orange, {obj40} and resn {resname}",
                "color yellow, resn TYR+VAL+SER and resi 434+438+441",
                "set stick_radius, 0.18",
                "set cartoon_transparency, 0.35",
                f"zoom (resn {resname}) or (resi 434+438+441), 8",
                f"orient (resn {resname}) or (resi 434+438+441)",
                "# Visualization only; not used for simulation.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def write_wrapper_script(ligand_id: str, resname: str, scripts_dir: Path) -> Path:
    wrapper_path = scripts_dir / f"compare_{ligand_id}_20chol_40chol.py"
    wrapper_path.write_text(
        "\n".join(
            [
                "from __future__ import annotations",
                "",
                "import importlib.util",
                "from pathlib import Path",
                "",
                "MODULE_PATH = Path(r\"C:\\TRKB_WP2\\ligand_bound_MD\\pilot_build\\batch_2.3_6.2_17.2\\scripts\\batch_compare_candidates.py\")",
                "SPEC = importlib.util.spec_from_file_location('batch_compare_candidates', MODULE_PATH)",
                "MODULE = importlib.util.module_from_spec(SPEC)",
                "assert SPEC.loader is not None",
                "SPEC.loader.exec_module(MODULE)",
                "",
                "if __name__ == '__main__':",
                f"    MODULE.run_single_ligand('{ligand_id}', '{resname}')",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return wrapper_path


def write_report(
    ligand_id: str,
    resname: str,
    input_paths: dict[str, dict[str, Path]],
    paths: dict[str, Path],
    summary20: dict[str, object],
    summary40: dict[str, object],
    occ20: dict[str, object],
    occ40: dict[str, object],
    contact_style20: str,
    contact_style40: str,
    contact_mode20: str,
    contact_mode40: str,
    response_label: str,
    interpretation: str,
) -> None:
    report = f"""# {ligand_id} / {resname} 20chol vs 40chol R2 Pilot Comparison

## 分析目的

比较 {ligand_id} / {resname} 在 20% CHOL 与 40% CHOL TRKB-TMD 体系中的 restrained MD pilot 稳定性。

## 输入数据

- 20chol monitor CSV: `{input_paths['20chol']['monitor_csv']}`
- 20chol DCD: `{input_paths['20chol']['dcd']}`
- 20chol final PDB: `{input_paths['20chol']['final_pdb']}`
- 20chol summary: `{input_paths['20chol']['summary_md']}`
- 40chol monitor CSV: `{input_paths['40chol']['monitor_csv']}`
- 40chol DCD: `{input_paths['40chol']['dcd']}`
- 40chol final PDB: `{input_paths['40chol']['final_pdb']}`
- 40chol summary: `{input_paths['40chol']['summary_md']}`

## 模拟状态说明

两个体系都属于 100 ps R2 PBC-aware restrained MD pilot。它们不是 production MD 模拟，本报告也不是 production MD 结论。

## 稳定性结果

- 20chol completed ps: {fmt_text(summary20['completed_ps'])}; 40chol completed ps: {fmt_text(summary40['completed_ps'])}
- 20chol temperature range: {fmt_text(summary20['temperature_min'])}-{fmt_text(summary20['temperature_max'])} K
- 40chol temperature range: {fmt_text(summary40['temperature_min'])}-{fmt_text(summary40['temperature_max'])} K
- NaN detected, 20chol/40chol: {bool_text(summary20['NaN_detected'])} / {bool_text(summary40['NaN_detected'])}
- severe clash ever, 20chol/40chol: {bool_text(summary20['severe_clash_ever'])} / {bool_text(summary40['severe_clash_ever'])}
- close clash ever, 20chol/40chol: {bool_text(summary20['close_clash_ever'])} / {bool_text(summary40['close_clash_ever'])}
- near chains remained 2 throughout, 20chol/40chol: {bool_text(summary20['near_chains_all_2'])} / {bool_text(summary40['near_chains_all_2'])}

## Ligand 保留情况

- 20chol ligand RMSD final/max: {fmt_text(summary20['ligand_rmsd_final'])}/{fmt_text(summary20['ligand_rmsd_max'])} A
- 40chol ligand RMSD final/max: {fmt_text(summary40['ligand_rmsd_final'])}/{fmt_text(summary40['ligand_rmsd_max'])} A
- 20chol TYR13 / VAL17 / SER20 final: {fmt_text(summary20['TYR13_distance_final'])} / {fmt_text(summary20['VAL17_distance_final'])} / {fmt_text(summary20['SER20_distance_final'])} A
- 40chol TYR13 / VAL17 / SER20 final: {fmt_text(summary40['TYR13_distance_final'])} / {fmt_text(summary40['VAL17_distance_final'])} / {fmt_text(summary40['SER20_distance_final'])} A
- 20chol TYR13 / VAL17 / SER20 mean: {fmt_text(summary20['TYR13_distance_mean'])} / {fmt_text(summary20['VAL17_distance_mean'])} / {fmt_text(summary20['SER20_distance_mean'])} A
- 40chol TYR13 / VAL17 / SER20 mean: {fmt_text(summary40['TYR13_distance_mean'])} / {fmt_text(summary40['VAL17_distance_mean'])} / {fmt_text(summary40['SER20_distance_mean'])} A
- 20chol key-min distance final/mean: {fmt_text(summary20['key_min_distance_final'])}/{fmt_text(summary20['key_min_distance_mean'])} A
- 40chol key-min distance final/mean: {fmt_text(summary40['key_min_distance_final'])}/{fmt_text(summary40['key_min_distance_mean'])} A

## Contact Occupancy 结果

- TYR13 `<4 / <5 / <6 / <8 A` occupancy
  - 20chol: {pct_text(occ20['TYR13_occ_lt4A'])} / {pct_text(occ20['TYR13_occ_lt5A'])} / {pct_text(occ20['TYR13_occ_lt6A'])} / {pct_text(occ20['TYR13_occ_lt8A'])}
  - 40chol: {pct_text(occ40['TYR13_occ_lt4A'])} / {pct_text(occ40['TYR13_occ_lt5A'])} / {pct_text(occ40['TYR13_occ_lt6A'])} / {pct_text(occ40['TYR13_occ_lt8A'])}
- VAL17 `<4 / <5 / <6 / <8 A` occupancy
  - 20chol: {pct_text(occ20['VAL17_occ_lt4A'])} / {pct_text(occ20['VAL17_occ_lt5A'])} / {pct_text(occ20['VAL17_occ_lt6A'])} / {pct_text(occ20['VAL17_occ_lt8A'])}
  - 40chol: {pct_text(occ40['VAL17_occ_lt4A'])} / {pct_text(occ40['VAL17_occ_lt5A'])} / {pct_text(occ40['VAL17_occ_lt6A'])} / {pct_text(occ40['VAL17_occ_lt8A'])}
- SER20 `<4 / <5 / <6 / <8 A` occupancy
  - 20chol: {pct_text(occ20['SER20_occ_lt4A'])} / {pct_text(occ20['SER20_occ_lt5A'])} / {pct_text(occ20['SER20_occ_lt6A'])} / {pct_text(occ20['SER20_occ_lt8A'])}
  - 40chol: {pct_text(occ40['SER20_occ_lt4A'])} / {pct_text(occ40['SER20_occ_lt5A'])} / {pct_text(occ40['SER20_occ_lt6A'])} / {pct_text(occ40['SER20_occ_lt8A'])}
- key-min-distance `<4 / <5 / <6 A` occupancy
  - 20chol: {pct_text(occ20['key_min_occ_lt4A'])} / {pct_text(occ20['key_min_occ_lt5A'])} / {pct_text(occ20['key_min_occ_lt6A'])}
  - 40chol: {pct_text(occ40['key_min_occ_lt4A'])} / {pct_text(occ40['key_min_occ_lt5A'])} / {pct_text(occ40['key_min_occ_lt6A'])}
- near chains = 2 fraction
  - 20chol: {pct_text(occ20['near_chains_eq2_fraction'])}
  - 40chol: {pct_text(occ40['near_chains_eq2_fraction'])}
- ligand RMSD mean / median / max / 95th percentile
  - 20chol: {fmt_text(occ20['ligand_rmsd_mean'])} / {fmt_text(occ20['ligand_rmsd_median'])} / {fmt_text(occ20['ligand_rmsd_max'])} / {fmt_text(occ20['ligand_rmsd_p95'])} A
  - 40chol: {fmt_text(occ40['ligand_rmsd_mean'])} / {fmt_text(occ40['ligand_rmsd_median'])} / {fmt_text(occ40['ligand_rmsd_max'])} / {fmt_text(occ40['ligand_rmsd_p95'])} A

## 20chol vs 40chol 接触模式差异

- 20chol contact mode label: `{contact_mode20}` ({contact_style20})
- 40chol contact mode label: `{contact_mode40}` ({contact_style40})
- Cholesterol response label: `{response_label}`
- Candidate interpretation: {interpretation}

## 局限性

- 100 ps restrained MD 不是 production MD。
- 两个体系都使用了 restraints，因此这更适合作为 pilot stability / retention 检查，而不是平衡态结论。
- CGenFF 参数 penalty 仍然需要和本分析一起解读。
- 目前仍只覆盖单个 ligand、单条 trajectory。
- 后续仍需要更多 ligand、更长时间、更少 restraints，以及 replicate trajectories 支持。

## 下一步建议

1. 如果 {ligand_id} 双体系 occupancy 保持良好，可作为后续跨胆固醇横向比较的重点候选。
2. 后续可进入 14.2 vs 19.1 vs {ligand_id} 的横向比较。
3. 暂时不要称为 production MD 结论。

## 输出

- Summary statistics table: `{paths['summary_table']}`
- Combined time series table: `{paths['combined_table']}`
- Contact occupancy summary table: `{paths['occupancy_table']}`
- Figures directory: `{paths['figures']}`
- PyMOL visualization script: `{paths['pymol_path']}`
"""
    paths["report_path"].write_text(report, encoding="utf-8")


def analyze_single_ligand(ligand_id: str, resname: str) -> dict[str, object]:
    paths = ensure_dir_tree(ligand_id)
    input_paths = discover_inputs_for_ligand(ligand_id, resname)

    missing = []
    for system in SYSTEMS:
        for key in ["monitor_csv", "summary_md", "final_pdb", "dcd", "result_json"]:
            if not input_paths[system][key].exists():
                missing.append(f"{system}:{key}")
    if missing:
        raise FileNotFoundError(f"Missing required inputs for {ligand_id}/{resname}: {', '.join(missing)}")

    raw_rows_by_system: dict[str, list[dict[str, str]]] = {}
    rows_by_system: dict[str, list[dict[str, object]]] = {}
    summaries: dict[str, dict[str, object]] = {}
    occupancies: dict[str, dict[str, object]] = {}
    result_json_by_system: dict[str, dict[str, object]] = {}
    for system in SYSTEMS:
        raw_rows = read_csv(input_paths[system]["monitor_csv"])
        rows = [normalize_monitor_row(system, row) for row in raw_rows]
        raw_rows_by_system[system] = raw_rows
        rows_by_system[system] = rows
        summaries[system] = summarize_comparison(system, rows, raw_rows)
        occupancies[system] = summarize_occupancy_system(system, rows)
        result_json_by_system[system] = read_json(input_paths[system]["result_json"])

    combined_rows = rows_by_system["20chol"] + rows_by_system["40chol"]
    combined_rows.sort(key=lambda row: (to_float(row.get("time_ps")), str(row.get("system"))))

    write_csv(paths["summary_table"], [summaries["20chol"], summaries["40chol"]], SUMMARY_FIELDS)
    write_csv(paths["combined_table"], combined_rows, COMBINED_FIELDS)
    write_csv(paths["occupancy_table"], [occupancies["20chol"], occupancies["40chol"]], OCCUPANCY_FIELDS)

    make_line_plot(rows_by_system, paths["figures"] / "temperature_20_vs_40.png", "Temperature vs time", "temperature (K)", "temperature_K")
    make_line_plot(
        rows_by_system,
        paths["figures"] / "potential_energy_20_vs_40.png",
        "Potential energy vs time",
        "potential energy (kJ/mol)",
        "potential_energy_kJ_mol",
    )
    make_line_plot(rows_by_system, paths["figures"] / "ligand_rmsd_20_vs_40.png", "Ligand RMSD vs time", "ligand RMSD (A)", "ligand_rmsd_A")
    make_line_plot(
        rows_by_system,
        paths["figures"] / "key_min_distance_20_vs_40.png",
        "Ligand minimum key-residue distance vs time",
        "key-min distance (A)",
        "key_min_distance_A",
    )
    make_per_residue_distance_plot(rows_by_system, paths["figures"] / "per_residue_distance_20_vs_40.png")
    make_line_plot(rows_by_system, paths["figures"] / "near_chains_20_vs_40.png", "Ligand near-chain count vs time", "near chains", "near_chains")
    make_clash_plot(rows_by_system, paths["figures"] / "clash_counts_20_vs_40.png")
    make_occupancy_barplot(occupancies, paths["figures"] / "contact_occupancy_barplot.png")
    make_per_residue_boxplot(rows_by_system, paths["figures"] / "per_residue_distance_boxplot.png")
    make_rmsd_distribution_plot(rows_by_system, paths["figures"] / "rmsd_distribution_20_vs_40.png")
    make_pymol_script(ligand_id, resname, input_paths, paths["pymol_path"])
    write_wrapper_script(ligand_id, resname, paths["scripts"])

    for system in SYSTEMS:
        occupancies[system]["contact_mode_label"] = classify_contact_mode(occupancies[system])
        occupancies[system]["near_chains_all_2"] = summaries[system]["near_chains_all_2"]

    response_label = classify_cholesterol_response(occupancies)
    interpretation = candidate_interpretation(ligand_id, resname, occupancies, response_label)
    contact_style20 = infer_contact_style(str(occupancies["20chol"]["contact_mode_label"]), bool(summaries["20chol"]["near_chains_all_2"]))
    contact_style40 = infer_contact_style(str(occupancies["40chol"]["contact_mode_label"]), bool(summaries["40chol"]["near_chains_all_2"]))

    write_report(
        ligand_id,
        resname,
        input_paths,
        paths,
        summaries["20chol"],
        summaries["40chol"],
        occupancies["20chol"],
        occupancies["40chol"],
        contact_style20,
        contact_style40,
        str(occupancies["20chol"]["contact_mode_label"]),
        str(occupancies["40chol"]["contact_mode_label"]),
        response_label,
        interpretation,
    )

    batch_rows: list[dict[str, object]] = []
    for system in SYSTEMS:
        summary = summaries[system]
        occ = occupancies[system]
        result = result_json_by_system[system]
        batch_rows.append(
            {
                "ligand_id": ligand_id,
                "resname": resname,
                "system": system,
                "completed_100ps": bool(result.get("completed_100ps")),
                "NaN_detected": summary["NaN_detected"],
                "near_chains_all_2": summary["near_chains_all_2"],
                "temperature_mean": summary["temperature_mean"],
                "ligand_rmsd_mean": occ["ligand_rmsd_mean"],
                "ligand_rmsd_final": summary["ligand_rmsd_final"],
                "ligand_rmsd_max": summary["ligand_rmsd_max"],
                "TYR13_occ_lt5A": occ["TYR13_occ_lt5A"],
                "VAL17_occ_lt5A": occ["VAL17_occ_lt5A"],
                "SER20_occ_lt5A": occ["SER20_occ_lt5A"],
                "TYR13_mean_distance": summary["TYR13_distance_mean"],
                "VAL17_mean_distance": summary["VAL17_distance_mean"],
                "SER20_mean_distance": summary["SER20_distance_mean"],
                "key_min_distance_mean": summary["key_min_distance_mean"],
                "key_min_occ_lt5A": occ["key_min_occ_lt5A"],
                "severe_clash_ever": summary["severe_clash_ever"],
                "close_clash_ever": summary["close_clash_ever"],
                "contact_mode_label": occ["contact_mode_label"],
                "cholesterol_response_label": response_label,
                "candidate_interpretation": interpretation,
            }
        )

    return {
        "ligand_id": ligand_id,
        "resname": resname,
        "paths": paths,
        "input_paths": input_paths,
        "summaries": summaries,
        "occupancies": occupancies,
        "response_label": response_label,
        "interpretation": interpretation,
        "batch_rows": batch_rows,
    }


def run_single_ligand(ligand_id: str, resname: str) -> dict[str, object]:
    return analyze_single_ligand(ligand_id, resname)


def main() -> None:
    batch_rows: list[dict[str, object]] = []
    results: list[dict[str, object]] = []
    for ligand_id, resname in CANDIDATES:
        result = analyze_single_ligand(ligand_id, resname)
        batch_rows.extend(result["batch_rows"])
        results.append(result)
    write_csv(BATCH_SUMMARY_CSV, batch_rows, BATCH_FIELDS)

    for result in results:
        ligand_id = result["ligand_id"]
        occ20 = result["occupancies"]["20chol"]
        occ40 = result["occupancies"]["40chol"]
        summary20 = result["summaries"]["20chol"]
        summary40 = result["summaries"]["40chol"]
        print(f"{ligand_id} / {result['resname']}")
        print(f"  20chol completed 100 ps: {bool_text(to_float(summary20['completed_ps']) >= 100)}")
        print(f"  40chol completed 100 ps: {bool_text(to_float(summary40['completed_ps']) >= 100)}")
        print(f"  20chol NaN: {bool_text(summary20['NaN_detected'])}")
        print(f"  40chol NaN: {bool_text(summary40['NaN_detected'])}")
        print(f"  20chol near chains all 2: {bool_text(summary20['near_chains_all_2'])}")
        print(f"  40chol near chains all 2: {bool_text(summary40['near_chains_all_2'])}")
        print(
            "  20chol TYR13/VAL17/SER20 <5 A occupancy: "
            f"{pct_text(occ20['TYR13_occ_lt5A'])} / {pct_text(occ20['VAL17_occ_lt5A'])} / {pct_text(occ20['SER20_occ_lt5A'])}"
        )
        print(
            "  40chol TYR13/VAL17/SER20 <5 A occupancy: "
            f"{pct_text(occ40['TYR13_occ_lt5A'])} / {pct_text(occ40['VAL17_occ_lt5A'])} / {pct_text(occ40['SER20_occ_lt5A'])}"
        )
        print(f"  contact mode in 20chol: {occ20['contact_mode_label']}")
        print(f"  contact mode in 40chol: {occ40['contact_mode_label']}")
        print(f"  cholesterol response label: {result['response_label']}")
        print(f"  recommended interpretation: {result['interpretation']}")

    def score_balanced(result: dict[str, object]) -> float:
        occ20 = result["occupancies"]["20chol"]
        occ40 = result["occupancies"]["40chol"]
        return min(
            to_float(occ20["TYR13_occ_lt5A"]),
            to_float(occ20["VAL17_occ_lt5A"]),
            to_float(occ20["SER20_occ_lt5A"]),
            to_float(occ40["TYR13_occ_lt5A"]),
            to_float(occ40["VAL17_occ_lt5A"]),
            to_float(occ40["SER20_occ_lt5A"]),
        )

    most_balanced = max(results, key=score_balanced)
    shifted = [r for r in results if r["response_label"] == "cholesterol_sensitive_contact_shift"]
    weaker = [r for r in results if r["response_label"] in {"interface_retained_but_weaker_contact", "weak_or_inconclusive"}]
    print(f"most balanced candidate: {most_balanced['ligand_id']} / {most_balanced['resname']}")
    print(
        "cholesterol-sensitive shift candidate: "
        + (", ".join(f"{r['ligand_id']} / {r['resname']}" for r in shifted) if shifted else "none")
    )
    print(
        "weaker-contact candidate: "
        + (", ".join(f"{r['ligand_id']} / {r['resname']}" for r in weaker) if weaker else "none")
    )
    print("report directories:")
    for result in results:
        print(f"  {result['paths']['reports']}")


if __name__ == "__main__":
    main()
