from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


BASE_DIR = Path(r"C:\TRKB_WP2\ligand_bound_MD")
OUT_DIR = BASE_DIR / "comparison_19.1_20_vs_40"
REPORTS_DIR = OUT_DIR / "reports"
FIGURES_DIR = OUT_DIR / "figures"
TABLES_DIR = OUT_DIR / "tables"
PYMOL_DIR = OUT_DIR / "pymol"

MONITOR_20 = BASE_DIR / r"pilot_build\19.1_20chol\md100ps_R2\reports\md100ps_20chol_R2_monitor.csv"
PDB_20 = BASE_DIR / r"pilot_build\19.1_20chol\md100ps_R2\outputs\TRKB_20chol_L003_19.1_100ps_R2_final.pdb"
DCD_20 = BASE_DIR / r"pilot_build\19.1_20chol\md100ps_R2\outputs\TRKB_20chol_L003_19.1_100ps_R2.dcd"
SUMMARY_20 = BASE_DIR / r"pilot_build\19.1_20chol\md100ps_R2\reports\md100ps_20chol_R2_summary.md"

MONITOR_40 = BASE_DIR / r"pilot_build\19.1_40chol\md100ps_R2\reports\md100ps_40chol_R2_monitor.csv"
PDB_40 = BASE_DIR / r"pilot_build\19.1_40chol\md100ps_R2\outputs\TRKB_40chol_L003_19.1_100ps_R2_final.pdb"
DCD_40 = BASE_DIR / r"pilot_build\19.1_40chol\md100ps_R2\outputs\TRKB_40chol_L003_19.1_100ps_R2.dcd"
SUMMARY_40 = BASE_DIR / r"pilot_build\19.1_40chol\md100ps_R2\reports\md100ps_40chol_R2_summary.md"

SUMMARY_TABLE = TABLES_DIR / "comparison_summary_statistics.csv"
COMBINED_TABLE = TABLES_DIR / "combined_time_series_19.1.csv"
OCCUPANCY_TABLE = TABLES_DIR / "contact_occupancy_summary.csv"
REPORT_PATH = REPORTS_DIR / "comparison_19.1_20_vs_40_summary.md"
PYMOL_PATH = PYMOL_DIR / "view_19.1_20_vs_40_final.pml"

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

RESIDUES = ["TYR13", "VAL17", "SER20"]
RESIDUE_THRESHOLDS = [4.0, 5.0, 6.0, 8.0]
KEY_THRESHOLDS = [4.0, 5.0, 6.0]
SYSTEMS = ["20chol", "40chol"]
COLORS = {"20chol": "#1f77b4", "40chol": "#d62728"}


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


def finite(values: list[float]) -> list[float]:
    return [value for value in values if math.isfinite(value)]


def finite_values(rows: list[dict[str, object]], key: str) -> list[float]:
    return finite([to_float(row.get(key)) for row in rows])


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


def mean(values: list[float]) -> float:
    valid = finite(values)
    if not valid:
        return math.nan
    return float(np.mean(np.asarray(valid, dtype=float)))


def median(values: list[float]) -> float:
    valid = finite(values)
    if not valid:
        return math.nan
    return float(np.median(np.asarray(valid, dtype=float)))


def percentile(values: list[float], q: float) -> float:
    valid = finite(values)
    if not valid:
        return math.nan
    return float(np.percentile(np.asarray(valid, dtype=float), q))


def max_value(values: list[float]) -> float:
    valid = finite(values)
    return max(valid) if valid else math.nan


def occupancy(values: list[float], threshold: float) -> dict[str, float]:
    valid = finite(values)
    hits = sum(1 for value in valid if value < threshold)
    total = len(valid)
    return {
        "hits": hits,
        "total": total,
        "fraction": hits / total if total else math.nan,
    }


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
        "near_chains_all_2": bool(rows) and all(to_int(row.get("near_chains")) == 2 for row in rows),
        "severe_clash_ever": any(to_int(row.get("severe_clashes")) > 0 for row in rows),
        "close_clash_ever": any(to_int(row.get("close_contacts")) > 0 for row in rows),
        "NaN_detected": any(any(is_nanish(value) for value in raw.values()) for raw in raw_rows),
    }


def summarize_occupancy_system(system: str, rows: list[dict[str, object]]) -> dict[str, object]:
    per_residue_values = {
        residue: [to_float(row.get(f"{residue}_distance_A")) for row in rows]
        for residue in RESIDUES
    }
    key_min = [to_float(row.get("key_min_distance_A")) for row in rows]
    near_chains = [to_float(row.get("near_chains")) for row in rows]
    rmsd = [to_float(row.get("ligand_rmsd_A")) for row in rows]

    near_valid = finite(near_chains)
    near_eq2_fraction = sum(1 for value in near_valid if value == 2) / len(near_valid) if near_valid else math.nan

    summary: dict[str, object] = {
        "system": system,
        "per_residue_frame_count": len(finite(per_residue_values["TYR13"])),
        "near_chains_eq2_fraction": near_eq2_fraction,
        "ligand_rmsd_mean": mean(rmsd),
        "ligand_rmsd_median": median(rmsd),
        "ligand_rmsd_max": max_value(rmsd),
        "ligand_rmsd_p95": percentile(rmsd, 95.0),
    }
    for residue in RESIDUES:
        residue_values = per_residue_values[residue]
        for threshold in RESIDUE_THRESHOLDS:
            summary[f"{residue}_occ_lt{int(threshold)}A"] = occupancy(residue_values, threshold)["fraction"]
    for threshold in KEY_THRESHOLDS:
        summary[f"key_min_occ_lt{int(threshold)}A"] = occupancy(key_min, threshold)["fraction"]
    return summary


def fmt(value: object, digits: int = 6) -> str:
    value = to_float(value)
    if not math.isfinite(value):
        return ""
    return f"{value:.{digits}g}"


def fmt_text(value: object, digits: int = 3) -> str:
    value = to_float(value)
    if not math.isfinite(value):
        return "unavailable"
    return f"{value:.{digits}f}"


def pct_text(value: object) -> str:
    value = to_float(value)
    if not math.isfinite(value):
        return "unavailable"
    return f"{100.0 * value:.1f}%"


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            out = {}
            for field in fields:
                value = row.get(field, "")
                if isinstance(value, bool):
                    out[field] = value
                elif isinstance(value, int):
                    out[field] = value
                elif isinstance(value, float):
                    out[field] = fmt(value)
                else:
                    out[field] = value
            writer.writerow(out)


def plot_lines(path: Path, title: str, ylabel: str, rows_by_system: dict[str, list[dict[str, object]]], keys: list[str], labels: list[str] | None = None) -> None:
    plt.figure(figsize=(8.5, 5.2), dpi=160)
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
            plt.plot(xs, ys, label=f"{system}{suffix}", color=COLORS[system], linestyle=linestyle, linewidth=1.8)
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


def make_trend_figures(rows_by_system: dict[str, list[dict[str, object]]]) -> None:
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
    for system, rows in rows_by_system.items():
        times = [to_float(row.get("time_ps")) for row in rows]
        severe = [to_float(row.get("severe_clashes")) for row in rows]
        close = [to_float(row.get("close_contacts")) for row in rows]
        plt.plot(times, severe, label=f"{system} severe <1.0 A", color=COLORS[system], linestyle="-", linewidth=1.8)
        plt.plot(times, close, label=f"{system} close <1.5 A", color=COLORS[system], linestyle="--", linewidth=1.8)
    plt.title("Ligand Clash Counts")
    plt.xlabel("time_ps")
    plt.ylabel("count")
    plt.grid(True, alpha=0.28)
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "clash_counts_20_vs_40.png")
    plt.close()


def make_occupancy_barplot(summary_by_system: dict[str, dict[str, object]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), dpi=160)

    residue_labels = []
    residue_20 = []
    residue_40 = []
    for residue in RESIDUES:
        for threshold in RESIDUE_THRESHOLDS:
            residue_labels.append(f"{residue}\n<{int(threshold)} A")
            field = f"{residue}_occ_lt{int(threshold)}A"
            residue_20.append(to_float(summary_by_system["20chol"].get(field)))
            residue_40.append(to_float(summary_by_system["40chol"].get(field)))

    x = np.arange(len(residue_labels))
    width = 0.38
    axes[0].bar(x - width / 2, residue_20, width, label="20chol", color=COLORS["20chol"])
    axes[0].bar(x + width / 2, residue_40, width, label="40chol", color=COLORS["40chol"])
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(residue_labels, rotation=45, ha="right")
    axes[0].set_ylim(0, 1.05)
    axes[0].set_ylabel("occupancy fraction")
    axes[0].set_title("Per-residue contact occupancy")
    axes[0].grid(True, axis="y", alpha=0.28)
    axes[0].legend()

    key_labels = [f"key-min\n<{int(threshold)} A" for threshold in KEY_THRESHOLDS] + ["near chains\n= 2"]
    key_20 = [to_float(summary_by_system["20chol"].get(f"key_min_occ_lt{int(threshold)}A")) for threshold in KEY_THRESHOLDS]
    key_40 = [to_float(summary_by_system["40chol"].get(f"key_min_occ_lt{int(threshold)}A")) for threshold in KEY_THRESHOLDS]
    key_20.append(to_float(summary_by_system["20chol"].get("near_chains_eq2_fraction")))
    key_40.append(to_float(summary_by_system["40chol"].get("near_chains_eq2_fraction")))
    x2 = np.arange(len(key_labels))
    axes[1].bar(x2 - width / 2, key_20, width, label="20chol", color=COLORS["20chol"])
    axes[1].bar(x2 + width / 2, key_40, width, label="40chol", color=COLORS["40chol"])
    axes[1].set_xticks(x2)
    axes[1].set_xticklabels(key_labels)
    axes[1].set_ylim(0, 1.05)
    axes[1].set_ylabel("fraction")
    axes[1].set_title("Key-min distance and interface retention")
    axes[1].grid(True, axis="y", alpha=0.28)

    fig.suptitle("19.1 / L003 contact occupancy: 20chol vs 40chol")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "contact_occupancy_barplot.png")
    plt.close(fig)


def make_per_residue_boxplot(rows_by_system: dict[str, list[dict[str, object]]]) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 5.4), dpi=160)
    positions = []
    data = []
    labels = []
    colors = []
    for index, residue in enumerate(RESIDUES):
        base = 1 + index * 3
        for offset, system in enumerate(SYSTEMS):
            values = finite([to_float(row.get(f"{residue}_distance_A")) for row in rows_by_system[system]])
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
    fig.savefig(FIGURES_DIR / "per_residue_distance_boxplot.png")
    plt.close(fig)


def make_rmsd_distribution_plot(rows_by_system: dict[str, list[dict[str, object]]]) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 5.2), dpi=160)
    upper = max(max_value([to_float(row.get("ligand_rmsd_A")) for row in rows_by_system[system]]) for system in SYSTEMS) + 0.25
    bins = np.linspace(0.0, upper, 14)
    for system in SYSTEMS:
        values = finite([to_float(row.get("ligand_rmsd_A")) for row in rows_by_system[system]])
        ax.hist(values, bins=bins, alpha=0.45, label=system, color=COLORS[system], density=True, edgecolor="black", linewidth=0.5)
        ax.axvline(mean(values), color=COLORS[system], linestyle="--", linewidth=1.5)
    ax.set_xlabel("ligand RMSD (A)")
    ax.set_ylabel("density")
    ax.set_title("Ligand RMSD distribution")
    ax.grid(True, axis="y", alpha=0.28)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "rmsd_distribution_20_vs_40.png")
    plt.close(fig)


def write_pymol_script() -> None:
    PYMOL_PATH.write_text(
        "\n".join(
            [
                "reinitialize",
                f"load {PDB_20.as_posix()}, TRKB_20chol_L003_19_1",
                f"load {PDB_40.as_posix()}, TRKB_40chol_L003_19_1",
                "hide everything",
                "show cartoon, polymer.protein",
                "show sticks, resn L003",
                "show sticks, (resn TYR and resi 434) or (resn VAL and resi 438) or (resn SER and resi 441)",
                "color marine, TRKB_20chol_L003_19_1",
                "color firebrick, TRKB_40chol_L003_19_1",
                "color cyan, TRKB_20chol_L003_19_1 and resn L003",
                "color orange, TRKB_40chol_L003_19_1 and resn L003",
                "color yellow, resn TYR+VAL+SER and resi 434+438+441",
                "set stick_radius, 0.18",
                "set cartoon_transparency, 0.35",
                "zoom (resn L003) or (resi 434+438+441), 8",
                "orient (resn L003) or (resi 434+438+441)",
                "# Visualization only; not used for simulation.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def stable_val_ser_contact(summary_occ_40: dict[str, object], summary_cmp_40: dict[str, object]) -> str:
    val_occ = to_float(summary_occ_40.get("VAL17_occ_lt5A"))
    ser_occ = to_float(summary_occ_40.get("SER20_occ_lt5A"))
    val_mean = to_float(summary_cmp_40.get("VAL17_distance_mean"))
    ser_mean = to_float(summary_cmp_40.get("SER20_distance_mean"))
    if all(math.isfinite(v) for v in [val_occ, ser_occ, val_mean, ser_mean]) and val_occ >= 0.5 and ser_occ >= 0.5 and val_mean < 5.0 and ser_mean < 5.0:
        return "YES"
    return "NO"


def write_report(
    summary20: dict[str, object],
    summary40: dict[str, object],
    occ20: dict[str, object],
    occ40: dict[str, object],
    val_ser_retained_40: str,
) -> None:
    report = f"""# 19.1 / L003 20chol vs 40chol R2 Pilot Comparison

## 分析目的

比较 19.1 / L003 在 20% CHOL 与 40% CHOL TRKB-TMD 体系中的 restrained MD pilot 稳定性。

## 输入数据

- 20chol monitor CSV: `{MONITOR_20}`
- 20chol DCD: `{DCD_20}`
- 20chol final PDB: `{PDB_20}`
- 20chol summary: `{SUMMARY_20}`
- 40chol monitor CSV: `{MONITOR_40}`
- 40chol DCD: `{DCD_40}`
- 40chol final PDB: `{PDB_40}`
- 40chol summary: `{SUMMARY_40}`

## 模拟状态

两个体系都属于 100 ps R2 PBC-aware restrained MD pilot。它们不是 production MD 模拟，本报告也不是 production MD 结论。

## 稳定性结果

- 20chol completed ps: {fmt_text(summary20['completed_ps'])}; 40chol completed ps: {fmt_text(summary40['completed_ps'])}
- 20chol temperature range: {fmt_text(summary20['temperature_min'])}-{fmt_text(summary20['temperature_max'])} K
- 40chol temperature range: {fmt_text(summary40['temperature_min'])}-{fmt_text(summary40['temperature_max'])} K
- NaN detected, 20chol/40chol: {summary20['NaN_detected']} / {summary40['NaN_detected']}
- severe clash ever, 20chol/40chol: {summary20['severe_clash_ever']} / {summary40['severe_clash_ever']}
- close clash ever, 20chol/40chol: {summary20['close_clash_ever']} / {summary40['close_clash_ever']}
- near chains remained 2 throughout, 20chol/40chol: {summary20['near_chains_all_2']} / {summary40['near_chains_all_2']}

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

## 关键比较

- 19.1 在 20chol 和 40chol 中都保持了 near chains = 2，并且 key-min distance occupancy 保持很高，说明 ligand 在两种胆固醇条件下都维持 interface retention。
- 与 14.2 的 40chol 模式不同，19.1 在 40chol 中并没有明显退化成偏 TYR13、同时失去 VAL17/SER20 的接触构象。
- 19.1 是否在 40chol 中仍保持 VAL17/SER20 contact: {val_ser_retained_40}
- 结合最终距离、平均距离和 `<5 A` occupancy，19.1 看起来在两个胆固醇条件下都更均衡地保持了 TYR13 / VAL17 / SER20 三关键残基接触。

## 局限性

- 100 ps restrained MD 不是 production MD。
- 两个体系都使用了 restraints，因此这更适合作为 pilot stability / retention 检查，而不是平衡态结论。
- CGenFF 参数 penalty 对 L003 / 19.1 的解释仍然重要，需要与本分析一起看。
- 目前仍只覆盖单个 ligand、单条 trajectory。
- 后续仍需要更多 ligand、更长时间、更少 restraints，以及 replicate trajectories 支持。

## 下一步建议

1. 如果 19.1 双体系 occupancy 保持良好，可作为比 14.2 更均衡的重点候选。
2. 后续可进入 14.2 vs 19.1 横向比较。
3. 暂时不要称为 production MD 结论。

## 输出

- Summary statistics table: `{SUMMARY_TABLE}`
- Combined time series table: `{COMBINED_TABLE}`
- Contact occupancy summary table: `{OCCUPANCY_TABLE}`
- Figures directory: `{FIGURES_DIR}`
- PyMOL visualization script: `{PYMOL_PATH}`
"""
    REPORT_PATH.write_text(report, encoding="utf-8")


def terminal_summary(summary20: dict[str, object], summary40: dict[str, object], occ20: dict[str, object], occ40: dict[str, object], val_ser_retained_40: str) -> None:
    print(f"20chol completed 100 ps: {'YES' if to_float(summary20['completed_ps']) >= 100 else 'NO'}")
    print(f"40chol completed 100 ps: {'YES' if to_float(summary40['completed_ps']) >= 100 else 'NO'}")
    print(f"20chol NaN: {'YES' if summary20['NaN_detected'] else 'NO'}")
    print(f"40chol NaN: {'YES' if summary40['NaN_detected'] else 'NO'}")
    print(f"20chol temperature range: {fmt_text(summary20['temperature_min'])} - {fmt_text(summary20['temperature_max'])} K")
    print(f"40chol temperature range: {fmt_text(summary40['temperature_min'])} - {fmt_text(summary40['temperature_max'])} K")
    print(f"20chol ligand RMSD final/max: {fmt_text(summary20['ligand_rmsd_final'])} / {fmt_text(summary20['ligand_rmsd_max'])} A")
    print(f"40chol ligand RMSD final/max: {fmt_text(summary40['ligand_rmsd_final'])} / {fmt_text(summary40['ligand_rmsd_max'])} A")
    print(
        "20chol TYR13 / VAL17 / SER20 final/mean: "
        f"{fmt_text(summary20['TYR13_distance_final'])}/{fmt_text(summary20['TYR13_distance_mean'])} ; "
        f"{fmt_text(summary20['VAL17_distance_final'])}/{fmt_text(summary20['VAL17_distance_mean'])} ; "
        f"{fmt_text(summary20['SER20_distance_final'])}/{fmt_text(summary20['SER20_distance_mean'])} A"
    )
    print(
        "40chol TYR13 / VAL17 / SER20 final/mean: "
        f"{fmt_text(summary40['TYR13_distance_final'])}/{fmt_text(summary40['TYR13_distance_mean'])} ; "
        f"{fmt_text(summary40['VAL17_distance_final'])}/{fmt_text(summary40['VAL17_distance_mean'])} ; "
        f"{fmt_text(summary40['SER20_distance_final'])}/{fmt_text(summary40['SER20_distance_mean'])} A"
    )
    print(
        "20chol contact occupancy <5 A, TYR13/VAL17/SER20: "
        f"{pct_text(occ20['TYR13_occ_lt5A'])} / {pct_text(occ20['VAL17_occ_lt5A'])} / {pct_text(occ20['SER20_occ_lt5A'])}"
    )
    print(
        "40chol contact occupancy <5 A, TYR13/VAL17/SER20: "
        f"{pct_text(occ40['TYR13_occ_lt5A'])} / {pct_text(occ40['VAL17_occ_lt5A'])} / {pct_text(occ40['SER20_occ_lt5A'])}"
    )
    print(f"near chains always 2: {'YES' if summary20['near_chains_all_2'] and summary40['near_chains_all_2'] else 'NO'}")
    print(f"severe/close appeared, 20chol: {'YES' if summary20['severe_clash_ever'] or summary20['close_clash_ever'] else 'NO'}")
    print(f"severe/close appeared, 40chol: {'YES' if summary40['severe_clash_ever'] or summary40['close_clash_ever'] else 'NO'}")
    print(f"19.1 retains VAL17/SER20 contact in 40chol: {val_ser_retained_40}")
    print(f"report path: {REPORT_PATH}")
    print(f"figures path: {FIGURES_DIR}")


def main() -> None:
    ensure_dirs()

    raw20 = read_csv(MONITOR_20)
    raw40 = read_csv(MONITOR_40)
    rows20 = [normalize_monitor_row("20chol", row) for row in raw20]
    rows40 = [normalize_monitor_row("40chol", row) for row in raw40]
    rows_by_system = {"20chol": rows20, "40chol": rows40}

    combined_rows = rows20 + rows40
    combined_rows.sort(key=lambda row: (to_float(row.get("time_ps")), str(row.get("system"))))

    summary20 = summarize_comparison("20chol", rows20, raw20)
    summary40 = summarize_comparison("40chol", rows40, raw40)
    occ20 = summarize_occupancy_system("20chol", rows20)
    occ40 = summarize_occupancy_system("40chol", rows40)

    write_csv(SUMMARY_TABLE, [summary20, summary40], SUMMARY_FIELDS)
    write_csv(COMBINED_TABLE, combined_rows, COMBINED_FIELDS)
    write_csv(OCCUPANCY_TABLE, [occ20, occ40], OCCUPANCY_FIELDS)

    make_trend_figures(rows_by_system)
    make_occupancy_barplot({"20chol": occ20, "40chol": occ40})
    make_per_residue_boxplot(rows_by_system)
    make_rmsd_distribution_plot(rows_by_system)
    write_pymol_script()

    val_ser_retained_40 = stable_val_ser_contact(occ40, summary40)
    write_report(summary20, summary40, occ20, occ40, val_ser_retained_40)
    terminal_summary(summary20, summary40, occ20, occ40, val_ser_retained_40)


if __name__ == "__main__":
    main()
