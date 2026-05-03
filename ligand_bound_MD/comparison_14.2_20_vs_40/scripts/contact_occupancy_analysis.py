from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


BASE_DIR = Path(r"C:\TRKB_WP2\ligand_bound_MD\comparison_14.2_20_vs_40")
TABLES_DIR = BASE_DIR / "tables"
FIGURES_DIR = BASE_DIR / "figures"
REPORTS_DIR = BASE_DIR / "reports"

COMBINED_TABLE = TABLES_DIR / "combined_time_series_14.2.csv"
COMPARISON_SUMMARY_TABLE = TABLES_DIR / "comparison_summary_statistics.csv"
OCCUPANCY_TABLE = TABLES_DIR / "contact_occupancy_summary.csv"
REPORT_PATH = REPORTS_DIR / "contact_occupancy_analysis.md"
BARPLOT_PATH = FIGURES_DIR / "contact_occupancy_barplot.png"
BOXPLOT_PATH = FIGURES_DIR / "per_residue_distance_boxplot.png"
RMSD_PLOT_PATH = FIGURES_DIR / "rmsd_distribution_20_vs_40.png"

SUMMARY_FIELDS = [
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


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def finite(values: list[float]) -> list[float]:
    return [value for value in values if math.isfinite(value)]


def occupancy(values: list[float], threshold: float) -> dict[str, float]:
    valid = finite(values)
    hits = sum(1 for value in valid if value < threshold)
    total = len(valid)
    return {
        "hits": hits,
        "total": total,
        "fraction": hits / total if total else math.nan,
    }


def percentile(values: list[float], q: float) -> float:
    valid = finite(values)
    if not valid:
        return math.nan
    return float(np.percentile(np.asarray(valid, dtype=float), q))


def median(values: list[float]) -> float:
    valid = finite(values)
    if not valid:
        return math.nan
    return float(np.median(np.asarray(valid, dtype=float)))


def mean(values: list[float]) -> float:
    valid = finite(values)
    if not valid:
        return math.nan
    return float(np.mean(np.asarray(valid, dtype=float)))


def max_value(values: list[float]) -> float:
    valid = finite(values)
    return max(valid) if valid else math.nan


def load_combined_rows() -> list[dict[str, object]]:
    rows = []
    for row in read_csv(COMBINED_TABLE):
        rows.append(
            {
                "time_ps": to_float(row.get("time_ps")),
                "system": row.get("system", ""),
                "temperature_K": to_float(row.get("temperature_K")),
                "potential_energy_kJ_mol": to_float(row.get("potential_energy_kJ_mol")),
                "ligand_rmsd_A": to_float(row.get("ligand_rmsd_A")),
                "key_min_distance_A": to_float(row.get("key_min_distance_A")),
                "TYR13_distance_A": to_float(row.get("TYR13_distance_A")),
                "VAL17_distance_A": to_float(row.get("VAL17_distance_A")),
                "SER20_distance_A": to_float(row.get("SER20_distance_A")),
                "near_chains": to_float(row.get("near_chains")),
                "severe_clashes": to_float(row.get("severe_clashes")),
                "close_contacts": to_float(row.get("close_contacts")),
                "status": row.get("status", ""),
            }
        )
    return rows


def load_comparison_summary() -> dict[str, dict[str, str]]:
    return {row["system"]: row for row in read_csv(COMPARISON_SUMMARY_TABLE)}


def summarize_system(system: str, rows: list[dict[str, object]]) -> dict[str, object]:
    per_residue_values = {
        residue: [to_float(row.get(f"{residue}_distance_A")) for row in rows]
        for residue in RESIDUES
    }
    key_min = [to_float(row.get("key_min_distance_A")) for row in rows]
    near_chains = [to_float(row.get("near_chains")) for row in rows]
    rmsd = [to_float(row.get("ligand_rmsd_A")) for row in rows]

    near_valid = finite(near_chains)
    near_eq2_fraction = (
        sum(1 for value in near_valid if value == 2) / len(near_valid)
        if near_valid
        else math.nan
    )

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
            key = f"{residue}_occ_lt{int(threshold)}A"
            summary[key] = occupancy(residue_values, threshold)["fraction"]

    for threshold in KEY_THRESHOLDS:
        key = f"key_min_occ_lt{int(threshold)}A"
        summary[key] = occupancy(key_min, threshold)["fraction"]

    return summary


def write_summary_table(rows: list[dict[str, object]]) -> None:
    with OCCUPANCY_TABLE.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        for row in rows:
            out = {}
            for field in SUMMARY_FIELDS:
                value = row.get(field, "")
                if isinstance(value, float):
                    out[field] = "" if not math.isfinite(value) else f"{value:.6g}"
                else:
                    out[field] = value
            writer.writerow(out)


def make_barplot(summary_by_system: dict[str, dict[str, object]]) -> None:
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

    fig.suptitle("14.2 / L002 contact occupancy: 20chol vs 40chol")
    fig.tight_layout()
    fig.savefig(BARPLOT_PATH)
    plt.close(fig)


def make_boxplot(rows_by_system: dict[str, list[dict[str, object]]]) -> None:
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
    fig.savefig(BOXPLOT_PATH)
    plt.close(fig)


def make_rmsd_plot(rows_by_system: dict[str, list[dict[str, object]]]) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 5.2), dpi=160)
    bins = np.linspace(0.0, max(max_value([to_float(row.get("ligand_rmsd_A")) for row in rows_by_system[system]]) for system in SYSTEMS) + 0.25, 14)
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
    fig.savefig(RMSD_PLOT_PATH)
    plt.close(fig)


def fmt(value: object, digits: int = 3) -> str:
    value = to_float(value)
    if not math.isfinite(value):
        return "unavailable"
    return f"{value:.{digits}f}"


def pct(value: object) -> str:
    value = to_float(value)
    if not math.isfinite(value):
        return "unavailable"
    return f"{100.0 * value:.1f}%"


def write_report(summary_by_system: dict[str, dict[str, object]], comparison_summary: dict[str, dict[str, str]]) -> None:
    s20 = summary_by_system["20chol"]
    s40 = summary_by_system["40chol"]
    val_drop_5 = to_float(s20["VAL17_occ_lt5A"]) - to_float(s40["VAL17_occ_lt5A"])
    ser_drop_5 = to_float(s20["SER20_occ_lt5A"]) - to_float(s40["SER20_occ_lt5A"])
    tyr_gain_5 = to_float(s40["TYR13_occ_lt5A"]) - to_float(s20["TYR13_occ_lt5A"])

    report = f"""# Contact Occupancy / Distance Occupancy Analysis

## Scope

This analysis uses only the existing 14.2 / L002 `20chol` and `40chol` 100 ps restrained MD pilot outputs from `{COMBINED_TABLE}` and `{COMPARISON_SUMMARY_TABLE}`.

## Important Interpretation

- This is restrained MD pilot analysis, not production MD.
- 40chol still shows interface retention because key-min distance remains low overall and near chains stay at 2 for the full analyzed trajectory.
- 20chol per-residue occupancies use only frames with real per-residue distances available from the prior DCD-based reconstruction. No missing per-residue values were fabricated.

## Occupancy Summary

- TYR13 occupancy is stronger in 40chol than in 20chol:
  40chol TYR13 `<5 A` occupancy = {pct(s40["TYR13_occ_lt5A"])}, versus 20chol = {pct(s20["TYR13_occ_lt5A"])}.
- VAL17 occupancy drops in 40chol:
  20chol VAL17 `<5 A` occupancy = {pct(s20["VAL17_occ_lt5A"])}, versus 40chol = {pct(s40["VAL17_occ_lt5A"])}.
- SER20 occupancy drops in 40chol:
  20chol SER20 `<5 A` occupancy = {pct(s20["SER20_occ_lt5A"])}, versus 40chol = {pct(s40["SER20_occ_lt5A"])}.
- Key-min-distance occupancy remains high in both systems:
  20chol `<5 A` = {pct(s20["key_min_occ_lt5A"])}, 40chol `<5 A` = {pct(s40["key_min_occ_lt5A"])}.
- Near chains = 2 fraction:
  20chol = {pct(s20["near_chains_eq2_fraction"])}, 40chol = {pct(s40["near_chains_eq2_fraction"])}.

## RMSD Summary

- 20chol RMSD mean / median / max / 95th percentile = {fmt(s20["ligand_rmsd_mean"])} / {fmt(s20["ligand_rmsd_median"])} / {fmt(s20["ligand_rmsd_max"])} / {fmt(s20["ligand_rmsd_p95"])} A.
- 40chol RMSD mean / median / max / 95th percentile = {fmt(s40["ligand_rmsd_mean"])} / {fmt(s40["ligand_rmsd_median"])} / {fmt(s40["ligand_rmsd_max"])} / {fmt(s40["ligand_rmsd_p95"])} A.

## Interpretation

- 40chol retains the ligand at the interface, but the contact pattern shifts.
- VAL17 and SER20 contact occupancies are clearly lower in 40chol than in 20chol, with `<5 A` occupancy drops of {100.0 * val_drop_5:.1f} and {100.0 * ser_drop_5:.1f} percentage points.
- TYR13 becomes more prominent in 40chol, with a `<5 A` occupancy gain of {100.0 * tyr_gain_5:.1f} percentage points.
- Together, these results suggest that 40chol may alter the interfacial contact pose of 14.2 rather than removing interface retention.

## Outputs

- Occupancy summary table: `{OCCUPANCY_TABLE}`
- Contact occupancy bar plot: `{BARPLOT_PATH}`
- Per-residue distance box plot: `{BOXPLOT_PATH}`
- RMSD distribution plot: `{RMSD_PLOT_PATH}`
"""
    REPORT_PATH.write_text(report, encoding="utf-8")


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    combined_rows = load_combined_rows()
    comparison_summary = load_comparison_summary()
    rows_by_system = {
        system: [row for row in combined_rows if row["system"] == system]
        for system in SYSTEMS
    }
    summary_rows = [summarize_system(system, rows_by_system[system]) for system in SYSTEMS]
    summary_by_system = {row["system"]: row for row in summary_rows}

    write_summary_table(summary_rows)
    make_barplot(summary_by_system)
    make_boxplot(rows_by_system)
    make_rmsd_plot(rows_by_system)
    write_report(summary_by_system, comparison_summary)

    print(f"20chol TYR13/VAL17/SER20 <5A occupancy: {pct(summary_by_system['20chol']['TYR13_occ_lt5A'])} / {pct(summary_by_system['20chol']['VAL17_occ_lt5A'])} / {pct(summary_by_system['20chol']['SER20_occ_lt5A'])}")
    print(f"40chol TYR13/VAL17/SER20 <5A occupancy: {pct(summary_by_system['40chol']['TYR13_occ_lt5A'])} / {pct(summary_by_system['40chol']['VAL17_occ_lt5A'])} / {pct(summary_by_system['40chol']['SER20_occ_lt5A'])}")
    print(f"20chol key-min <4/<5/<6A occupancy: {pct(summary_by_system['20chol']['key_min_occ_lt4A'])} / {pct(summary_by_system['20chol']['key_min_occ_lt5A'])} / {pct(summary_by_system['20chol']['key_min_occ_lt6A'])}")
    print(f"40chol key-min <4/<5/<6A occupancy: {pct(summary_by_system['40chol']['key_min_occ_lt4A'])} / {pct(summary_by_system['40chol']['key_min_occ_lt5A'])} / {pct(summary_by_system['40chol']['key_min_occ_lt6A'])}")
    print(f"near chains = 2 fraction, 20chol/40chol: {pct(summary_by_system['20chol']['near_chains_eq2_fraction'])} / {pct(summary_by_system['40chol']['near_chains_eq2_fraction'])}")
    print(f"20chol RMSD mean/median/max/p95: {fmt(summary_by_system['20chol']['ligand_rmsd_mean'])} / {fmt(summary_by_system['20chol']['ligand_rmsd_median'])} / {fmt(summary_by_system['20chol']['ligand_rmsd_max'])} / {fmt(summary_by_system['20chol']['ligand_rmsd_p95'])} A")
    print(f"40chol RMSD mean/median/max/p95: {fmt(summary_by_system['40chol']['ligand_rmsd_mean'])} / {fmt(summary_by_system['40chol']['ligand_rmsd_median'])} / {fmt(summary_by_system['40chol']['ligand_rmsd_max'])} / {fmt(summary_by_system['40chol']['ligand_rmsd_p95'])} A")
    print(f"report path: {REPORT_PATH}")


if __name__ == "__main__":
    main()
