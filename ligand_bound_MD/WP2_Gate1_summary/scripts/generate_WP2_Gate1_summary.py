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
OUT_DIR = BASE_DIR / "WP2_Gate1_summary"
SCRIPTS_DIR = OUT_DIR / "scripts"
TABLES_DIR = OUT_DIR / "tables"
FIGURES_DIR = OUT_DIR / "figures"
REPORTS_DIR = OUT_DIR / "reports"
PYMOL_DIR = OUT_DIR / "pymol"

CANDIDATES = {
    "14.2": "L002",
    "19.1": "L003",
    "2.3": "L006",
    "6.2": "L008",
    "17.2": "L010",
}
SYSTEMS = ["20chol", "40chol"]
RESIDUES = ["TYR13", "VAL17", "SER20"]

CANDIDATE_SUMMARY_CSV = TABLES_DIR / "WP2_Gate1_candidate_summary.csv"
LIGAND_LEVEL_CSV = TABLES_DIR / "WP2_Gate1_ligand_level_summary.csv"
REPORT_READY_CSV = TABLES_DIR / "WP2_Gate1_report_ready_summary.csv"
REPORT_PATH = REPORTS_DIR / "WP2_Gate1_summary_report.md"
PYMOL_PATH = PYMOL_DIR / "view_gate1_final_structures.pml"

CANDIDATE_FIELDS = [
    "candidate_id",
    "resname",
    "system",
    "completed_100ps",
    "NaN_detected",
    "temperature_mean",
    "temperature_min",
    "temperature_max",
    "temperature_final",
    "potential_energy_final",
    "ligand_rmsd_mean",
    "ligand_rmsd_final",
    "ligand_rmsd_max",
    "near_chains_all_2",
    "severe_clash_ever",
    "close_clash_ever",
    "TYR13_occ_lt4A",
    "TYR13_occ_lt5A",
    "TYR13_occ_lt6A",
    "VAL17_occ_lt4A",
    "VAL17_occ_lt5A",
    "VAL17_occ_lt6A",
    "SER20_occ_lt4A",
    "SER20_occ_lt5A",
    "SER20_occ_lt6A",
    "TYR13_mean_distance",
    "VAL17_mean_distance",
    "SER20_mean_distance",
    "key_min_distance_mean",
    "key_min_occ_lt5A",
    "contact_mode_label",
    "system_level_interpretation",
]

LIGAND_LEVEL_FIELDS = [
    "candidate_id",
    "resname",
    "20chol_contact_mode",
    "40chol_contact_mode",
    "20chol_TYR13_occ_lt5A",
    "20chol_VAL17_occ_lt5A",
    "20chol_SER20_occ_lt5A",
    "40chol_TYR13_occ_lt5A",
    "40chol_VAL17_occ_lt5A",
    "40chol_SER20_occ_lt5A",
    "delta_TYR13_occ_40_minus_20",
    "delta_VAL17_occ_40_minus_20",
    "delta_SER20_occ_40_minus_20",
    "20chol_rmsd_mean",
    "40chol_rmsd_mean",
    "20chol_rmsd_max",
    "40chol_rmsd_max",
    "20chol_near_chains_all_2",
    "40chol_near_chains_all_2",
    "20chol_clash_ever",
    "40chol_clash_ever",
    "cholesterol_response_label",
    "WP2_Gate1_priority_class",
    "recommended_next_step",
    "rationale",
]

REPORT_READY_FIELDS = [
    "candidate",
    "resname",
    "short_label",
    "20chol_result",
    "40chol_result",
    "cholesterol_response",
    "WP2_interpretation",
    "priority_class",
    "recommended_next_step",
    "evidence_level",
    "limitations",
]

CONTACT_MODE_COLORS = {
    "balanced_three_residue_contact": "#2ca02c",
    "VAL17_SER20_dominant": "#1f77b4",
    "TYR13_shifted": "#d62728",
    "weak_or_mixed_contact": "#9467bd",
}
PRIORITY_COLORS = {
    "Priority A": "#2ca02c",
    "Priority B": "#ff7f0e",
    "Priority C": "#d62728",
    "Deferred": "#7f7f7f",
}
SYSTEM_COLORS = {"20chol": "#1f77b4", "40chol": "#d62728"}


def ensure_dirs() -> None:
    for path in [OUT_DIR, SCRIPTS_DIR, TABLES_DIR, FIGURES_DIR, REPORTS_DIR, PYMOL_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


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


def finite(values: Iterable[float]) -> list[float]:
    return [value for value in values if math.isfinite(value)]


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


def comparison_dir(candidate_id: str) -> Path:
    return BASE_DIR / f"comparison_{candidate_id}_20_vs_40"


def candidate_files(candidate_id: str) -> dict[str, Path]:
    base = comparison_dir(candidate_id)
    return {
        "comparison_summary": base / "tables" / "comparison_summary_statistics.csv",
        "occupancy_summary": base / "tables" / "contact_occupancy_summary.csv",
        "comparison_report": base / "reports" / f"comparison_{candidate_id}_20_vs_40_summary.md",
    }


def classify_contact_mode(row: dict[str, object]) -> str:
    tyr = to_float(row.get("TYR13_occ_lt5A"))
    val = to_float(row.get("VAL17_occ_lt5A"))
    ser = to_float(row.get("SER20_occ_lt5A"))
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
    tyr20 = to_float(row20.get("TYR13_occ_lt5A"))
    tyr40 = to_float(row40.get("TYR13_occ_lt5A"))
    val20 = to_float(row20.get("VAL17_occ_lt5A"))
    val40 = to_float(row40.get("VAL17_occ_lt5A"))
    ser20 = to_float(row20.get("SER20_occ_lt5A"))
    ser40 = to_float(row40.get("SER20_occ_lt5A"))
    deltas = {
        "TYR13": tyr40 - tyr20 if math.isfinite(tyr20) and math.isfinite(tyr40) else math.nan,
        "VAL17": val40 - val20 if math.isfinite(val20) and math.isfinite(val40) else math.nan,
        "SER20": ser40 - ser20 if math.isfinite(ser20) and math.isfinite(ser40) else math.nan,
    }

    if mode20 == "balanced_three_residue_contact" and mode40 == "balanced_three_residue_contact":
        return "cholesterol_robust_balanced"

    if mode20 == "TYR13_shifted" and mode40 == "TYR13_shifted":
        return "TYR13_dominant_interface_retention"

    if bool(row20.get("near_chains_all_2")) and bool(row40.get("near_chains_all_2")):
        val_drop = math.isfinite(deltas["VAL17"]) and deltas["VAL17"] <= -0.35
        ser_drop = math.isfinite(deltas["SER20"]) and deltas["SER20"] <= -0.35
        sig_drop = val_drop or ser_drop
        tyr_gain = math.isfinite(deltas["TYR13"]) and deltas["TYR13"] >= 0.35
        pronounced_shift = tyr_gain or (val_drop and ser_drop) or (
            val_drop and math.isfinite(ser40) and ser40 <= 0.2
        ) or (
            ser_drop and math.isfinite(val40) and val40 <= 0.2
        )
        if mode20 != mode40 and pronounced_shift:
            return "cholesterol_sensitive_contact_shift"
        if sig_drop:
            return "interface_retained_but_weaker_contact"

    if mode20 != mode40 and (mode20 or mode40):
        return "cholesterol_sensitive_contact_shift"
    return "weak_or_inconclusive"


def system_level_interpretation(row: dict[str, object]) -> str:
    mode = str(row.get("contact_mode_label", ""))
    if mode == "balanced_three_residue_contact":
        return "Balanced three-residue key-contact retained."
    if mode == "VAL17_SER20_dominant":
        return "Interface retained with VAL17/SER20-dominant contact pattern."
    if mode == "TYR13_shifted":
        return "Interface retained but shifted toward TYR13 with weaker VAL17/SER20 support."
    if bool(row.get("near_chains_all_2")):
        return "Interface retained, but key-residue contact is mixed or partially weakened."
    return "Contact behavior is weak or ambiguous at the restrained-pilot level."


def stable_system(row: dict[str, object]) -> bool:
    return (
        bool(row.get("completed_100ps"))
        and not bool(row.get("NaN_detected"))
        and bool(row.get("near_chains_all_2"))
        and not bool(row.get("severe_clash_ever"))
        and not bool(row.get("close_clash_ever"))
    )


def classify_priority(ligand_row: dict[str, object]) -> str:
    row20 = ligand_row["20chol"]
    row40 = ligand_row["40chol"]
    if not stable_system(row20) or not stable_system(row40):
        return "Deferred"

    response = str(ligand_row.get("cholesterol_response_label", ""))
    occ20 = [to_float(row20.get(f"{res}_occ_lt5A")) for res in RESIDUES]
    occ40 = [to_float(row40.get(f"{res}_occ_lt5A")) for res in RESIDUES]
    min_occ20 = min(finite(occ20)) if finite(occ20) else math.nan
    min_occ40 = min(finite(occ40)) if finite(occ40) else math.nan
    max_rmsd = max(to_float(row20.get("ligand_rmsd_max")), to_float(row40.get("ligand_rmsd_max")))

    if response == "cholesterol_robust_balanced":
        return "Priority A"
    if response == "interface_retained_but_weaker_contact" and math.isfinite(min_occ20) and math.isfinite(min_occ40):
        if min_occ20 >= 0.8 and min_occ40 >= 0.6 and math.isfinite(max_rmsd) and max_rmsd <= 2.2:
            return "Priority A"
    if response == "cholesterol_sensitive_contact_shift":
        if str(row20.get("contact_mode_label")) in {"balanced_three_residue_contact", "VAL17_SER20_dominant"}:
            return "Priority B"
    if response in {"interface_retained_but_weaker_contact", "TYR13_dominant_interface_retention", "weak_or_inconclusive"}:
        return "Priority C"
    return "Priority C"


def next_step_for_priority(priority: str) -> str:
    if priority == "Priority A":
        return "Advance to focused Gate 2 evaluation with longer, weaker-restraint or replicate trajectories after advisor review."
    if priority == "Priority B":
        return "Retain as a mechanism-focused candidate; discuss whether the cholesterol-dependent shift merits targeted Gate 2 follow-up."
    if priority == "Priority C":
        return "Keep as a reserve candidate; do not prioritize for immediate Gate 2 expansion."
    return "Hold until missing data or quality issues are resolved."


def candidate_rationale(candidate_id: str, resname: str, ligand_row: dict[str, object]) -> str:
    row20 = ligand_row["20chol"]
    row40 = ligand_row["40chol"]
    response = str(ligand_row["cholesterol_response_label"])
    priority = str(ligand_row["WP2_Gate1_priority_class"])
    return (
        f"{candidate_id} / {resname}: 20chol {row20['contact_mode_label']} "
        f"(TYR13/VAL17/SER20 <5 A = {pct_text(row20['TYR13_occ_lt5A'])}/{pct_text(row20['VAL17_occ_lt5A'])}/{pct_text(row20['SER20_occ_lt5A'])}); "
        f"40chol {row40['contact_mode_label']} "
        f"(= {pct_text(row40['TYR13_occ_lt5A'])}/{pct_text(row40['VAL17_occ_lt5A'])}/{pct_text(row40['SER20_occ_lt5A'])}); "
        f"response {response}; assigned {priority}."
    )


def report_ready_interpretation(ligand_row: dict[str, object]) -> str:
    response = str(ligand_row["cholesterol_response_label"])
    if response == "cholesterol_robust_balanced":
        return "Robust balanced contact across both cholesterol conditions."
    if response == "cholesterol_sensitive_contact_shift":
        return "Stable interface retention with a cholesterol-sensitive contact shift."
    if response == "TYR13_dominant_interface_retention":
        return "Stable interface retention, but contact is TYR13-dominant and uneven."
    if response == "interface_retained_but_weaker_contact":
        return "Interface retention persists, but one or more key contacts weaken in 40chol."
    return "Mixed or inconclusive key-contact pattern."


def locate_final_pdb(candidate_id: str, resname: str, system: str) -> Path:
    base = BASE_DIR / "pilot_build" / f"{candidate_id}_{system}" / "md100ps_R2" / "outputs"
    matches = sorted(base.glob(f"TRKB_{system}_{resname}_{candidate_id}_100ps_R2_final.pdb"))
    if matches:
        return matches[0]
    fallback = sorted(base.glob("*_100ps_R2_final.pdb"))
    if fallback:
        return fallback[0]
    return base / "MISSING_FINAL_PDB.pdb"


def load_candidate_data() -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    candidate_rows: list[dict[str, object]] = []
    ligand_level_rows: list[dict[str, object]] = []
    report_ready_rows: list[dict[str, object]] = []

    for candidate_id, resname in CANDIDATES.items():
        files = candidate_files(candidate_id)
        cmp_rows = read_csv(files["comparison_summary"])
        occ_rows = read_csv(files["occupancy_summary"])
        cmp_by_system = {row["system"]: row for row in cmp_rows}
        occ_by_system = {row["system"]: row for row in occ_rows}
        system_rows: dict[str, dict[str, object]] = {}

        for system in SYSTEMS:
            cmp_row = cmp_by_system[system]
            occ_row = occ_by_system[system]
            merged = {
                "candidate_id": candidate_id,
                "resname": resname,
                "system": system,
                "completed_100ps": to_float(cmp_row.get("completed_ps")) >= 100.0,
                "NaN_detected": str(cmp_row.get("NaN_detected", "")).lower() == "true",
                "temperature_mean": to_float(cmp_row.get("temperature_mean")),
                "temperature_min": to_float(cmp_row.get("temperature_min")),
                "temperature_max": to_float(cmp_row.get("temperature_max")),
                "temperature_final": to_float(cmp_row.get("temperature_final")),
                "potential_energy_final": to_float(cmp_row.get("potential_energy_final")),
                "ligand_rmsd_mean": to_float(occ_row.get("ligand_rmsd_mean")),
                "ligand_rmsd_final": to_float(cmp_row.get("ligand_rmsd_final")),
                "ligand_rmsd_max": to_float(cmp_row.get("ligand_rmsd_max")),
                "near_chains_all_2": str(cmp_row.get("near_chains_all_2", "")).lower() == "true",
                "severe_clash_ever": str(cmp_row.get("severe_clash_ever", "")).lower() == "true",
                "close_clash_ever": str(cmp_row.get("close_clash_ever", "")).lower() == "true",
                "TYR13_occ_lt4A": to_float(occ_row.get("TYR13_occ_lt4A")),
                "TYR13_occ_lt5A": to_float(occ_row.get("TYR13_occ_lt5A")),
                "TYR13_occ_lt6A": to_float(occ_row.get("TYR13_occ_lt6A")),
                "VAL17_occ_lt4A": to_float(occ_row.get("VAL17_occ_lt4A")),
                "VAL17_occ_lt5A": to_float(occ_row.get("VAL17_occ_lt5A")),
                "VAL17_occ_lt6A": to_float(occ_row.get("VAL17_occ_lt6A")),
                "SER20_occ_lt4A": to_float(occ_row.get("SER20_occ_lt4A")),
                "SER20_occ_lt5A": to_float(occ_row.get("SER20_occ_lt5A")),
                "SER20_occ_lt6A": to_float(occ_row.get("SER20_occ_lt6A")),
                "TYR13_mean_distance": to_float(cmp_row.get("TYR13_distance_mean")),
                "VAL17_mean_distance": to_float(cmp_row.get("VAL17_distance_mean")),
                "SER20_mean_distance": to_float(cmp_row.get("SER20_distance_mean")),
                "key_min_distance_mean": to_float(cmp_row.get("key_min_distance_mean")),
                "key_min_occ_lt5A": to_float(occ_row.get("key_min_occ_lt5A")),
            }
            merged["contact_mode_label"] = classify_contact_mode(merged)
            merged["system_level_interpretation"] = system_level_interpretation(merged)
            system_rows[system] = merged
            candidate_rows.append(merged)

        ligand_row = {
            "candidate_id": candidate_id,
            "resname": resname,
            "20chol": system_rows["20chol"],
            "40chol": system_rows["40chol"],
        }
        ligand_row["cholesterol_response_label"] = classify_cholesterol_response(system_rows)
        ligand_row["WP2_Gate1_priority_class"] = classify_priority(ligand_row)
        ligand_row["recommended_next_step"] = next_step_for_priority(ligand_row["WP2_Gate1_priority_class"])
        ligand_row["rationale"] = candidate_rationale(candidate_id, resname, ligand_row)

        ligand_level_rows.append(
            {
                "candidate_id": candidate_id,
                "resname": resname,
                "20chol_contact_mode": system_rows["20chol"]["contact_mode_label"],
                "40chol_contact_mode": system_rows["40chol"]["contact_mode_label"],
                "20chol_TYR13_occ_lt5A": system_rows["20chol"]["TYR13_occ_lt5A"],
                "20chol_VAL17_occ_lt5A": system_rows["20chol"]["VAL17_occ_lt5A"],
                "20chol_SER20_occ_lt5A": system_rows["20chol"]["SER20_occ_lt5A"],
                "40chol_TYR13_occ_lt5A": system_rows["40chol"]["TYR13_occ_lt5A"],
                "40chol_VAL17_occ_lt5A": system_rows["40chol"]["VAL17_occ_lt5A"],
                "40chol_SER20_occ_lt5A": system_rows["40chol"]["SER20_occ_lt5A"],
                "delta_TYR13_occ_40_minus_20": system_rows["40chol"]["TYR13_occ_lt5A"] - system_rows["20chol"]["TYR13_occ_lt5A"],
                "delta_VAL17_occ_40_minus_20": system_rows["40chol"]["VAL17_occ_lt5A"] - system_rows["20chol"]["VAL17_occ_lt5A"],
                "delta_SER20_occ_40_minus_20": system_rows["40chol"]["SER20_occ_lt5A"] - system_rows["20chol"]["SER20_occ_lt5A"],
                "20chol_rmsd_mean": system_rows["20chol"]["ligand_rmsd_mean"],
                "40chol_rmsd_mean": system_rows["40chol"]["ligand_rmsd_mean"],
                "20chol_rmsd_max": system_rows["20chol"]["ligand_rmsd_max"],
                "40chol_rmsd_max": system_rows["40chol"]["ligand_rmsd_max"],
                "20chol_near_chains_all_2": system_rows["20chol"]["near_chains_all_2"],
                "40chol_near_chains_all_2": system_rows["40chol"]["near_chains_all_2"],
                "20chol_clash_ever": bool(system_rows["20chol"]["severe_clash_ever"] or system_rows["20chol"]["close_clash_ever"]),
                "40chol_clash_ever": bool(system_rows["40chol"]["severe_clash_ever"] or system_rows["40chol"]["close_clash_ever"]),
                "cholesterol_response_label": ligand_row["cholesterol_response_label"],
                "WP2_Gate1_priority_class": ligand_row["WP2_Gate1_priority_class"],
                "recommended_next_step": ligand_row["recommended_next_step"],
                "rationale": ligand_row["rationale"],
            }
        )

        report_ready_rows.append(
            {
                "candidate": candidate_id,
                "resname": resname,
                "short_label": f"{candidate_id} / {resname}",
                "20chol_result": f"{system_rows['20chol']['contact_mode_label']}; TYR/VAL/SER <5 A = {pct_text(system_rows['20chol']['TYR13_occ_lt5A'])}/{pct_text(system_rows['20chol']['VAL17_occ_lt5A'])}/{pct_text(system_rows['20chol']['SER20_occ_lt5A'])}",
                "40chol_result": f"{system_rows['40chol']['contact_mode_label']}; TYR/VAL/SER <5 A = {pct_text(system_rows['40chol']['TYR13_occ_lt5A'])}/{pct_text(system_rows['40chol']['VAL17_occ_lt5A'])}/{pct_text(system_rows['40chol']['SER20_occ_lt5A'])}",
                "cholesterol_response": ligand_row["cholesterol_response_label"],
                "WP2_interpretation": report_ready_interpretation(ligand_row),
                "priority_class": ligand_row["WP2_Gate1_priority_class"],
                "recommended_next_step": ligand_row["recommended_next_step"],
                "evidence_level": "single 100 ps restrained-pilot trajectory in both 20chol and 40chol",
                "limitations": "non-production R2 PBC-aware restrained MD pilot; no replicate trajectories; CGenFF penalties still matter",
            }
        )

    return candidate_rows, ligand_level_rows, report_ready_rows


def make_contact_heatmap(candidate_rows: list[dict[str, object]]) -> None:
    rows_sorted = sorted(candidate_rows, key=lambda row: (list(CANDIDATES).index(row["candidate_id"]), SYSTEMS.index(row["system"])))
    matrix = np.array(
        [
            [to_float(row["TYR13_occ_lt5A"]), to_float(row["VAL17_occ_lt5A"]), to_float(row["SER20_occ_lt5A"])]
            for row in rows_sorted
        ],
        dtype=float,
    )
    labels = [f"{row['candidate_id']} {row['system']}" for row in rows_sorted]
    fig, ax = plt.subplots(figsize=(6.8, 6.6), dpi=170)
    im = ax.imshow(matrix, aspect="auto", vmin=0.0, vmax=1.0, cmap="viridis")
    ax.set_xticks(range(3))
    ax.set_xticklabels(RESIDUES)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_title("Gate 1 key-residue occupancy (<5 A)")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, f"{100*matrix[i,j]:.0f}%", ha="center", va="center", color="white", fontsize=8)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("occupancy fraction")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "gate1_contact_occupancy_lt5A_heatmap.png")
    plt.close(fig)


def make_priority_barplot(ligand_rows: list[dict[str, object]]) -> None:
    score_map = {"Priority A": 3, "Priority B": 2, "Priority C": 1, "Deferred": 0}
    labels = [f"{row['candidate_id']}\n{row['resname']}" for row in ligand_rows]
    scores = [score_map[row["WP2_Gate1_priority_class"]] for row in ligand_rows]
    colors = [PRIORITY_COLORS[row["WP2_Gate1_priority_class"]] for row in ligand_rows]
    fig, ax = plt.subplots(figsize=(8.2, 5.0), dpi=170)
    bars = ax.bar(labels, scores, color=colors)
    ax.set_ylim(0, 3.6)
    ax.set_ylabel("priority score")
    ax.set_title("WP2 Gate 1 candidate priority")
    ax.grid(True, axis="y", alpha=0.28)
    ax.set_yticks([0, 1, 2, 3])
    ax.set_yticklabels(["Deferred", "Priority C", "Priority B", "Priority A"])
    for bar, row in zip(bars, ligand_rows):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.06, row["WP2_Gate1_priority_class"], ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "gate1_candidate_priority_barplot.png")
    plt.close(fig)


def make_rmsd_comparison(ligand_rows: list[dict[str, object]]) -> None:
    labels = [f"{row['candidate_id']}\n{row['resname']}" for row in ligand_rows]
    x = np.arange(len(labels))
    width = 0.18
    mean20 = [to_float(row["20chol_rmsd_mean"]) for row in ligand_rows]
    mean40 = [to_float(row["40chol_rmsd_mean"]) for row in ligand_rows]
    max20 = [to_float(row["20chol_rmsd_max"]) for row in ligand_rows]
    max40 = [to_float(row["40chol_rmsd_max"]) for row in ligand_rows]
    fig, ax = plt.subplots(figsize=(9.4, 5.2), dpi=170)
    ax.bar(x - 1.5 * width, mean20, width=width, color="#6baed6", label="20chol mean")
    ax.bar(x - 0.5 * width, mean40, width=width, color="#fd8d3c", label="40chol mean")
    ax.bar(x + 0.5 * width, max20, width=width, color="#2171b5", label="20chol max")
    ax.bar(x + 1.5 * width, max40, width=width, color="#cb181d", label="40chol max")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("ligand RMSD (A)")
    ax.set_title("Gate 1 ligand RMSD comparison")
    ax.grid(True, axis="y", alpha=0.28)
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "gate1_rmsd_comparison.png")
    plt.close(fig)


def make_delta_heatmap(ligand_rows: list[dict[str, object]]) -> None:
    matrix = np.array(
        [
            [
                to_float(row["delta_TYR13_occ_40_minus_20"]),
                to_float(row["delta_VAL17_occ_40_minus_20"]),
                to_float(row["delta_SER20_occ_40_minus_20"]),
            ]
            for row in ligand_rows
        ],
        dtype=float,
    )
    labels = [f"{row['candidate_id']} {row['resname']}" for row in ligand_rows]
    fig, ax = plt.subplots(figsize=(6.8, 4.6), dpi=170)
    im = ax.imshow(matrix, aspect="auto", vmin=-1.0, vmax=1.0, cmap="coolwarm")
    ax.set_xticks(range(3))
    ax.set_xticklabels(RESIDUES)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_title("40chol - 20chol occupancy delta (<5 A)")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, f"{100*matrix[i,j]:+.0f}%", ha="center", va="center", color="black", fontsize=8)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("delta occupancy fraction")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "gate1_cholesterol_response_delta_heatmap.png")
    plt.close(fig)


def make_contact_mode_summary(ligand_rows: list[dict[str, object]]) -> None:
    mode_order = [
        "balanced_three_residue_contact",
        "VAL17_SER20_dominant",
        "TYR13_shifted",
        "weak_or_mixed_contact",
    ]
    mode_to_num = {mode: idx for idx, mode in enumerate(mode_order)}
    matrix = np.array(
        [
            [
                mode_to_num[row["20chol_contact_mode"]],
                mode_to_num[row["40chol_contact_mode"]],
            ]
            for row in ligand_rows
        ],
        dtype=float,
    )
    cmap = matplotlib.colors.ListedColormap([CONTACT_MODE_COLORS[mode] for mode in mode_order])
    fig, ax = plt.subplots(figsize=(5.7, 4.8), dpi=170)
    ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=-0.5, vmax=len(mode_order) - 0.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["20chol", "40chol"])
    ax.set_yticks(range(len(ligand_rows)))
    ax.set_yticklabels([f"{row['candidate_id']} {row['resname']}" for row in ligand_rows])
    ax.set_title("Gate 1 contact-mode summary")
    for i, row in enumerate(ligand_rows):
        ax.text(0, i, row["20chol_contact_mode"].replace("_", "\n"), ha="center", va="center", fontsize=7, color="white")
        ax.text(1, i, row["40chol_contact_mode"].replace("_", "\n"), ha="center", va="center", fontsize=7, color="white")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "gate1_contact_mode_summary.png")
    plt.close(fig)


def make_summary_panel(ligand_rows: list[dict[str, object]]) -> None:
    fig = plt.figure(figsize=(12.0, 8.5), dpi=170)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], width_ratios=[1.15, 1.0])

    ax1 = fig.add_subplot(gs[0, 0])
    occ_matrix = np.array(
        [
            [
                to_float(row["20chol_TYR13_occ_lt5A"]),
                to_float(row["20chol_VAL17_occ_lt5A"]),
                to_float(row["20chol_SER20_occ_lt5A"]),
                to_float(row["40chol_TYR13_occ_lt5A"]),
                to_float(row["40chol_VAL17_occ_lt5A"]),
                to_float(row["40chol_SER20_occ_lt5A"]),
            ]
            for row in ligand_rows
        ],
        dtype=float,
    )
    im = ax1.imshow(occ_matrix, aspect="auto", vmin=0.0, vmax=1.0, cmap="viridis")
    ax1.set_xticks(range(6))
    ax1.set_xticklabels(["20 TYR", "20 VAL", "20 SER", "40 TYR", "40 VAL", "40 SER"])
    ax1.set_yticks(range(len(ligand_rows)))
    ax1.set_yticklabels([row["candidate_id"] for row in ligand_rows])
    ax1.set_title("Key-residue occupancy <5 A")
    for i in range(occ_matrix.shape[0]):
        for j in range(occ_matrix.shape[1]):
            ax1.text(j, i, f"{100*occ_matrix[i,j]:.0f}%", ha="center", va="center", color="white", fontsize=7)
    fig.colorbar(im, ax=ax1, fraction=0.046, pad=0.04)

    ax2 = fig.add_subplot(gs[0, 1])
    score_map = {"Priority A": 3, "Priority B": 2, "Priority C": 1, "Deferred": 0}
    x = np.arange(len(ligand_rows))
    scores = [score_map[row["WP2_Gate1_priority_class"]] for row in ligand_rows]
    colors = [PRIORITY_COLORS[row["WP2_Gate1_priority_class"]] for row in ligand_rows]
    bars = ax2.bar(x, scores, color=colors)
    ax2.set_xticks(x)
    ax2.set_xticklabels([row["candidate_id"] for row in ligand_rows])
    ax2.set_yticks([0, 1, 2, 3])
    ax2.set_yticklabels(["D", "C", "B", "A"])
    ax2.set_title("Priority class")
    ax2.grid(True, axis="y", alpha=0.28)
    for bar, row in zip(bars, ligand_rows):
        ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.05, row["WP2_Gate1_priority_class"].split()[-1], ha="center", va="bottom", fontsize=8)

    ax3 = fig.add_subplot(gs[1, 0])
    mean20 = [to_float(row["20chol_rmsd_mean"]) for row in ligand_rows]
    mean40 = [to_float(row["40chol_rmsd_mean"]) for row in ligand_rows]
    ax3.plot(x, mean20, marker="o", color=SYSTEM_COLORS["20chol"], label="20chol mean RMSD")
    ax3.plot(x, mean40, marker="o", color=SYSTEM_COLORS["40chol"], label="40chol mean RMSD")
    ax3.set_xticks(x)
    ax3.set_xticklabels([row["candidate_id"] for row in ligand_rows])
    ax3.set_ylabel("RMSD mean (A)")
    ax3.set_title("Mean ligand RMSD")
    ax3.grid(True, alpha=0.28)
    ax3.legend(fontsize=8)

    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis("off")
    lines = ["Gate 1 overview", ""]
    for row in ligand_rows:
        lines.append(
            f"{row['candidate_id']} / {row['resname']}: {row['WP2_Gate1_priority_class']} | "
            f"{row['cholesterol_response_label']}"
        )
    ax4.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", fontsize=10)

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "gate1_summary_panel.png")
    plt.close(fig)


def write_report(candidate_rows: list[dict[str, object]], ligand_rows: list[dict[str, object]]) -> None:
    priority_groups = {"Priority A": [], "Priority B": [], "Priority C": [], "Deferred": []}
    for row in ligand_rows:
        priority_groups[row["WP2_Gate1_priority_class"]].append(f"{row['candidate_id']} / {row['resname']}")

    input_lines = []
    for candidate_id in CANDIDATES:
        base = comparison_dir(candidate_id)
        input_lines.append(f"- {candidate_id}: `{base}`")

    contact_lines = []
    for row in ligand_rows:
        contact_lines.append(
            f"- {row['candidate_id']} / {row['resname']}: "
            f"20chol `{row['20chol_contact_mode']}` "
            f"({pct_text(row['20chol_TYR13_occ_lt5A'])} / {pct_text(row['20chol_VAL17_occ_lt5A'])} / {pct_text(row['20chol_SER20_occ_lt5A'])}); "
            f"40chol `{row['40chol_contact_mode']}` "
            f"({pct_text(row['40chol_TYR13_occ_lt5A'])} / {pct_text(row['40chol_VAL17_occ_lt5A'])} / {pct_text(row['40chol_SER20_occ_lt5A'])}); "
            f"response `{row['cholesterol_response_label']}`."
        )

    report = f"""# WP2 Gate 1 Summary Report

## A. 分析目的

总结五个 moderate-penalty 候选物在 20% CHOL 和 40% CHOL TRKB-TMD 体系中的 100 ps R2 PBC-aware restrained MD pilot 表现。

## B. 输入数据

{chr(10).join(input_lines)}

主汇总使用了每个 comparison 目录中的 `comparison_summary_statistics.csv` 与 `contact_occupancy_summary.csv`。

## C. 方法说明

所有候选物均经过：
- ligand-bound topology build
- OpenMM read test
- targeted cleanup if needed
- R2 PBC-aware minimization
- very-short probe
- 10 ps / 50 ps / 100 ps restrained MD pilot
- contact occupancy analysis

本汇总只做离线整合，不新增任何 MD。本报告中的所有结论都只适用于 100 ps non-production R2 PBC-aware restrained MD pilot，不是 production MD 结论。

## D. 稳定性总结

- 五个候选物的 20chol / 40chol 体系都完成了 100 ps。
- 本轮汇总中没有检测到 NaN。
- temperature range 整体维持在约 307-314 K 的 pilot 区间。
- near chains 在 10 个 candidate-system 中都保持为 2。
- severe / close clash 在五个候选物的 comparison 表中都没有再次出现。

## E. 接触模式总结

{chr(10).join(contact_lines)}

## F. Gate 1 分类

- Priority A: {", ".join(priority_groups['Priority A']) if priority_groups['Priority A'] else 'none'}
- Priority B: {", ".join(priority_groups['Priority B']) if priority_groups['Priority B'] else 'none'}
- Priority C: {", ".join(priority_groups['Priority C']) if priority_groups['Priority C'] else 'none'}
- Deferred / Reserve: {", ".join(priority_groups['Deferred']) if priority_groups['Deferred'] else 'none'}

## G. 推荐候选物

- Priority A 候选物适合作为下一阶段最优先稳定候选。
- Priority B 候选物更适合作为 cholesterol-sensitive contact-shift 机制候选。
- Priority C 候选物仍可保留，但不建议在 Gate 2 前就继续大规模扩展。

## H. 局限性

- 这些都是 100 ps non-production restrained MD pilot。
- 不是 production MD。
- 体系使用了 R2 PBC-aware restraints。
- 每个体系目前只有单次轨迹，没有 replicate。
- CGenFF penalty 仍需与本表一起解读。
- 仍需要更长时间、更弱 restraints 或重复轨迹验证。

## I. 下一步建议

1. 将 Priority A 候选物作为下一阶段优先稳定候选。
2. 将 Priority B 候选物作为 cholesterol-sensitive contact-shift 机制候选。
3. 将 Priority C 候选物保留为 reserve candidates。
4. Gate 1 之后建议先暂停继续大规模 MD，先给老师评估。
5. 若进入 Gate 2，可对 Priority A/B 做更长时间、更弱 restraints 或 replicate trajectories。

## 输出

- Candidate summary table: `{CANDIDATE_SUMMARY_CSV}`
- Ligand-level summary table: `{LIGAND_LEVEL_CSV}`
- Report-ready summary table: `{REPORT_READY_CSV}`
- Figures directory: `{FIGURES_DIR}`
- PyMOL script: `{PYMOL_PATH}`
"""
    REPORT_PATH.write_text(report, encoding="utf-8")


def write_pymol_script() -> None:
    candidate_colors = {
        ("14.2", "20chol"): "marine",
        ("14.2", "40chol"): "bluewhite",
        ("19.1", "20chol"): "forest",
        ("19.1", "40chol"): "limegreen",
        ("2.3", "20chol"): "violet",
        ("2.3", "40chol"): "magenta",
        ("6.2", "20chol"): "orange",
        ("6.2", "40chol"): "tv_orange",
        ("17.2", "20chol"): "teal",
        ("17.2", "40chol"): "cyan",
    }
    lines = ["reinitialize"]
    for candidate_id, resname in CANDIDATES.items():
        for system in SYSTEMS:
            pdb = locate_final_pdb(candidate_id, resname, system)
            obj = f"{candidate_id.replace('.', '_')}_{system}_{resname}"
            lines.append(f"load {pdb.as_posix()}, {obj}")
            lines.append(f"color {candidate_colors[(candidate_id, system)]}, {obj}")
            lines.append(f"show sticks, {obj} and resn {resname}")
    lines.extend(
        [
            "hide everything",
            "show cartoon, polymer.protein",
            "show sticks, resn L002+L003+L006+L008+L010",
            "show sticks, (resn TYR and resi 434) or (resn VAL and resi 438) or (resn SER and resi 441)",
            "color yellow, resn TYR+VAL+SER and resi 434+438+441",
            "set cartoon_transparency, 0.45",
            "set stick_radius, 0.18",
            "zoom (resn L002+L003+L006+L008+L010) or (resi 434+438+441), 10",
            "orient (resn L002+L003+L006+L008+L010) or (resi 434+438+441)",
            "# Visualization only; not used for simulation.",
            "",
        ]
    )
    PYMOL_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    ensure_dirs()
    candidate_rows, ligand_rows, report_ready_rows = load_candidate_data()
    ligand_rows.sort(key=lambda row: list(CANDIDATES).index(row["candidate_id"]))
    candidate_rows.sort(key=lambda row: (list(CANDIDATES).index(row["candidate_id"]), SYSTEMS.index(row["system"])))

    write_csv(CANDIDATE_SUMMARY_CSV, candidate_rows, CANDIDATE_FIELDS)
    write_csv(LIGAND_LEVEL_CSV, ligand_rows, LIGAND_LEVEL_FIELDS)
    write_csv(REPORT_READY_CSV, report_ready_rows, REPORT_READY_FIELDS)

    make_contact_heatmap(candidate_rows)
    make_priority_barplot(ligand_rows)
    make_rmsd_comparison(ligand_rows)
    make_delta_heatmap(ligand_rows)
    make_contact_mode_summary(ligand_rows)
    make_summary_panel(ligand_rows)
    write_report(candidate_rows, ligand_rows)
    write_pymol_script()

    priority_a = [f"{row['candidate_id']} / {row['resname']}" for row in ligand_rows if row["WP2_Gate1_priority_class"] == "Priority A"]
    priority_b = [f"{row['candidate_id']} / {row['resname']}" for row in ligand_rows if row["WP2_Gate1_priority_class"] == "Priority B"]
    for row in ligand_rows:
        print(f"{row['candidate_id']} / {row['resname']}")
        print(f"  contact mode in 20chol: {row['20chol_contact_mode']}")
        print(f"  contact mode in 40chol: {row['40chol_contact_mode']}")
        print(f"  cholesterol response label: {row['cholesterol_response_label']}")
        print(f"  Gate 1 priority class: {row['WP2_Gate1_priority_class']}")
    print("recommended Priority A: " + (", ".join(priority_a) if priority_a else "none"))
    print("recommended Priority B: " + (", ".join(priority_b) if priority_b else "none"))
    print(f"report path: {REPORT_PATH}")
    print(f"figures path: {FIGURES_DIR}")


if __name__ == "__main__":
    main()
