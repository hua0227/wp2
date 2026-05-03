from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any


TRKB_ROOT = Path(r"C:\TRKB_WP2")
MD_ROOT = TRKB_ROOT / "ligand_bound_MD"
CGENFF_ROOT = MD_ROOT / "cgenff_parameterization"
RETURNED_STR_DIR = CGENFF_ROOT / "returned_str"
REPORTS_DIR = CGENFF_ROOT / "reports"
CHECKLIST_CSV = MD_ROOT / "parameterization_needed" / "ligand_parameterization_checklist.csv"

PENALTY_SUMMARY_CSV = REPORTS_DIR / "cgenff_penalty_summary.csv"
PENALTY_REVIEW_MD = REPORTS_DIR / "cgenff_penalty_review.md"

LIGANDS = [
    ("1", "8.1", "L001"),
    ("2", "14.2", "L002"),
    ("3", "19.1", "L003"),
    ("4", "9.2", "L004"),
    ("5", "12.2", "L005"),
    ("6", "2.3", "L006"),
    ("7", "4.3", "L007"),
    ("8", "6.2", "L008"),
    ("9", "8.3", "L009"),
    ("10", "17.2", "L010"),
]

SUMMARY_FIELDS = [
    "ligand_id",
    "resname",
    "str_file",
    "str_found",
    "param_penalty",
    "charge_penalty",
    "bond_penalty",
    "angle_penalty",
    "dihedral_penalty",
    "improper_penalty",
    "max_penalty",
    "n_penalty_ge_10",
    "n_penalty_ge_50",
    "n_penalty_ge_100",
    "penalty_level",
    "needs_manual_review",
    "notes",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.3f}".rstrip("0").rstrip(".")


def max_or_none(values: list[float]) -> float | None:
    return max(values) if values else None


def penalty_level(max_penalty: float | None) -> str:
    if max_penalty is None:
        return "unknown"
    if max_penalty < 10:
        return "low"
    if max_penalty < 50:
        return "moderate"
    if max_penalty < 100:
        return "high"
    return "very_high"


def charge_penalty_from_atom_line(line: str) -> float | None:
    if not line.lstrip().upper().startswith("ATOM ") or "!" not in line:
        return None
    comment = line.split("!", 1)[1]
    match = re.search(r"[-+]?\d+(?:\.\d+)?", comment)
    return float(match.group(0)) if match else None


def parse_str_penalties(text: str) -> dict[str, Any]:
    param_penalty: float | None = None
    charge_penalty: float | None = None
    charge_values: list[float] = []
    section_values = {
        "BONDS": [],
        "ANGLES": [],
        "DIHEDRALS": [],
        "IMPROPERS": [],
    }

    section: str | None = None
    in_param_block = False

    for raw_line in text.splitlines():
        line = raw_line.strip()
        upper = line.upper()

        summary = re.search(
            r"param\s+penalty\s*=\s*([-+]?\d+(?:\.\d+)?)\s*;\s*charge\s+penalty\s*=\s*([-+]?\d+(?:\.\d+)?)",
            line,
            flags=re.IGNORECASE,
        )
        if summary:
            param_penalty = float(summary.group(1))
            charge_penalty = float(summary.group(2))
            continue

        if upper.startswith("READ PARAM"):
            in_param_block = True
            section = None
            continue

        if upper in section_values:
            section = upper if in_param_block else None
            continue

        if upper in {"END", "RETURN"} and in_param_block:
            section = None

        charge_value = charge_penalty_from_atom_line(raw_line)
        if charge_value is not None:
            charge_values.append(charge_value)

        if in_param_block and section in section_values:
            match = re.search(r"\bpenalty\s*=\s*([-+]?\d+(?:\.\d+)?)", line, flags=re.IGNORECASE)
            if match:
                section_values[section].append(float(match.group(1)))

    all_values: list[float] = []
    if param_penalty is not None:
        all_values.append(param_penalty)
    if charge_penalty is not None:
        all_values.append(charge_penalty)
    all_values.extend(charge_values)
    for values in section_values.values():
        all_values.extend(values)

    max_penalty = max_or_none(all_values)
    return {
        "param_penalty": param_penalty,
        "charge_penalty": charge_penalty,
        "bond_penalty": max_or_none(section_values["BONDS"]),
        "angle_penalty": max_or_none(section_values["ANGLES"]),
        "dihedral_penalty": max_or_none(section_values["DIHEDRALS"]),
        "improper_penalty": max_or_none(section_values["IMPROPERS"]),
        "max_penalty": max_penalty,
        "n_penalty_ge_10": sum(1 for value in all_values if value >= 10),
        "n_penalty_ge_50": sum(1 for value in all_values if value >= 50),
        "n_penalty_ge_100": sum(1 for value in all_values if value >= 100),
    }


def summary_row(rank: str, ligand_id: str, resname: str) -> dict[str, Any]:
    str_path = RETURNED_STR_DIR / f"{resname}_{ligand_id}.str"
    if not str_path.exists():
        return {
            "ligand_id": ligand_id,
            "resname": resname,
            "str_file": str(str_path),
            "str_found": "0",
            "param_penalty": "",
            "charge_penalty": "",
            "bond_penalty": "",
            "angle_penalty": "",
            "dihedral_penalty": "",
            "improper_penalty": "",
            "max_penalty": "",
            "n_penalty_ge_10": "",
            "n_penalty_ge_50": "",
            "n_penalty_ge_100": "",
            "penalty_level": "unknown",
            "needs_manual_review": "1",
            "notes": "missing returned str file",
        }

    parsed = parse_str_penalties(str_path.read_text(encoding="utf-8", errors="replace"))
    level = penalty_level(parsed["max_penalty"])
    notes = []
    if parsed["param_penalty"] is None or parsed["charge_penalty"] is None:
        notes.append("summary param/charge penalty not parsed")
    if level in {"high", "very_high"}:
        notes.append("high CGenFF analogy penalty; inspect or optimize parameters before MD")
    elif level == "moderate":
        notes.append("basic validation recommended")
    elif level == "low":
        notes.append("low preliminary penalty")
    else:
        notes.append("penalty could not be parsed")

    return {
        "ligand_id": ligand_id,
        "resname": resname,
        "str_file": str(str_path),
        "str_found": "1",
        "param_penalty": fmt(parsed["param_penalty"]),
        "charge_penalty": fmt(parsed["charge_penalty"]),
        "bond_penalty": fmt(parsed["bond_penalty"]),
        "angle_penalty": fmt(parsed["angle_penalty"]),
        "dihedral_penalty": fmt(parsed["dihedral_penalty"]),
        "improper_penalty": fmt(parsed["improper_penalty"]),
        "max_penalty": fmt(parsed["max_penalty"]),
        "n_penalty_ge_10": parsed["n_penalty_ge_10"],
        "n_penalty_ge_50": parsed["n_penalty_ge_50"],
        "n_penalty_ge_100": parsed["n_penalty_ge_100"],
        "penalty_level": level,
        "needs_manual_review": "1" if level in {"high", "very_high", "unknown"} else "0",
        "notes": "; ".join(notes),
    }


def status_from_row(row: dict[str, Any]) -> str:
    if row["str_found"] != "1":
        return "pending_cgenff"
    if row["penalty_level"] in {"low", "moderate"}:
        return "cgenff_ready_prelim"
    return "cgenff_needs_review"


def update_checklist(summary_rows: list[dict[str, Any]]) -> None:
    checklist = read_csv(CHECKLIST_CSV)
    existing_fields = list(checklist[0].keys())
    wanted_fields = [
        "str_ready",
        "str_path",
        "cgenff_penalty_checked",
        "param_penalty",
        "charge_penalty",
        "max_penalty",
        "penalty_level",
        "needs_manual_review",
        "parameter_status",
    ]
    fieldnames = existing_fields[:]
    for field in wanted_fields:
        if field not in fieldnames:
            fieldnames.append(field)

    by_ligand = {row["ligand_id"]: row for row in summary_rows}
    for row in checklist:
        ligand_id = row.get("ligand_id", "")
        if ligand_id not in by_ligand:
            continue
        summary = by_ligand[ligand_id]
        row["str_ready"] = summary["str_found"]
        row["str_path"] = summary["str_file"]
        row["cgenff_penalty_checked"] = "1" if summary["str_found"] == "1" and summary["penalty_level"] != "unknown" else "0"
        row["param_penalty"] = summary["param_penalty"]
        row["charge_penalty"] = summary["charge_penalty"]
        row["max_penalty"] = summary["max_penalty"]
        row["penalty_level"] = summary["penalty_level"]
        row["needs_manual_review"] = summary["needs_manual_review"]
        row["parameter_status"] = status_from_row(summary)

    write_csv(CHECKLIST_CSV, checklist, fieldnames)


def write_review_markdown(summary_rows: list[dict[str, Any]]) -> None:
    ready = [row for row in summary_rows if row["penalty_level"] in {"low", "moderate"}]
    cautious = [row for row in summary_rows if row["penalty_level"] in {"high", "very_high", "unknown"}]

    lines = [
        "# CGenFF Penalty Review",
        "",
        "This review marks risk based on CGenFF/ParamChem penalty values only. No ligand is automatically discarded.",
        "High penalties do not necessarily mean a ligand cannot be used, but they do require manual inspection and possibly parameter optimization before MD.",
        "",
        "## Penalty Summary",
        "",
        "| Ligand | Resname | Param | Charge | Max | Level | Manual review | Notes |",
        "|---|---|---:|---:|---:|---|---:|---|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {row['ligand_id']} | {row['resname']} | {row['param_penalty']} | {row['charge_penalty']} | "
            f"{row['max_penalty']} | {row['penalty_level']} | {row['needs_manual_review']} | {row['notes']} |"
        )

    lines.extend(["", "## Preliminary OpenMM Read-Test Candidates", ""])
    if ready:
        for row in ready:
            lines.append(f"- {row['ligand_id']} ({row['resname']}): {row['penalty_level']}, max penalty {row['max_penalty']}")
    else:
        lines.append("- None based on current low/moderate cutoff.")

    lines.extend(["", "## Ligands Requiring Caution", ""])
    if cautious:
        for row in cautious:
            lines.append(f"- {row['ligand_id']} ({row['resname']}): {row['penalty_level']}, max penalty {row['max_penalty']}")
    else:
        lines.append("- None.")

    lines.extend(
        [
            "",
            "## Interpretation Notes",
            "",
            "- `low`: max penalty below 10.",
            "- `moderate`: max penalty from 10 to below 50.",
            "- `high`: max penalty from 50 to below 100.",
            "- `very_high`: max penalty 100 or higher.",
            "- `unknown`: penalty could not be parsed or `.str` was missing.",
            "- High penalty is a risk marker, not an automatic exclusion rule.",
        ]
    )
    PENALTY_REVIEW_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    rows = [summary_row(rank, ligand_id, resname) for rank, ligand_id, resname in LIGANDS]
    write_csv(PENALTY_SUMMARY_CSV, rows, SUMMARY_FIELDS)
    update_checklist(rows)
    write_review_markdown(rows)

    found_count = sum(1 for row in rows if row["str_found"] == "1")
    missing = [f"{row['resname']}_{row['ligand_id']}.str" for row in rows if row["str_found"] != "1"]
    high_risk = [row["ligand_id"] for row in rows if row["penalty_level"] in {"high", "very_high", "unknown"}]

    print("Found .str file count:", found_count)
    print("Missing .str files:")
    if missing:
        for item in missing:
            print(item)
    else:
        print("NONE")

    print("\nPenalty summary:")
    for row in rows:
        print(
            row["ligand_id"],
            row["resname"],
            "param=", row["param_penalty"],
            "charge=", row["charge_penalty"],
            "max=", row["max_penalty"],
            "level=", row["penalty_level"],
            "review=", row["needs_manual_review"],
        )

    print("\nhigh / very_high / unknown ligand list:")
    if high_risk:
        for ligand_id in high_risk:
            print(ligand_id)
    else:
        print("NONE")

    print("\nUpdated checklist:")
    print(CHECKLIST_CSV)
    print("\nReview markdown:")
    print(PENALTY_REVIEW_MD)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
