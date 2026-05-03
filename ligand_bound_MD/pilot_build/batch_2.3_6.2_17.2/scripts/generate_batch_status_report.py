from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


BASE_DIR = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build")
BATCH_DIR = BASE_DIR / "batch_2.3_6.2_17.2"
REPORTS_DIR = BATCH_DIR / "reports"
STATUS_CSV = REPORTS_DIR / "batch_gate_status.csv"
STATUS_MD = REPORTS_DIR / "batch_gate_status.md"

CANDIDATES = [
    ("2.3", "L006"),
    ("6.2", "L008"),
    ("17.2", "L010"),
]
SYSTEMS = ["20chol", "40chol"]

STATUS_FIELDS = [
    "ligand_id",
    "resname",
    "system",
    "preflight_ok",
    "psf_pdb_generated",
    "read_test_succeeded",
    "cleanup_executed",
    "minimization_succeeded",
    "very_short_probe_passed",
    "completed_10ps",
    "recommend_50ps",
    "completed_50ps",
    "recommend_100ps",
    "completed_100ps",
    "recommend_comparison_analysis",
    "final_stage",
    "stop_reason",
    "error",
]


def load_json_if_exists(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def result_paths(ligand_id: str, system: str) -> dict[str, Path]:
    root = BASE_DIR / f"{ligand_id}_{system}"
    return {
        "gate": root / "reports" / "pilot_gate_summary.json",
        "md50": root / "md50ps_R2" / "reports" / f"md50ps_{system}_R2_result.json",
        "md100": root / "md100ps_R2" / "reports" / f"md100ps_{system}_R2_result.json",
    }


def summarize_system_status(ligand_id: str, resname: str, system: str, payloads: dict[str, dict[str, Any]]) -> dict[str, Any]:
    gate = payloads.get("gate", {})
    md50 = payloads.get("md50", {})
    md100 = payloads.get("md100", {})

    final_stage = "not_started"
    stop_reason = ""
    error = ""

    if md100:
        if md100.get("completed_100ps"):
            final_stage = "completed_100ps"
        else:
            final_stage = "md50ps_failed_during_md100ps"
        stop_reason = str(md100.get("stop_reason", "") or "")
        error = str(md100.get("error", "") or "")
    elif md50:
        if md50.get("completed_50ps"):
            final_stage = "completed_50ps_only"
        else:
            final_stage = "10ps_failed_during_md50ps"
        stop_reason = str(md50.get("stop_reason", "") or "")
        error = str(md50.get("error", "") or "")
    elif gate:
        if gate.get("completed_10ps"):
            final_stage = "completed_10ps_only"
        else:
            final_stage = "pre_50ps_gate_failed"
        stop_reason = str(gate.get("stop_reason", "") or "")
        error = str(gate.get("error", "") or "")

    return {
        "ligand_id": ligand_id,
        "resname": resname,
        "system": system,
        "preflight_ok": bool(gate.get("preflight_ok", False)),
        "psf_pdb_generated": bool(gate.get("psf_pdb_generated", False)),
        "read_test_succeeded": bool(gate.get("read_test_succeeded", False)),
        "cleanup_executed": bool(gate.get("cleanup_executed", False)),
        "minimization_succeeded": bool(gate.get("minimization_succeeded", False)),
        "very_short_probe_passed": bool(gate.get("very_short_probe_passed", False)),
        "completed_10ps": bool(gate.get("completed_10ps", False)),
        "recommend_50ps": bool(gate.get("recommend_50ps", False)),
        "completed_50ps": bool(md50.get("completed_50ps", False)),
        "recommend_100ps": bool(md50.get("recommend_100ps", False)),
        "completed_100ps": bool(md100.get("completed_100ps", False)),
        "recommend_comparison_analysis": bool(md100.get("recommend_comparison_analysis", False)),
        "final_stage": final_stage,
        "stop_reason": stop_reason,
        "error": error,
    }


def write_csv(rows: list[dict[str, Any]]) -> None:
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    with STATUS_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def write_report(rows: list[dict[str, Any]]) -> None:
    lines = [
        "# Batch Gate Status Report",
        "",
        "This report tracks the gate-based ligand-bound non-production pilot pipeline for 2.3 / L006, 6.2 / L008, and 17.2 / L010 across 20chol and 40chol.",
        "",
        "- Pipeline order: preflight -> psfgen build -> OpenMM read test -> clash diagnosis -> targeted cleanup -> R2 minimization -> very-short probe -> 10 ps -> 50 ps -> 100 ps",
        "- Any ligand-system stops at the first failing or unsafe stage.",
        "- Other ligand-systems may continue independently.",
        "",
        "## Final Status Table",
        "",
        "| ligand_id | resname | system | final_stage | completed_10ps | completed_50ps | completed_100ps | stop_reason | error |",
        "|---|---|---|---|---:|---:|---:|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['ligand_id']} | {row['resname']} | {row['system']} | {row['final_stage']} | "
            f"{row['completed_10ps']} | {row['completed_50ps']} | {row['completed_100ps']} | "
            f"{row['stop_reason']} | {row['error']} |"
        )
    STATUS_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    rows = []
    for ligand_id, resname in CANDIDATES:
        for system in SYSTEMS:
            paths = result_paths(ligand_id, system)
            payloads = {name: load_json_if_exists(path) for name, path in paths.items()}
            rows.append(summarize_system_status(ligand_id, resname, system, payloads))
    write_csv(rows)
    write_report(rows)
    print(f"status_csv: {STATUS_CSV}")
    print(f"status_md: {STATUS_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
