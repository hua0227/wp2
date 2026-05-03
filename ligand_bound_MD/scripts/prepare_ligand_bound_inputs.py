from __future__ import annotations

import csv
import shutil
from pathlib import Path


TRKB_ROOT = Path(r"C:\TRKB_WP2")
DOCKING_ROOT = TRKB_ROOT / "docking"
MD_ROOT = TRKB_ROOT / "ligand_bound_MD"

INTEGRATED_CANDIDATE_SUMMARY = DOCKING_ROOT / "integrated_candidate_summary.csv"

SELECTED_INPUTS_DIR = MD_ROOT / "selected_inputs"
PARAMETERIZATION_DIR = MD_ROOT / "parameterization_needed"
SYSTEMS_20CHOL_DIR = MD_ROOT / "systems_20chol"
SYSTEMS_40CHOL_DIR = MD_ROOT / "systems_40chol"
SCRIPTS_DIR = MD_ROOT / "scripts"
NOTES_DIR = MD_ROOT / "notes"

BEST_POSE_PDBQT_DIR = DOCKING_ROOT / "best_pose_review_pdbqt"
BEST_POSE_PDB_DIR = DOCKING_ROOT / "best_pose_review_pdb"
PYMOL_PACKAGES_DIR = DOCKING_ROOT / "pymol_packages"

LIGAND_SDF_DIRS = [
    DOCKING_ROOT / "ligands" / "sdf_3d",
    Path(r"C:\Users\14566\Desktop\大创\2026\小分子\sdf_3d"),
]
LIGAND_PDBQT_DIRS = [
    DOCKING_ROOT / "ligands" / "pdbqt",
    Path(r"C:\Users\14566\Desktop\大创\2026\小分子\pdbqt"),
]

OUTPUT_TOP10 = MD_ROOT / "selected_candidates_top10.csv"
OUTPUT_TOP6 = MD_ROOT / "selected_candidates_top6.csv"
OUTPUT_TOP3 = MD_ROOT / "selected_candidates_top3.csv"
CHECKLIST_CSV = PARAMETERIZATION_DIR / "ligand_parameterization_checklist.csv"
PLAN_MD = NOTES_DIR / "next_step_ligand_bound_MD_plan.md"
WORKFLOW_MD = SCRIPTS_DIR / "README_ligand_bound_workflow.md"

SYSTEMS = ["20chol", "40chol"]


def to_float(value: str | None, default: float = 0.0) -> float:
    try:
        if value is None or value == "":
            return default
        return float(value)
    except ValueError:
        return default


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def find_exact_ligand_file(ligand_id: str, directories: list[Path], suffix: str) -> Path | None:
    target_name = f"{ligand_id}{suffix}".lower()
    for directory in directories:
        if not directory.exists():
            continue
        for item in directory.iterdir():
            if item.is_file() and item.name.lower() == target_name:
                return item
    return None


def find_best_pose(ligand_id: str, system: str, mode: str, suffix: str) -> Path | None:
    directory = BEST_POSE_PDBQT_DIR if suffix == ".pdbqt" else BEST_POSE_PDB_DIR
    if not directory.exists():
        return None

    exact = directory / f"candidate_{ligand_id}_{system}_mode{mode}{suffix}"
    if exact.exists():
        return exact

    matches = sorted(directory.glob(f"candidate_{ligand_id}_{system}_mode*{suffix}"))
    return matches[0] if matches else None


def copy_if_found(source: Path | None, destination_dir: Path, copied: list[str]) -> bool:
    if source is None or not source.exists():
        return False
    destination = destination_dir / source.name
    shutil.copy2(source, destination)
    copied.append(str(destination))
    return True


def pymol_open_path(ligand_id: str, system: str) -> Path:
    return PYMOL_PACKAGES_DIR / system / f"{ligand_id}-{system}" / f"open_{ligand_id}-{system}.pml"


def write_selected_readme(
    path: Path,
    ligand_id: str,
    system: str,
    row: dict[str, str],
    copied_files: list[str],
    missing_items: list[str],
    pml_path: Path,
) -> None:
    path.write_text(
        "\n".join(
            [
                f"Ligand-bound MD input staging for candidate {ligand_id} / {system}",
                "",
                "This directory contains copied docking-derived ligand inputs only.",
                "Do not run ligand-bound MD until CHARMM-compatible ligand parameters are prepared.",
                "",
                f"final_score_100: {row.get('final_score_100', '')}",
                f"priority_class: {row.get('priority_class', '')}",
                f"mode_20chol: {row.get('mode_20chol', '')}",
                f"mode_40chol: {row.get('mode_40chol', '')}",
                f"best_system: {row.get('best_system', '')}",
                f"best_mode: {row.get('best_mode', '')}",
                "",
                f"PyMOL package open script: {pml_path}",
                f"PyMOL package exists: {pml_path.exists()}",
                "",
                "Copied files:",
                *(f"- {item}" for item in copied_files),
                "",
                "Missing files/items:",
                *(f"- {item}" for item in missing_items),
                "",
            ]
        ),
        encoding="utf-8",
    )


def write_plan_doc() -> None:
    PLAN_MD.write_text(
        """# Ligand-Bound MD Preparation Plan

## Current Stage Goal

Prepare selected ligand-bound MD input folders for the top-ranked TRKB-TMD docking candidates without running ligand-bound MD and without modifying the original apo OpenMM systems.

## Completed Upstream Work

- Apo TRKB-TMD short OpenMM MD has been completed for 20% CHOL and 40% CHOL systems.
- AutoDock Vina docking has been completed for candidate ligands and reference ligands.
- All-mode contact analysis has been completed.
- Integrated scoring has been generated from docking affinity, TYR13/VAL17/SER20 contacts, two-chain proximity, pose mode reliability, 20%/40% CHOL robustness, and reference anchors.
- PyMOL packages are available for candidate and reference pose inspection.

## Scope of This Stage

This stage only prepares ligand-bound MD input folders. It does not run ligand-bound MD, does not regenerate receptors, and does not modify the original apo system directories.

## Ligand Parameterization Requirement

Each ligand needs CGenFF/CHARMM-compatible parameters before OpenMM MD. At minimum, each selected ligand should have:

- A ligand `.mol2` file with checked atom names and chemistry.
- A CGenFF/CHARMM `.str` stream file containing ligand topology and parameters.

Do not run OpenMM ligand-bound MD directly without ligand force-field parameters. Missing parameter files should remain marked as `pending` in the checklist.

## Later Formal Ligand-Bound MD Workflow

1. Ligand parameterization.
2. Merge the docked ligand pose back into the 20%/40% CHOL apo membrane system.
3. Check topology, residue name, atom naming, protonation, charge, and coordinate consistency.
4. Run an OpenMM read/load test.
5. Run energy minimization.
6. Run a 5000-step short MD sanity check.
7. Analyze whether the ligand remains near TYR13/VAL17/SER20.
""",
        encoding="utf-8",
    )


def write_workflow_doc() -> None:
    WORKFLOW_MD.write_text(
        """# Ligand-Bound MD Workflow Template

This is a workflow note only. It does not run MD automatically.

## 1. Select Ligands

Use `selected_candidates_top10.csv`, `selected_candidates_top6.csv`, or `selected_candidates_top3.csv` as the selection source. The current staging folders are under `selected_inputs/<ligand_id>/<system>`.

## 2. Prepare Ligand Parameters

For each ligand, prepare CHARMM-compatible ligand parameters before any MD run:

- Confirm the source SDF/PDBQT corresponds to the intended ligand.
- Generate or curate a CGenFF-compatible `.mol2`.
- Generate the matching CGenFF/CHARMM `.str`.
- Keep parameterization status in `parameterization_needed/ligand_parameterization_checklist.csv`.

## 3. Build Ligand-Bound Systems

For each ligand and system:

- Start from a copied working system, not the original apo OpenMM directory.
- Merge the docked ligand pose into the apo 20% or 40% CHOL membrane coordinates.
- Preserve receptor, membrane, solvent, ions, and ligand coordinate consistency.
- Check residue names and atom names against the ligand `.str` file.

## 4. Pre-MD Validation

- Confirm OpenMM can read the combined topology and coordinates.
- Check total charge, ligand atom count, and missing parameters.
- Visually inspect the merged system in PyMOL/VMD.

## 5. Short Sanity Run

After parameters and topology pass validation:

- Energy minimization.
- 5000-step short MD.
- Confirm ligand remains near TYR13/VAL17/SER20 and does not create severe clashes.

## 6. Longer Production MD

Only consider longer ligand-bound MD after the short sanity run is stable and the ligand parameters have been reviewed.
""",
        encoding="utf-8",
    )


def selected_rows(rows: list[dict[str, str]], n: int) -> list[dict[str, str]]:
    sorted_rows = sorted(rows, key=lambda r: to_float(r.get("final_score_100")), reverse=True)
    result = []
    for rank, row in enumerate(sorted_rows[:n], start=1):
        copied = dict(row)
        copied["rank"] = str(rank)
        result.append(copied)
    return result


def main() -> int:
    if not INTEGRATED_CANDIDATE_SUMMARY.exists():
        raise FileNotFoundError(INTEGRATED_CANDIDATE_SUMMARY)

    for directory in [
        MD_ROOT,
        SELECTED_INPUTS_DIR,
        PARAMETERIZATION_DIR,
        SYSTEMS_20CHOL_DIR,
        SYSTEMS_40CHOL_DIR,
        SCRIPTS_DIR,
        NOTES_DIR,
    ]:
        directory.mkdir(parents=True, exist_ok=True)

    all_candidates = read_csv(INTEGRATED_CANDIDATE_SUMMARY)
    top10 = selected_rows(all_candidates, 10)
    top6 = selected_rows(all_candidates, 6)
    top3 = selected_rows(all_candidates, 3)

    selected_fieldnames = ["rank"] + [name for name in all_candidates[0].keys() if name != "rank"]
    write_csv(OUTPUT_TOP10, top10, selected_fieldnames)
    write_csv(OUTPUT_TOP6, top6, selected_fieldnames)
    write_csv(OUTPUT_TOP3, top3, selected_fieldnames)

    checklist_rows: list[dict[str, str]] = []
    missing_records: list[str] = []

    for row in top10:
        ligand_id = row["ligand_id"]
        ligand_root = SELECTED_INPUTS_DIR / ligand_id
        ligand_root.mkdir(parents=True, exist_ok=True)

        source_sdf = find_exact_ligand_file(ligand_id, LIGAND_SDF_DIRS, ".sdf")
        source_pdbqt = find_exact_ligand_file(ligand_id, LIGAND_PDBQT_DIRS, ".pdbqt")

        pose_found = {}

        for system in SYSTEMS:
            system_dir = ligand_root / system
            system_dir.mkdir(parents=True, exist_ok=True)

            copied_files: list[str] = []
            missing_items: list[str] = []
            mode = row.get(f"mode_{system}", "")

            pose_pdbqt = find_best_pose(ligand_id, system, mode, ".pdbqt")
            pose_pdb = find_best_pose(ligand_id, system, mode, ".pdb")

            pose_found[(system, "pdbqt")] = copy_if_found(pose_pdbqt, system_dir, copied_files)
            pose_found[(system, "pdb")] = copy_if_found(pose_pdb, system_dir, copied_files)

            if not pose_found[(system, "pdbqt")]:
                item = f"{ligand_id} {system} best pose PDBQT"
                missing_items.append(item)
                missing_records.append(item)
            if not pose_found[(system, "pdb")]:
                item = f"{ligand_id} {system} best pose PDB"
                missing_items.append(item)
                missing_records.append(item)

            sdf_ok = copy_if_found(source_sdf, system_dir, copied_files)
            pdbqt_ok = copy_if_found(source_pdbqt, system_dir, copied_files)

            if not sdf_ok:
                item = f"{ligand_id} source SDF"
                missing_items.append(item)
            if not pdbqt_ok:
                item = f"{ligand_id} source PDBQT"
                missing_items.append(item)

            pml_path = pymol_open_path(ligand_id, system)
            if not pml_path.exists():
                item = f"{ligand_id} {system} PyMOL open_*.pml"
                missing_items.append(item)
                missing_records.append(item)

            write_selected_readme(
                system_dir / "README.txt",
                ligand_id,
                system,
                row,
                copied_files,
                missing_items,
                pml_path,
            )

        checklist_rows.append(
            {
                "ligand_id": ligand_id,
                "rank": row["rank"],
                "final_score_100": row.get("final_score_100", ""),
                "best_system": row.get("best_system", ""),
                "mode_20chol": row.get("mode_20chol", ""),
                "mode_40chol": row.get("mode_40chol", ""),
                "sdf_found": int(source_sdf is not None),
                "pdbqt_found": int(source_pdbqt is not None),
                "best_pose_20chol_pdb_found": int(pose_found.get(("20chol", "pdb"), False)),
                "best_pose_40chol_pdb_found": int(pose_found.get(("40chol", "pdb"), False)),
                "required_cgenff_mol2": f"{ligand_id}.mol2",
                "required_cgenff_str": f"{ligand_id}.str",
                "parameter_status": "pending",
                "notes": "Awaiting CGenFF/CHARMM-compatible ligand parameters.",
            }
        )

        if source_sdf is None:
            missing_records.append(f"{ligand_id} source SDF")
        if source_pdbqt is None:
            missing_records.append(f"{ligand_id} source PDBQT")

    checklist_fields = [
        "ligand_id",
        "rank",
        "final_score_100",
        "best_system",
        "mode_20chol",
        "mode_40chol",
        "sdf_found",
        "pdbqt_found",
        "best_pose_20chol_pdb_found",
        "best_pose_40chol_pdb_found",
        "required_cgenff_mol2",
        "required_cgenff_str",
        "parameter_status",
        "notes",
    ]
    write_csv(CHECKLIST_CSV, checklist_rows, checklist_fields)

    write_plan_doc()
    write_workflow_doc()

    print("Top10 candidates:")
    for row in top10:
        print(
            row["rank"],
            row["ligand_id"],
            "final_score_100=", row.get("final_score_100", ""),
            "best_system=", row.get("best_system", ""),
        )

    print("\nSelected input directories:")
    for row in top10:
        ligand_id = row["ligand_id"]
        for system in SYSTEMS:
            print(SELECTED_INPUTS_DIR / ligand_id / system)

    print("\nMissing SDF/PDBQT/pose/PyMOL files:")
    if missing_records:
        for item in sorted(set(missing_records)):
            print(item)
    else:
        print("None")

    print("\nParameterization checklist:")
    print(CHECKLIST_CSV)

    print("\nNext-step plan:")
    print(PLAN_MD)

    print("\nGenerated files:")
    print(OUTPUT_TOP10)
    print(OUTPUT_TOP6)
    print(OUTPUT_TOP3)
    print(WORKFLOW_MD)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
