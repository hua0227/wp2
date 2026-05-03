from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path


TRKB_ROOT = Path(r"C:\TRKB_WP2")
DOCKING_ROOT = TRKB_ROOT / "docking"
MD_ROOT = TRKB_ROOT / "ligand_bound_MD"

BEST_POSE_PER_LIGAND = DOCKING_ROOT / "best_pose_per_ligand.csv"
BEST_POSE_PDBQT_DIR = DOCKING_ROOT / "best_pose_review_pdbqt"
BEST_POSE_PDB_DIR = DOCKING_ROOT / "best_pose_review_pdb"

SELECTED_INPUTS_DIR = MD_ROOT / "selected_inputs"
CHECKLIST_CSV = MD_ROOT / "parameterization_needed" / "ligand_parameterization_checklist.csv"

LIGAND_SDF_DIRS = [
    DOCKING_ROOT / "ligands" / "sdf_3d",
    Path(r"C:\Users\14566\Desktop\大创\2026\小分子\sdf_3d"),
]
LIGAND_PDBQT_DIRS = [
    DOCKING_ROOT / "ligands" / "pdbqt",
    Path(r"C:\Users\14566\Desktop\大创\2026\小分子\pdbqt"),
]

TOP10_LIGANDS = ["8.1", "14.2", "19.1", "9.2", "12.2", "2.3", "4.3", "6.2", "8.3", "17.2"]
SYSTEMS = ["20chol", "40chol"]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def extract_pdbqt_mode(source: Path, mode: int) -> str:
    capture = False
    model_seen = False
    lines: list[str] = []
    target_model = f"MODEL {mode}"

    for line in source.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("MODEL"):
            capture = stripped == target_model
            if capture:
                model_seen = True
                lines = [line]
            continue

        if capture:
            lines.append(line)
            if stripped == "ENDMDL":
                return "\n".join(lines) + "\n"

    if model_seen:
        return "\n".join(lines) + "\n"

    raise ValueError(f"Mode {mode} not found in {source}")


def find_exact_ligand_file(ligand_id: str, directories: list[Path], suffix: str) -> Path | None:
    target_name = f"{ligand_id}{suffix}".lower()
    for directory in directories:
        if not directory.exists():
            continue
        for item in directory.iterdir():
            if item.is_file() and item.name.lower() == target_name:
                return item
    return None


def best_pose_output_paths(ligand_id: str, system: str, mode: str) -> tuple[Path, Path]:
    stem = f"candidate_{ligand_id}_{system}_mode{mode}"
    return BEST_POSE_PDBQT_DIR / f"{stem}.pdbqt", BEST_POSE_PDB_DIR / f"{stem}.pdb"


def convert_pdbqt_to_pdb(pdbqt_path: Path, pdb_path: Path) -> None:
    result = subprocess.run(
        ["obabel", "-ipdbqt", str(pdbqt_path), "-opdb", "-O", str(pdb_path)],
        cwd=str(DOCKING_ROOT),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Open Babel conversion failed for {pdbqt_path}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )


def copy_file(source: Path | None, destination_dir: Path) -> bool:
    if source is None or not source.exists():
        return False
    destination_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination_dir / source.name)
    return True


def top10_pose_rows() -> list[dict[str, str]]:
    rows = read_csv(BEST_POSE_PER_LIGAND)
    by_key = {
        (row["ligand_id"], row["system"]): row
        for row in rows
        if row.get("group") == "candidate"
        and row.get("ligand_id") in TOP10_LIGANDS
        and row.get("system") in SYSTEMS
    }

    ordered = []
    for ligand_id in TOP10_LIGANDS:
        for system in SYSTEMS:
            row = by_key.get((ligand_id, system))
            if row is None:
                raise ValueError(f"Missing best_pose_per_ligand row for {ligand_id} {system}")
            ordered.append(row)
    return ordered


def write_selected_readme(
    destination: Path,
    ligand_id: str,
    system: str,
    mode: str,
    pose_pdbqt: Path,
    pose_pdb: Path,
    source_sdf: Path | None,
    source_pdbqt: Path | None,
    missing: list[str],
) -> None:
    status = "NONE" if not missing else "\n".join(f"- {item}" for item in missing)
    destination.write_text(
        "\n".join(
            [
                f"Ligand-bound MD staged inputs for candidate {ligand_id} / {system}",
                "",
                "This directory is input staging only. Do not run ligand-bound MD until ligand parameters are ready.",
                "",
                f"selected mode: {mode}",
                f"best pose PDBQT: {pose_pdbqt}",
                f"best pose PDB: {pose_pdb}",
                f"source ligand SDF: {source_sdf or ''}",
                f"source ligand PDBQT: {source_pdbqt or ''}",
                "",
                "Missing files/items:",
                status,
                "",
            ]
        ),
        encoding="utf-8",
    )


def sync_selected_inputs(generated: dict[tuple[str, str], tuple[Path, Path]]) -> list[str]:
    missing: list[str] = []

    for ligand_id in TOP10_LIGANDS:
        source_sdf = find_exact_ligand_file(ligand_id, LIGAND_SDF_DIRS, ".sdf")
        source_pdbqt = find_exact_ligand_file(ligand_id, LIGAND_PDBQT_DIRS, ".pdbqt")

        for system in SYSTEMS:
            system_missing: list[str] = []
            destination_dir = SELECTED_INPUTS_DIR / ligand_id / system
            destination_dir.mkdir(parents=True, exist_ok=True)

            pose_pdbqt, pose_pdb = generated[(ligand_id, system)]
            if not copy_file(source_sdf, destination_dir):
                system_missing.append(f"{ligand_id} source SDF")
            if not copy_file(source_pdbqt, destination_dir):
                system_missing.append(f"{ligand_id} source PDBQT")
            if not copy_file(pose_pdbqt, destination_dir):
                system_missing.append(f"{ligand_id} {system} best pose PDBQT")
            if not copy_file(pose_pdb, destination_dir):
                system_missing.append(f"{ligand_id} {system} best pose PDB")

            mode = pose_pdbqt.stem.split("_mode")[-1]
            write_selected_readme(
                destination_dir / "README.txt",
                ligand_id,
                system,
                mode,
                pose_pdbqt,
                pose_pdb,
                source_sdf,
                source_pdbqt,
                system_missing,
            )
            missing.extend(system_missing)

    return missing


def update_checklist(generated: dict[tuple[str, str], tuple[Path, Path]]) -> None:
    rows = read_csv(CHECKLIST_CSV)
    fieldnames = list(rows[0].keys())

    for row in rows:
        ligand_id = row["ligand_id"]
        if ligand_id not in TOP10_LIGANDS:
            continue
        source_sdf = find_exact_ligand_file(ligand_id, LIGAND_SDF_DIRS, ".sdf")
        source_pdbqt = find_exact_ligand_file(ligand_id, LIGAND_PDBQT_DIRS, ".pdbqt")
        row["sdf_found"] = str(int(source_sdf is not None))
        row["pdbqt_found"] = str(int(source_pdbqt is not None))
        row["best_pose_20chol_pdb_found"] = str(int(generated[(ligand_id, "20chol")][1].exists()))
        row["best_pose_40chol_pdb_found"] = str(int(generated[(ligand_id, "40chol")][1].exists()))
        row["parameter_status"] = "pending"

    write_csv(CHECKLIST_CSV, rows, fieldnames)


def verify_outputs(generated: dict[tuple[str, str], tuple[Path, Path]]) -> list[str]:
    missing: list[str] = []
    for ligand_id in TOP10_LIGANDS:
        source_sdf = find_exact_ligand_file(ligand_id, LIGAND_SDF_DIRS, ".sdf")
        source_pdbqt = find_exact_ligand_file(ligand_id, LIGAND_PDBQT_DIRS, ".pdbqt")
        for system in SYSTEMS:
            pose_pdbqt, pose_pdb = generated[(ligand_id, system)]
            selected_dir = SELECTED_INPUTS_DIR / ligand_id / system

            checks = [
                (pose_pdbqt, f"{ligand_id} {system} best pose PDBQT"),
                (pose_pdb, f"{ligand_id} {system} best pose PDB"),
                (selected_dir / (source_sdf.name if source_sdf else f"{ligand_id}.sdf"), f"{ligand_id} {system} selected original SDF"),
                (
                    selected_dir / (source_pdbqt.name if source_pdbqt else f"{ligand_id}.pdbqt"),
                    f"{ligand_id} {system} selected original PDBQT",
                ),
                (selected_dir / pose_pdbqt.name, f"{ligand_id} {system} selected best pose PDBQT"),
                (selected_dir / pose_pdb.name, f"{ligand_id} {system} selected best pose PDB"),
                (selected_dir / "README.txt", f"{ligand_id} {system} README"),
            ]
            for path, label in checks:
                if not path.exists():
                    missing.append(label)
    return missing


def main() -> int:
    BEST_POSE_PDBQT_DIR.mkdir(parents=True, exist_ok=True)
    BEST_POSE_PDB_DIR.mkdir(parents=True, exist_ok=True)

    generated: dict[tuple[str, str], tuple[Path, Path]] = {}
    for row in top10_pose_rows():
        ligand_id = row["ligand_id"]
        system = row["system"]
        mode = row["mode"]
        pose_file = Path(row["pose_file"])
        if not pose_file.exists():
            raise FileNotFoundError(pose_file)

        output_pdbqt, output_pdb = best_pose_output_paths(ligand_id, system, mode)
        output_pdbqt.write_text(extract_pdbqt_mode(pose_file, int(mode)), encoding="utf-8")
        convert_pdbqt_to_pdb(output_pdbqt, output_pdb)
        generated[(ligand_id, system)] = (output_pdbqt, output_pdb)

    sync_missing = sync_selected_inputs(generated)
    update_checklist(generated)
    verify_missing = verify_outputs(generated)
    missing = sorted(set(sync_missing + verify_missing))

    all_best_pose_pdbqt = all(generated[(ligand_id, system)][0].exists() for ligand_id in TOP10_LIGANDS for system in SYSTEMS)
    all_best_pose_pdb = all(generated[(ligand_id, system)][1].exists() for ligand_id in TOP10_LIGANDS for system in SYSTEMS)
    all_selected_synced = not verify_missing

    print("Top10 x 2 systems all have best pose PDBQT:", all_best_pose_pdbqt)
    print("Top10 x 2 systems all have best pose PDB:", all_best_pose_pdb)
    print("selected_inputs sync complete:", all_selected_synced)
    print("\nGenerated/checked best poses:")
    for ligand_id in TOP10_LIGANDS:
        for system in SYSTEMS:
            pdbqt, pdb = generated[(ligand_id, system)]
            print(ligand_id, system, "PDBQT=", pdbqt, "PDB=", pdb)

    print("\nStill missing files:")
    if missing:
        for item in missing:
            print(item)
    else:
        print("NONE")

    print("\nUpdated checklist:")
    print(CHECKLIST_CSV)

    return 0 if not missing else 1


if __name__ == "__main__":
    raise SystemExit(main())
