from __future__ import annotations

import csv
import subprocess
from pathlib import Path


TRKB_ROOT = Path(r"C:\TRKB_WP2")
MD_ROOT = TRKB_ROOT / "ligand_bound_MD"

CGENFF_ROOT = MD_ROOT / "cgenff_parameterization"
UPLOAD_MOL2_DIR = CGENFF_ROOT / "upload_mol2"
RETURNED_STR_DIR = CGENFF_ROOT / "returned_str"
LOGS_DIR = CGENFF_ROOT / "logs"
REPORTS_DIR = CGENFF_ROOT / "reports"

SELECTED_INPUTS_DIR = MD_ROOT / "selected_inputs"
CHECKLIST_CSV = MD_ROOT / "parameterization_needed" / "ligand_parameterization_checklist.csv"
MAPPING_CSV = REPORTS_DIR / "ligand_resname_mapping.csv"
UPLOAD_README = CGENFF_ROOT / "README_upload_to_cgenff.md"

TOP10_LIGANDS = ["8.1", "14.2", "19.1", "9.2", "12.2", "2.3", "4.3", "6.2", "8.3", "17.2"]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def resname_for_rank(rank: int) -> str:
    return f"L{rank:03d}"


def ligand_plan() -> list[dict[str, str]]:
    return [
        {
            "rank": str(rank),
            "ligand_id": ligand_id,
            "resname": resname_for_rank(rank),
        }
        for rank, ligand_id in enumerate(TOP10_LIGANDS, start=1)
    ]


def find_source_sdf(ligand_id: str) -> Path | None:
    preferred = SELECTED_INPUTS_DIR / ligand_id / "20chol" / f"{ligand_id}.sdf"
    if preferred.exists():
        return preferred

    for system in ["20chol", "40chol"]:
        folder = SELECTED_INPUTS_DIR / ligand_id / system
        if not folder.exists():
            continue
        exact = folder / f"{ligand_id}.sdf"
        if exact.exists():
            return exact
        matches = sorted(folder.glob("*.sdf"))
        if matches:
            return matches[0]
    return None


def check_mol2(path: Path) -> tuple[bool, bool, bool]:
    if not path.exists() or path.stat().st_size == 0:
        return False, False, False

    section = None
    has_atoms = False
    has_bonds = False
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        upper = line.upper()
        if upper.startswith("@<TRIPOS>"):
            section = upper
            continue
        if not line:
            continue
        if section == "@<TRIPOS>ATOM":
            has_atoms = True
        elif section == "@<TRIPOS>BOND":
            has_bonds = True

    return True, has_atoms, has_bonds


def rewrite_mol2_resname(path: Path, resname: str) -> None:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    rewritten: list[str] = []
    section = None
    replace_molecule_name = False

    for line in lines:
        stripped = line.strip()
        upper = stripped.upper()
        if upper.startswith("@<TRIPOS>"):
            section = upper
            replace_molecule_name = section == "@<TRIPOS>MOLECULE"
            rewritten.append(line)
            continue

        if replace_molecule_name and stripped:
            rewritten.append(resname)
            replace_molecule_name = False
            continue

        if section == "@<TRIPOS>ATOM" and stripped:
            parts = line.split()
            if len(parts) >= 8:
                parts[7] = resname
                rewritten.append(" ".join(parts))
                continue

        if section == "@<TRIPOS>SUBSTRUCTURE" and stripped:
            parts = line.split()
            if len(parts) >= 2:
                parts[1] = resname
                rewritten.append(" ".join(parts))
                continue

        rewritten.append(line)

    path.write_text("\n".join(rewritten) + "\n", encoding="utf-8")


def convert_sdf_to_mol2(source_sdf: Path, output_mol2: Path, resname: str, log_path: Path) -> bool:
    output_mol2.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    result = subprocess.run(
        ["obabel", str(source_sdf), "-O", str(output_mol2), "-h"],
        cwd=str(MD_ROOT),
        capture_output=True,
        text=True,
    )
    log_path.write_text(
        "\n".join(
            [
                f"command: obabel {source_sdf} -O {output_mol2} -h",
                f"returncode: {result.returncode}",
                "",
                "stdout:",
                result.stdout,
                "",
                "stderr:",
                result.stderr,
                "",
            ]
        ),
        encoding="utf-8",
    )

    if result.returncode != 0 or not output_mol2.exists() or output_mol2.stat().st_size == 0:
        return False

    rewrite_mol2_resname(output_mol2, resname)
    return True


def write_upload_readme(mapping_rows: list[dict[str, str]]) -> None:
    lines = [
        "# CGenFF / ParamChem Upload Package",
        "",
        "Submit the `.mol2` files in `upload_mol2` to CGenFF/ParamChem.",
        "",
        "Each ligand needs a returned `.str` file. Recommended returned file names:",
        "",
        "| Rank | Ligand ID | Residue name | Upload mol2 | Expected returned str |",
        "|---:|---|---|---|---|",
    ]
    for row in mapping_rows:
        ligand_id = row["ligand_id"]
        resname = row["resname"]
        mol2_name = Path(row["output_mol2"]).name if row["output_mol2"] else ""
        str_name = f"{resname}_{ligand_id}.str"
        lines.append(f"| {row['rank']} | {ligand_id} | {resname} | `{mol2_name}` | `{str_name}` |")

    lines.extend(
        [
            "",
            "Place returned `.str` files in:",
            "",
            f"`{RETURNED_STR_DIR}`",
            "",
            "Do not manually fabricate `.str` files.",
            "",
            "After ParamChem/CGenFF returns the `.str` files, inspect CGenFF penalties carefully.",
            "Pay special attention to high-penalty bond, angle, and dihedral parameters before any MD setup.",
        ]
    )
    UPLOAD_README.write_text("\n".join(lines) + "\n", encoding="utf-8")


def update_checklist(mapping_rows: list[dict[str, str]]) -> None:
    rows = read_csv(CHECKLIST_CSV)
    existing_fields = list(rows[0].keys())
    new_fields = [
        "resname",
        "mol2_ready",
        "mol2_path",
        "str_ready",
        "str_path",
        "cgenff_penalty_checked",
        "parameter_status",
    ]
    fieldnames = existing_fields[:]
    for field in new_fields:
        if field not in fieldnames:
            fieldnames.append(field)

    mapping_by_ligand = {row["ligand_id"]: row for row in mapping_rows}
    for row in rows:
        ligand_id = row.get("ligand_id", "")
        if ligand_id not in mapping_by_ligand:
            continue
        mapped = mapping_by_ligand[ligand_id]
        str_path = RETURNED_STR_DIR / f"{mapped['resname']}_{ligand_id}.str"
        row["resname"] = mapped["resname"]
        row["mol2_ready"] = "1" if mapped["mol2_found"] == "1" and mapped["mol2_has_atoms"] == "1" and mapped["mol2_has_bonds"] == "1" else "0"
        row["mol2_path"] = mapped["output_mol2"]
        row["str_ready"] = "1" if str_path.exists() else "0"
        row["str_path"] = str(str_path)
        row["cgenff_penalty_checked"] = "0"
        row["parameter_status"] = "pending_cgenff" if row["mol2_ready"] == "1" and row["str_ready"] == "0" else "pending"

    write_csv(CHECKLIST_CSV, rows, fieldnames)


def main() -> int:
    for directory in [CGENFF_ROOT, UPLOAD_MOL2_DIR, RETURNED_STR_DIR, LOGS_DIR, REPORTS_DIR]:
        directory.mkdir(parents=True, exist_ok=True)

    mapping_rows: list[dict[str, str]] = []
    failed: list[str] = []

    for item in ligand_plan():
        rank = item["rank"]
        ligand_id = item["ligand_id"]
        resname = item["resname"]
        source_sdf = find_source_sdf(ligand_id)
        output_mol2 = UPLOAD_MOL2_DIR / f"{resname}_{ligand_id}.mol2"
        log_path = LOGS_DIR / f"{resname}_{ligand_id}_obabel.log"
        notes = []

        if source_sdf is None:
            found, has_atoms, has_bonds = False, False, False
            notes.append("source SDF not found")
            failed.append(ligand_id)
        else:
            converted = convert_sdf_to_mol2(source_sdf, output_mol2, resname, log_path)
            found, has_atoms, has_bonds = check_mol2(output_mol2)
            if not converted:
                notes.append("Open Babel conversion failed")
            if found and not has_atoms:
                notes.append("mol2 missing ATOM records")
            if found and not has_bonds:
                notes.append("mol2 missing BOND records")
            if not (converted and found and has_atoms and has_bonds):
                failed.append(ligand_id)
            if converted and found and has_atoms and has_bonds:
                notes.append("ready for CGenFF upload; inspect returned penalties")

        mapping_rows.append(
            {
                "rank": rank,
                "ligand_id": ligand_id,
                "resname": resname,
                "source_sdf": str(source_sdf) if source_sdf else "",
                "output_mol2": str(output_mol2),
                "mol2_found": str(int(found)),
                "mol2_has_atoms": str(int(has_atoms)),
                "mol2_has_bonds": str(int(has_bonds)),
                "notes": "; ".join(notes),
            }
        )

    mapping_fields = [
        "rank",
        "ligand_id",
        "resname",
        "source_sdf",
        "output_mol2",
        "mol2_found",
        "mol2_has_atoms",
        "mol2_has_bonds",
        "notes",
    ]
    write_csv(MAPPING_CSV, mapping_rows, mapping_fields)
    write_upload_readme(mapping_rows)
    update_checklist(mapping_rows)

    success_count = sum(
        1
        for row in mapping_rows
        if row["mol2_found"] == "1" and row["mol2_has_atoms"] == "1" and row["mol2_has_bonds"] == "1"
    )

    print("Top10 ligand to residue name mapping:")
    for row in mapping_rows:
        print(row["rank"], row["ligand_id"], "->", row["resname"], row["output_mol2"])

    print("\nSuccessfully generated mol2 count:", success_count)
    print("Failed ligands:")
    if failed:
        for ligand_id in sorted(set(failed), key=failed.index):
            print(ligand_id)
    else:
        print("NONE")

    print("\nupload_mol2 folder:")
    print(UPLOAD_MOL2_DIR)

    print("\nREADME_upload_to_cgenff.md:")
    print(UPLOAD_README)

    print("\nUpdated checklist:")
    print(CHECKLIST_CSV)

    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
