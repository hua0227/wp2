from __future__ import annotations

import csv
from pathlib import Path
from typing import Any


TRKB_ROOT = Path(r"C:\TRKB_WP2")
MD_ROOT = TRKB_ROOT / "ligand_bound_MD"
DOCKING_ROOT = TRKB_ROOT / "docking"

PREFLIGHT_ROOT = MD_ROOT / "preflight"
LOGS_DIR = PREFLIGHT_ROOT / "logs"
COMPLEX_PREVIEWS_DIR = PREFLIGHT_ROOT / "complex_previews"
LIGAND_POSE_NAMED_DIR = PREFLIGHT_ROOT / "ligand_pose_named"
REPORTS_DIR = PREFLIGHT_ROOT / "reports"

SELECTED_INPUTS_DIR = MD_ROOT / "selected_inputs"
MOL2_DIR = MD_ROOT / "cgenff_parameterization" / "upload_mol2"
STR_DIR = MD_ROOT / "cgenff_parameterization" / "returned_str"
PENALTY_SUMMARY_CSV = MD_ROOT / "cgenff_parameterization" / "reports" / "cgenff_penalty_summary.csv"
CHECKLIST_CSV = MD_ROOT / "parameterization_needed" / "ligand_parameterization_checklist.csv"

TOPOLOGY_DIR = TRKB_ROOT / "TRKB_20CHOL" / "toppar"
BASE_PARAMETER_FILES = [
    TOPOLOGY_DIR / "top_all36_prot.rtf",
    TOPOLOGY_DIR / "top_all36_lipid.rtf",
    TOPOLOGY_DIR / "top_all36_cgenff.rtf",
    TOPOLOGY_DIR / "par_all36m_prot.prm",
    TOPOLOGY_DIR / "par_all36_lipid.prm",
    TOPOLOGY_DIR / "par_all36_cgenff.prm",
    TOPOLOGY_DIR / "toppar_water_ions.str",
    TOPOLOGY_DIR / "toppar_all36_lipid_cholesterol.str",
]

APO_PDB_BY_SYSTEM = {
    "20chol": TRKB_ROOT / "TRKB_20CHOL" / "openmm_short_md_output" / "short_md_final.pdb",
    "40chol": TRKB_ROOT / "TRKB_40CHOL" / "openmm_short_md_output" / "short_md_final.pdb",
}

LIGANDS = [
    ("8.1", "L001"),
    ("14.2", "L002"),
    ("19.1", "L003"),
    ("9.2", "L004"),
    ("12.2", "L005"),
    ("2.3", "L006"),
    ("4.3", "L007"),
    ("6.2", "L008"),
    ("8.3", "L009"),
    ("17.2", "L010"),
]
SYSTEMS = ["20chol", "40chol"]

REPORT_CSV = REPORTS_DIR / "preflight_ligand_parameter_check.csv"
SUMMARY_MD = REPORTS_DIR / "preflight_summary.md"

REPORT_FIELDS = [
    "ligand_id",
    "resname",
    "system",
    "str_found",
    "str_load_status",
    "str_load_error",
    "sdf_found",
    "pdbqt_found",
    "best_pose_pdb_found",
    "best_pose_pdbqt_found",
    "mol2_found",
    "mol2_atom_count",
    "pose_pdb_atom_count",
    "pose_pdbqt_atom_count",
    "atom_count_match",
    "named_pose_generated",
    "complex_preview_generated",
    "penalty_level",
    "parameter_status",
    "notes",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def mol2_element(atom_type: str, atom_name: str) -> str:
    token = atom_type.split(".", 1)[0].strip()
    letters = "".join(ch for ch in token if ch.isalpha())
    if not letters:
        letters = "".join(ch for ch in atom_name if ch.isalpha())
    if not letters:
        return ""
    if len(letters) >= 2 and letters[:2].title() in {"Cl", "Br"}:
        return letters[:2].title()
    return letters[0].upper()


def parse_mol2_atoms(path: Path) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    if not path.exists():
        return atoms

    in_atoms = False
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        upper = line.upper()
        if upper.startswith("@<TRIPOS>"):
            in_atoms = upper == "@<TRIPOS>ATOM"
            continue
        if not in_atoms or not line:
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        atoms.append(
            {
                "index": int(parts[0]),
                "name": parts[1],
                "x": float(parts[2]),
                "y": float(parts[3]),
                "z": float(parts[4]),
                "type": parts[5],
                "element": mol2_element(parts[5], parts[1]),
            }
        )
    return atoms


def parse_pdb_atom_coords(path: Path) -> list[dict[str, Any]]:
    coords: list[dict[str, Any]] = []
    if not path.exists():
        return coords
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        element = line[76:78].strip()
        coords.append({"x": x, "y": y, "z": z, "element": element, "line": line})
    return coords


def parse_pdbqt_atom_coords(path: Path) -> list[dict[str, Any]]:
    coords: list[dict[str, Any]] = []
    if not path.exists():
        return coords
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        atom_type = line[77:].strip()
        element = mol2_element(atom_type, line[12:16].strip()) if atom_type else ""
        coords.append({"x": x, "y": y, "z": z, "element": element, "line": line})
    return coords


def format_atom_name(name: str, element: str) -> str:
    clean = name[:4]
    if len(clean) == 4:
        return clean
    if len(element) == 1:
        return f" {clean:<3}"[:4]
    return f"{clean:<4}"[:4]


def build_named_pose_pdb_lines(mol2_atoms: list[dict[str, Any]], coords: list[dict[str, Any]], resname: str) -> list[str]:
    lines: list[str] = []
    for serial, (atom, coord) in enumerate(zip(mol2_atoms, coords), start=1):
        atom_name = format_atom_name(str(atom["name"]), str(atom.get("element", "")))
        element = str(atom.get("element", "")).rjust(2)[:2]
        lines.append(
            f"HETATM{serial:5d} {atom_name} {resname:>4s} Z{1:4d}    "
            f"{coord['x']:8.3f}{coord['y']:8.3f}{coord['z']:8.3f}"
            f"{1.00:6.2f}{0.00:6.2f}          {element}"
        )
    return lines


def write_named_pose(path: Path, mol2_atoms: list[dict[str, Any]], coords: list[dict[str, Any]], resname: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = build_named_pose_pdb_lines(mol2_atoms, coords, resname)
    path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")


def write_complex_preview(apo_pdb: Path, named_pose: Path, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    apo_lines: list[str] = []
    for line in apo_pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            apo_lines.append(line)
    ligand_lines = [
        line
        for line in named_pose.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith("ATOM") or line.startswith("HETATM")
    ]
    output_path.write_text("\n".join(apo_lines + ligand_lines + ["END", ""]) , encoding="utf-8")


def find_file(folder: Path, exact_name: str | None = None, pattern: str | None = None) -> Path | None:
    if not folder.exists():
        return None
    if exact_name:
        exact = folder / exact_name
        if exact.exists():
            return exact
    if pattern:
        matches = sorted(folder.glob(pattern))
        return matches[0] if matches else None
    return None


def load_penalty_info() -> dict[str, dict[str, str]]:
    return {row["ligand_id"]: row for row in read_csv(PENALTY_SUMMARY_CSV)}


def load_checklist_info() -> dict[str, dict[str, str]]:
    return {row["ligand_id"]: row for row in read_csv(CHECKLIST_CSV)}


def import_charmm_parameter_set():
    try:
        from openmm.app import CharmmParameterSet

        return CharmmParameterSet, ""
    except Exception as first:
        try:
            from simtk.openmm.app import CharmmParameterSet

            return CharmmParameterSet, ""
        except Exception as second:
            return None, f"OpenMM import failed: {first!r}; legacy simtk import failed: {second!r}"


def check_str_loads() -> dict[str, dict[str, str]]:
    CharmmParameterSet, import_error = import_charmm_parameter_set()
    results: dict[str, dict[str, str]] = {}
    missing_base = [str(path) for path in BASE_PARAMETER_FILES if not path.exists()]

    for ligand_id, resname in LIGANDS:
        str_path = STR_DIR / f"{resname}_{ligand_id}.str"
        key = ligand_id
        if not str_path.exists():
            results[key] = {"status": "missing", "error": f"missing {str_path}"}
            continue
        if missing_base:
            results[key] = {"status": "failed", "error": "missing base files: " + "; ".join(missing_base)}
            continue
        if CharmmParameterSet is None:
            results[key] = {"status": "failed", "error": import_error}
            continue

        try:
            CharmmParameterSet(*[str(path) for path in BASE_PARAMETER_FILES + [str_path]])
            results[key] = {"status": "ok", "error": ""}
        except Exception as exc:
            results[key] = {"status": "failed", "error": repr(exc)}

    log_path = LOGS_DIR / "str_load_status.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(
        "\n".join(f"{ligand_id}\t{value['status']}\t{value['error']}" for ligand_id, value in results.items()) + "\n",
        encoding="utf-8",
    )
    return results


def system_row(
    ligand_id: str,
    resname: str,
    system: str,
    str_status: dict[str, str],
    penalty: dict[str, str],
    checklist: dict[str, str],
) -> dict[str, Any]:
    selected_dir = SELECTED_INPUTS_DIR / ligand_id / system
    sdf = find_file(selected_dir, exact_name=f"{ligand_id}.sdf", pattern="*.sdf")
    pdbqt = find_file(selected_dir, exact_name=f"{ligand_id}.pdbqt", pattern="*.pdbqt")
    best_pose_pdb = find_file(selected_dir, pattern=f"candidate_{ligand_id}_{system}_mode*.pdb")
    best_pose_pdbqt = find_file(selected_dir, pattern=f"candidate_{ligand_id}_{system}_mode*.pdbqt")
    mol2 = MOL2_DIR / f"{resname}_{ligand_id}.mol2"
    str_file = STR_DIR / f"{resname}_{ligand_id}.str"

    mol2_atoms = parse_mol2_atoms(mol2)
    pose_pdb_coords = parse_pdb_atom_coords(best_pose_pdb) if best_pose_pdb else []
    pose_pdbqt_coords = parse_pdbqt_atom_coords(best_pose_pdbqt) if best_pose_pdbqt else []
    mol2_count = len(mol2_atoms)
    pdb_count = len(pose_pdb_coords)
    pdbqt_count = len(pose_pdbqt_coords)
    atom_count_match = mol2_count > 0 and mol2_count == pdb_count and mol2_count == pdbqt_count

    notes: list[str] = []
    named_pose_generated = False
    complex_preview_generated = False

    if not selected_dir.exists():
        notes.append("selected input directory missing")
    if not atom_count_match:
        notes.append("atom count mismatch; named pose not generated")
    else:
        coords = pose_pdb_coords if pdb_count == mol2_count else pose_pdbqt_coords
        named_pose = (
            LIGAND_POSE_NAMED_DIR
            / ligand_id
            / system
            / f"{resname}_{ligand_id}_{system}_pose_named.pdb"
        )
        write_named_pose(named_pose, mol2_atoms, coords, resname)
        named_pose_generated = named_pose.exists()

        apo_pdb = APO_PDB_BY_SYSTEM[system]
        if apo_pdb.exists() and named_pose_generated:
            preview = (
                COMPLEX_PREVIEWS_DIR
                / ligand_id
                / system
                / f"complex_preview_{ligand_id}_{system}.pdb"
            )
            write_complex_preview(apo_pdb, named_pose, preview)
            complex_preview_generated = preview.exists()
        elif not apo_pdb.exists():
            notes.append(f"apo PDB missing for {system}")

    return {
        "ligand_id": ligand_id,
        "resname": resname,
        "system": system,
        "str_found": int(str_file.exists()),
        "str_load_status": str_status["status"],
        "str_load_error": str_status["error"],
        "sdf_found": int(sdf is not None and sdf.exists()),
        "pdbqt_found": int(pdbqt is not None and pdbqt.exists()),
        "best_pose_pdb_found": int(best_pose_pdb is not None and best_pose_pdb.exists()),
        "best_pose_pdbqt_found": int(best_pose_pdbqt is not None and best_pose_pdbqt.exists()),
        "mol2_found": int(mol2.exists()),
        "mol2_atom_count": mol2_count,
        "pose_pdb_atom_count": pdb_count,
        "pose_pdbqt_atom_count": pdbqt_count,
        "atom_count_match": int(atom_count_match),
        "named_pose_generated": int(named_pose_generated),
        "complex_preview_generated": int(complex_preview_generated),
        "penalty_level": penalty.get("penalty_level", ""),
        "parameter_status": checklist.get("parameter_status", ""),
        "notes": "; ".join(notes),
    }


def write_summary(rows: list[dict[str, Any]], str_results: dict[str, dict[str, str]]) -> None:
    str_ok = [ligand_id for ligand_id, value in str_results.items() if value["status"] == "ok"]
    str_failed = [ligand_id for ligand_id, value in str_results.items() if value["status"] != "ok"]
    by_ligand: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        by_ligand.setdefault(str(row["ligand_id"]), []).append(row)

    atom_ok = [
        ligand_id
        for ligand_id, ligand_rows in by_ligand.items()
        if all(int(row["atom_count_match"]) == 1 for row in ligand_rows)
    ]
    named_ok = [
        ligand_id
        for ligand_id, ligand_rows in by_ligand.items()
        if all(int(row["named_pose_generated"]) == 1 for row in ligand_rows)
    ]
    preview_ok = [
        ligand_id
        for ligand_id, ligand_rows in by_ligand.items()
        if all(int(row["complex_preview_generated"]) == 1 for row in ligand_rows)
    ]
    cautious = [
        ligand_id
        for ligand_id, ligand_rows in by_ligand.items()
        if ligand_rows[0]["penalty_level"] in {"high", "very_high", "unknown"}
    ]
    recommended_priority = ["14.2", "19.1", "2.3", "6.2", "17.2"]
    recommended = [
        ligand_id
        for ligand_id in recommended_priority
        if ligand_id in by_ligand
        and ligand_id in str_ok
        and by_ligand[ligand_id][0]["penalty_level"] == "moderate"
    ]

    lines = [
        "# Ligand-Bound MD Preflight Summary",
        "",
        "## Purpose",
        "",
        "This preflight checks ligand parameter file readability, staged input consistency, atom counts, CGenFF naming for docked poses, and visual-only complex previews before any ligand-bound MD setup.",
        "",
        "No MD, minimization, receptor regeneration, docking rerun, `.str` modification, PSF generation, or original apo system modification was performed.",
        "",
        "Complex preview PDB files are only for PyMOL visualization. They are not simulation-ready topologies and must not be used directly as proof that OpenMM ligand-bound MD setup is complete.",
        "",
        "## CharmmParameterSet Load Check",
        "",
        f"- Successful ligand `.str` loads: {len(str_ok)}",
        f"- Failed ligand `.str` loads: {', '.join(str_failed) if str_failed else 'NONE'}",
        "",
        "## Atom Count Consistency",
        "",
        f"- Ligands with atom counts consistent in both systems: {', '.join(atom_ok) if atom_ok else 'NONE'}",
        "",
        "## Named Pose Generation",
        "",
        f"- Ligands with named poses generated for both systems: {', '.join(named_ok) if named_ok else 'NONE'}",
        "",
        "## Complex Preview Generation",
        "",
        f"- Ligands with complex previews generated for both systems: {', '.join(preview_ok) if preview_ok else 'NONE'}",
        "",
        "## High-Penalty Ligands Requiring Caution",
        "",
        f"- {', '.join(cautious) if cautious else 'NONE'}",
        "",
        "High penalty does not automatically exclude a ligand, but it flags parameters for careful manual review or optimization before MD.",
        "",
        "## Recommended First-Round OpenMM Read-Test Ligands",
        "",
        f"- {', '.join(recommended) if recommended else 'NONE'}",
        "",
        "These recommendations prioritize moderate-penalty ligands requested for first-round parameter/read testing and successful `.str` loading.",
        "Because atom-count matching failed for the current docked poses, complex-level OpenMM read tests still require a corrected ligand pose/topology alignment first.",
    ]
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    for directory in [LOGS_DIR, COMPLEX_PREVIEWS_DIR, LIGAND_POSE_NAMED_DIR, REPORTS_DIR]:
        directory.mkdir(parents=True, exist_ok=True)

    str_results = check_str_loads()
    penalty_info = load_penalty_info()
    checklist_info = load_checklist_info()

    rows: list[dict[str, Any]] = []
    for ligand_id, resname in LIGANDS:
        for system in SYSTEMS:
            rows.append(
                system_row(
                    ligand_id,
                    resname,
                    system,
                    str_results.get(ligand_id, {"status": "failed", "error": "missing str load result"}),
                    penalty_info.get(ligand_id, {}),
                    checklist_info.get(ligand_id, {}),
                )
            )

    write_csv(REPORT_CSV, rows, REPORT_FIELDS)
    write_summary(rows, str_results)

    str_success = [ligand_id for ligand_id, value in str_results.items() if value["status"] == "ok"]
    str_failed = [ligand_id for ligand_id, value in str_results.items() if value["status"] != "ok"]
    mismatches = [
        f"{row['ligand_id']} {row['system']}"
        for row in rows
        if int(row["atom_count_match"]) != 1
    ]
    named_count = sum(1 for row in rows if int(row["named_pose_generated"]) == 1)
    preview_count = sum(1 for row in rows if int(row["complex_preview_generated"]) == 1)
    recommended_priority = ["14.2", "19.1", "2.3", "6.2", "17.2"]
    by_ligand: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        by_ligand.setdefault(str(row["ligand_id"]), []).append(row)
    recommended = [
        ligand_id
        for ligand_id in recommended_priority
        if ligand_id in by_ligand
        and str_results.get(ligand_id, {}).get("status") == "ok"
        and by_ligand[ligand_id][0]["penalty_level"] == "moderate"
    ]

    print("str load success count:", len(str_success))
    print("str load failed ligands:")
    if str_failed:
        for ligand_id in str_failed:
            print(ligand_id, str_results[ligand_id]["error"])
    else:
        print("NONE")

    print("\natom count mismatch list:")
    if mismatches:
        for item in mismatches:
            print(item)
    else:
        print("NONE")

    print("\nnamed pose generated count:", named_count)
    print("complex preview generated count:", preview_count)
    print("recommended first-round OpenMM read-test ligands:")
    print(", ".join(recommended) if recommended else "NONE")
    print("\npreflight CSV:")
    print(REPORT_CSV)
    print("\npreflight Markdown:")
    print(SUMMARY_MD)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
