from __future__ import annotations

import csv
import math
import shutil
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any


PILOT_ROOT = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol")
CLASH_ROOT = PILOT_ROOT / "clash_resolution"
REPORTS_DIR = CLASH_ROOT / "reports"
SCRIPTS_DIR = CLASH_ROOT / "scripts"
OUTPUTS_DIR = CLASH_ROOT / "outputs"

BOUND_PSF = PILOT_ROOT / "outputs" / "TRKB_20chol_L002_14.2_bound.psf"
BOUND_PDB = PILOT_ROOT / "outputs" / "TRKB_20chol_L002_14.2_bound.pdb"
CLEANED_PSF = OUTPUTS_DIR / "TRKB_20chol_L002_14.2_bound_cleaned.psf"
CLEANED_PDB = OUTPUTS_DIR / "TRKB_20chol_L002_14.2_bound_cleaned.pdb"

CSV_OUT = REPORTS_DIR / "bound_clash_classification.csv"
MD_OUT = REPORTS_DIR / "bound_clash_classification.md"
CLEAN_TCL = SCRIPTS_DIR / "make_cleaned_bound_system.tcl"

SEVERE_CUTOFF_A = 1.0
CLOSE_CUTOFF_A = 1.5

PROTEIN_RESNAMES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "ASH",
    "CYS",
    "CYM",
    "CYX",
    "GLN",
    "GLU",
    "GLH",
    "GLY",
    "HIS",
    "HSD",
    "HSE",
    "HSP",
    "ILE",
    "LEU",
    "LYS",
    "LYN",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}
LIPID_RESNAMES = {"POP", "POPC", "POPE", "POPS", "POPG", "DOPC", "DPPC", "DMPC"}
CHOLESTEROL_RESNAMES = {"CHL", "CHL1", "CHOL", "CLR"}
WATER_RESNAMES = {"TIP3", "HOH", "WAT", "SOL"}
ION_RESNAMES = {"SOD", "CLA", "POT", "NA", "CL", "K", "CAL", "MG", "ZN"}


@dataclass(frozen=True)
class PdbAtom:
    serial: int
    name: str
    resname: str
    segid: str
    resid: str
    element: str
    x: float
    y: float
    z: float


@dataclass(frozen=True)
class Contact:
    contact_type: str
    distance: float
    ligand_atom: PdbAtom
    other_atom: PdbAtom
    other_category: str


def normalize_element(raw: str, atom_name: str = "") -> str:
    token = (raw or atom_name).strip()
    letters = "".join(ch for ch in token if ch.isalpha())
    if not letters:
        return ""
    upper = letters.upper()
    if upper.startswith("CL"):
        return "Cl"
    if upper.startswith("BR"):
        return "Br"
    if upper.startswith("NA"):
        return "Na"
    return upper[0]


def parse_pdb_atom_line(line: str) -> PdbAtom | None:
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None
    try:
        serial = int(line[6:11])
        name = line[12:16].strip()
        resname = line[17:21].strip()
        resid = line[22:26].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        return None
    segid = line[72:76].strip() if len(line) >= 76 else ""
    element = normalize_element(line[76:78].strip() if len(line) >= 78 else "", name)
    return PdbAtom(serial, name, resname, segid, resid, element, x, y, z)


def read_pdb_atoms(path: Path) -> list[PdbAtom]:
    atoms: list[PdbAtom] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        atom = parse_pdb_atom_line(line)
        if atom is not None:
            atoms.append(atom)
    return atoms


def is_ligand(atom: PdbAtom) -> bool:
    return atom.segid == "LIG" or atom.resname == "L002"


def is_heavy(atom: PdbAtom) -> bool:
    return atom.element != "H"


def classify_atom(atom: PdbAtom) -> str:
    resname = atom.resname.upper()
    if resname in PROTEIN_RESNAMES or atom.segid.upper().startswith("PRO"):
        return "protein"
    if resname in CHOLESTEROL_RESNAMES:
        return "cholesterol"
    if resname in LIPID_RESNAMES or atom.segid.upper() == "MEMB":
        return "lipid"
    if resname in WATER_RESNAMES or atom.segid.upper().startswith("WT"):
        return "water"
    if resname in ION_RESNAMES or atom.segid.upper().startswith("ION"):
        return "ion"
    return "other/unknown"


def find_contacts(atoms: list[PdbAtom]) -> list[Contact]:
    ligand_heavy = [atom for atom in atoms if is_ligand(atom) and is_heavy(atom)]
    other_heavy = [atom for atom in atoms if not is_ligand(atom) and is_heavy(atom)]
    close_sq = CLOSE_CUTOFF_A * CLOSE_CUTOFF_A
    severe_sq = SEVERE_CUTOFF_A * SEVERE_CUTOFF_A
    contacts: list[Contact] = []

    for ligand_atom in ligand_heavy:
        for other_atom in other_heavy:
            dx = ligand_atom.x - other_atom.x
            dy = ligand_atom.y - other_atom.y
            dz = ligand_atom.z - other_atom.z
            d2 = dx * dx + dy * dy + dz * dz
            if d2 < close_sq:
                contact_type = "severe_clash" if d2 < severe_sq else "close_contact"
                contacts.append(
                    Contact(
                        contact_type=contact_type,
                        distance=math.sqrt(d2),
                        ligand_atom=ligand_atom,
                        other_atom=other_atom,
                        other_category=classify_atom(other_atom),
                    )
                )
    contacts.sort(key=lambda item: item.distance)
    return contacts


def deletable_water_ion_residues(contacts: list[Contact]) -> list[tuple[str, str, str, str]]:
    deletions: set[tuple[str, str, str, str]] = set()
    for contact in contacts:
        if contact.contact_type != "severe_clash":
            continue
        if contact.other_category not in {"water", "ion"}:
            continue
        atom = contact.other_atom
        deletions.add((atom.segid, atom.resid, atom.resname, contact.other_category))
    return sorted(deletions)


def contact_rows(contacts: list[Contact]) -> list[dict[str, Any]]:
    return [
        {
            "contact_type": contact.contact_type,
            "distance_A": f"{contact.distance:.3f}",
            "category": contact.other_category,
            "ligand_serial": contact.ligand_atom.serial,
            "ligand_atom_name": contact.ligand_atom.name,
            "ligand_resname": contact.ligand_atom.resname,
            "ligand_segid": contact.ligand_atom.segid,
            "ligand_resid": contact.ligand_atom.resid,
            "other_serial": contact.other_atom.serial,
            "other_atom_name": contact.other_atom.name,
            "other_resname": contact.other_atom.resname,
            "other_segid": contact.other_atom.segid,
            "other_resid": contact.other_atom.resid,
            "other_element": contact.other_atom.element,
        }
        for contact in contacts
    ]


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames = [
        "contact_type",
        "distance_A",
        "category",
        "ligand_serial",
        "ligand_atom_name",
        "ligand_resname",
        "ligand_segid",
        "ligand_resid",
        "other_serial",
        "other_atom_name",
        "other_resname",
        "other_segid",
        "other_resid",
        "other_element",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def category_table(contacts: list[Contact]) -> list[str]:
    counts = Counter(contact.other_category for contact in contacts)
    categories = ["protein", "lipid", "cholesterol", "water", "ion", "other/unknown"]
    lines = ["| Category | Count |", "|---|---:|"]
    for category in categories:
        lines.append(f"| {category} | {counts.get(category, 0)} |")
    return lines


def write_markdown(path: Path, atoms: list[PdbAtom], contacts: list[Contact], deletions: list[tuple[str, str, str, str]]) -> None:
    severe = [contact for contact in contacts if contact.contact_type == "severe_clash"]
    close = [contact for contact in contacts if contact.distance < CLOSE_CUTOFF_A]
    ligand_heavy_count = sum(1 for atom in atoms if is_ligand(atom) and is_heavy(atom))
    other_heavy_count = sum(1 for atom in atoms if not is_ligand(atom) and is_heavy(atom))
    lines = [
        "# Bound-System Clash Classification",
        "",
        "## Scope",
        "",
        "Pilot: ligand 14.2 / L002 / 20chol bound PSF/PDB. This report classifies ligand heavy atom contacts against non-ligand heavy atoms.",
        "",
        "## Summary",
        "",
        f"- Bound PDB: `{BOUND_PDB}`",
        f"- Ligand heavy atoms checked: {ligand_heavy_count}",
        f"- Non-ligand heavy atoms checked: {other_heavy_count}",
        f"- Severe clashes < {SEVERE_CUTOFF_A:.1f} A: {len(severe)}",
        f"- Close contacts < {CLOSE_CUTOFF_A:.1f} A: {len(close)}",
        "",
        "Severe clashes are included in the close-contact count.",
        "",
        "## Severe Clash Categories",
        "",
        *category_table(severe),
        "",
        "## Close Contact Categories",
        "",
        *category_table(close),
        "",
        "## Water/Ion Deletion Decision",
        "",
    ]
    if deletions:
        lines.append("The cleaned-system Tcl will delete these severe-clashing water/ion residues as whole residues:")
        for segid, resid, resname, category in deletions:
            lines.append(f"- {category}: segid {segid}, resid {resid}, resname {resname}")
    else:
        lines.append("No severe-clashing water or ion residues were found. No atoms will be deleted; cleaned outputs will be pilot-local copies of the bound PSF/PDB.")

    lines.extend(
        [
            "",
            "## Closest Contacts",
            "",
            "| Type | Distance A | Category | Ligand atom | Other atom | Other residue |",
            "|---|---:|---|---|---|---|",
        ]
    )
    for contact in contacts[:50]:
        other = contact.other_atom
        lig = contact.ligand_atom
        lines.append(
            f"| {contact.contact_type} | {contact.distance:.3f} | {contact.other_category} | "
            f"{lig.name} {lig.segid}:{lig.resid} | {other.name} {other.serial} | "
            f"{other.resname} {other.segid}:{other.resid} |"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def tcl_quote(value: str) -> str:
    return value.replace("\\", "/").replace('"', '\\"')


def write_cleaned_system_tcl(path: Path, deletions: list[tuple[str, str, str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Pilot-local cleaned system builder for 14.2 / L002 / 20chol.",
        "# Deletes only severe-clashing whole water/ion residues when present.",
        "# Does not modify original apo files or ligand .str.",
        f"set input_psf \"{tcl_quote(str(BOUND_PSF))}\"",
        f"set input_pdb \"{tcl_quote(str(BOUND_PDB))}\"",
        f"set output_psf \"{tcl_quote(str(CLEANED_PSF))}\"",
        f"set output_pdb \"{tcl_quote(str(CLEANED_PDB))}\"",
        "",
    ]
    if not deletions:
        lines.extend(
            [
                'puts "No severe water/ion clashes found. Copying bound PSF/PDB to cleaned outputs."',
                "file copy -force $input_psf $output_psf",
                "file copy -force $input_pdb $output_pdb",
                'puts "Wrote cleaned PSF: $output_psf"',
                'puts "Wrote cleaned PDB: $output_pdb"',
                "quit",
                "",
            ]
        )
    else:
        terms = [
            f'(segid "{segid}" and resid "{resid}" and resname "{resname}")'
            for segid, resid, resname, _category in deletions
        ]
        selection = " or ".join(terms)
        lines.extend(
            [
                "mol new $input_psf type psf waitfor all",
                "mol addfile $input_pdb type pdb waitfor all",
                f"set del [atomselect top \"{selection}\"]",
                'puts "Deleting severe-clashing water/ion residues:"',
                "puts [$del get {segid resid resname name index}]",
                f"set keep [atomselect top \"not ({selection})\"]",
                "$keep writepsf $output_psf",
                "$keep writepdb $output_pdb",
                "$del delete",
                "$keep delete",
                'puts "Wrote cleaned PSF: $output_psf"',
                'puts "Wrote cleaned PDB: $output_pdb"',
                "quit",
                "",
            ]
        )
    path.write_text("\n".join(lines), encoding="utf-8")


def run_classification() -> tuple[int, int, int, list[Contact], list[tuple[str, str, str, str]]]:
    atoms = read_pdb_atoms(BOUND_PDB)
    contacts = find_contacts(atoms)
    deletions = deletable_water_ion_residues(contacts)
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(CSV_OUT, contact_rows(contacts))
    write_markdown(MD_OUT, atoms, contacts, deletions)
    write_cleaned_system_tcl(CLEAN_TCL, deletions)
    severe_count = sum(1 for contact in contacts if contact.contact_type == "severe_clash")
    close_count = sum(1 for contact in contacts if contact.distance < CLOSE_CUTOFF_A)
    return severe_count, close_count, len(deletions), contacts, deletions


def main() -> int:
    severe_count, close_count, deletion_count, _contacts, _deletions = run_classification()
    print(f"severe_clashes_lt_{SEVERE_CUTOFF_A:.1f}A: {severe_count}")
    print(f"close_contacts_lt_{CLOSE_CUTOFF_A:.1f}A: {close_count}")
    print(f"water_ion_residues_to_delete: {deletion_count}")
    print(f"classification CSV: {CSV_OUT}")
    print(f"classification Markdown: {MD_OUT}")
    print(f"cleaned-system Tcl: {CLEAN_TCL}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
