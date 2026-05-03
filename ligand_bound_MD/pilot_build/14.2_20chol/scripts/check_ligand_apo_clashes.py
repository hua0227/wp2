from __future__ import annotations

import csv
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any


PILOT_ROOT = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol")
REPORTS_DIR = PILOT_ROOT / "reports"

APO_PDB = Path(r"C:\TRKB_WP2\TRKB_20CHOL\openmm_short_md_output\short_md_final.pdb")
LIGAND_PDB = Path(
    r"C:\TRKB_WP2\ligand_bound_MD\preflight\ligand_pose_allatom_named\14.2\20chol\L002_14.2_20chol_allatom_named.pdb"
)

CSV_OUT = REPORTS_DIR / "clash_precheck.csv"
MD_OUT = REPORTS_DIR / "clash_precheck.md"

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
    "HID",
    "HIE",
    "HIP",
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
WATER_RESNAMES = {"HOH", "TIP3", "WAT", "SOL"}
ION_RESNAMES = {
    "CLA",
    "CL",
    "SOD",
    "NA",
    "POT",
    "K",
    "CAL",
    "CA",
    "MG",
    "ZN",
    "CES",
}
CHOLESTEROL_RESNAMES = {"CHL", "CHOL", "CHL1", "CLR"}
LIPID_RESNAMES = {
    "POP",
    "POPC",
    "POPE",
    "POPS",
    "POPG",
    "DPPC",
    "DOPC",
    "DOPE",
    "DOPS",
    "DOPG",
    "DMPC",
    "DLPC",
    "DSPC",
}


@dataclass(frozen=True)
class PdbAtom:
    serial: int
    name: str
    resname: str
    chain: str
    resid: str
    element: str
    x: float
    y: float
    z: float
    line: str


@dataclass(frozen=True)
class Contact:
    contact_type: str
    distance: float
    ligand_atom: PdbAtom
    apo_atom: PdbAtom
    apo_category: str


def normalize_element(raw: str, atom_name: str = "") -> str:
    token = raw.strip()
    if not token:
        token = atom_name.strip()
    letters = "".join(ch for ch in token if ch.isalpha())
    if not letters:
        return ""
    upper = letters.upper()
    if upper.startswith("CL"):
        return "Cl"
    if upper.startswith("BR"):
        return "Br"
    return upper[0]


def parse_pdb_atom_line(line: str) -> PdbAtom | None:
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None
    try:
        serial = int(line[6:11])
        name = line[12:16].strip()
        resname = line[17:21].strip()
        chain = line[21].strip()
        resid = line[22:26].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        return None
    element = normalize_element(line[76:78].strip(), name)
    return PdbAtom(serial, name, resname, chain, resid, element, x, y, z, line.rstrip("\n"))


def read_pdb_atoms(path: Path) -> list[PdbAtom]:
    atoms: list[PdbAtom] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        atom = parse_pdb_atom_line(line)
        if atom is not None:
            atoms.append(atom)
    return atoms


def classify_residue(resname: str) -> str:
    clean = resname.strip().upper()
    if clean in PROTEIN_RESNAMES:
        return "protein"
    if clean in CHOLESTEROL_RESNAMES:
        return "cholesterol"
    if clean in LIPID_RESNAMES:
        return "lipid"
    if clean in WATER_RESNAMES:
        return "water"
    if clean in ION_RESNAMES:
        return "ion"
    return "unknown"


def heavy_only(atoms: list[PdbAtom]) -> list[PdbAtom]:
    return [atom for atom in atoms if atom.element != "H"]


def distance(a: PdbAtom, b: PdbAtom) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def find_contacts(
    ligand_atoms: list[PdbAtom],
    apo_atoms: list[PdbAtom],
    severe_cutoff: float = SEVERE_CUTOFF_A,
    close_cutoff: float = CLOSE_CUTOFF_A,
) -> list[Contact]:
    contacts: list[Contact] = []
    ligand_heavy = heavy_only(ligand_atoms)
    apo_heavy = heavy_only(apo_atoms)
    close_sq = close_cutoff * close_cutoff
    severe_sq = severe_cutoff * severe_cutoff

    for ligand_atom in ligand_heavy:
        for apo_atom in apo_heavy:
            dist_sq = (ligand_atom.x - apo_atom.x) ** 2 + (ligand_atom.y - apo_atom.y) ** 2 + (
                ligand_atom.z - apo_atom.z
            ) ** 2
            if dist_sq < close_sq:
                contact_type = "severe_clash" if dist_sq < severe_sq else "close_contact"
                contacts.append(
                    Contact(
                        contact_type=contact_type,
                        distance=math.sqrt(dist_sq),
                        ligand_atom=ligand_atom,
                        apo_atom=apo_atom,
                        apo_category=classify_residue(apo_atom.resname),
                    )
                )
    contacts.sort(key=lambda contact: contact.distance)
    return contacts


def category_counts(contacts: list[Contact]) -> Counter[str]:
    return Counter(contact.apo_category for contact in contacts)


def contact_rows(contacts: list[Contact]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for contact in contacts:
        rows.append(
            {
                "contact_type": contact.contact_type,
                "distance_A": f"{contact.distance:.3f}",
                "apo_category": contact.apo_category,
                "ligand_serial": contact.ligand_atom.serial,
                "ligand_atom": contact.ligand_atom.name,
                "ligand_resname": contact.ligand_atom.resname,
                "ligand_chain": contact.ligand_atom.chain,
                "ligand_resid": contact.ligand_atom.resid,
                "apo_serial": contact.apo_atom.serial,
                "apo_atom": contact.apo_atom.name,
                "apo_resname": contact.apo_atom.resname,
                "apo_chain": contact.apo_atom.chain,
                "apo_resid": contact.apo_atom.resid,
                "apo_element": contact.apo_atom.element,
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames = [
        "contact_type",
        "distance_A",
        "apo_category",
        "ligand_serial",
        "ligand_atom",
        "ligand_resname",
        "ligand_chain",
        "ligand_resid",
        "apo_serial",
        "apo_atom",
        "apo_resname",
        "apo_chain",
        "apo_resid",
        "apo_element",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def count_table(counter: Counter[str]) -> list[str]:
    categories = ["protein", "lipid", "water", "ion", "cholesterol", "unknown"]
    lines = ["| Category | Count |", "|---|---:|"]
    for category in categories:
        lines.append(f"| {category} | {counter.get(category, 0)} |")
    return lines


def write_markdown(
    path: Path,
    apo_atoms: list[PdbAtom],
    ligand_atoms: list[PdbAtom],
    contacts: list[Contact],
) -> None:
    severe = [contact for contact in contacts if contact.contact_type == "severe_clash"]
    close = [contact for contact in contacts if contact.distance < CLOSE_CUTOFF_A]
    lines = [
        "# Ligand-Apo Clash Precheck",
        "",
        "## Scope",
        "",
        "Pilot: ligand 14.2 / L002 / 20chol. This check reports heavy-atom contacts only and does not delete or alter any apo atoms.",
        "",
        "## Inputs",
        "",
        f"- Apo PDB: `{APO_PDB}`",
        f"- Ligand all-atom PDB: `{LIGAND_PDB}`",
        "",
        "## Summary",
        "",
        f"- Apo atoms parsed: {len(apo_atoms)}",
        f"- Apo heavy atoms checked: {len(heavy_only(apo_atoms))}",
        f"- Ligand atoms parsed: {len(ligand_atoms)}",
        f"- Ligand heavy atoms checked: {len(heavy_only(ligand_atoms))}",
        f"- Severe clashes, ligand heavy to apo heavy < {SEVERE_CUTOFF_A:.1f} A: {len(severe)}",
        f"- Close contacts, ligand heavy to apo heavy < {CLOSE_CUTOFF_A:.1f} A: {len(close)}",
        "",
        "Severe clashes are a subset of the close-contact count.",
        "",
        "## Severe Clash Categories",
        "",
        *count_table(category_counts(severe)),
        "",
        "## Close Contact Categories",
        "",
        *count_table(category_counts(close)),
    ]
    if contacts:
        lines.extend(
            [
                "",
                "## Closest Contacts",
                "",
                "| Type | Distance A | Apo category | Ligand atom | Apo atom | Apo residue |",
                "|---|---:|---|---|---|---|",
            ]
        )
        for contact in contacts[:20]:
            lines.append(
                f"| {contact.contact_type} | {contact.distance:.3f} | {contact.apo_category} | "
                f"{contact.ligand_atom.name} {contact.ligand_atom.serial} | "
                f"{contact.apo_atom.name} {contact.apo_atom.serial} | "
                f"{contact.apo_atom.resname} {contact.apo_atom.chain}:{contact.apo_atom.resid} |"
            )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_precheck() -> tuple[int, int, list[Contact]]:
    apo_atoms = read_pdb_atoms(APO_PDB)
    ligand_atoms = read_pdb_atoms(LIGAND_PDB)
    contacts = find_contacts(ligand_atoms, apo_atoms)
    rows = contact_rows(contacts)
    write_csv(CSV_OUT, rows)
    write_markdown(MD_OUT, apo_atoms, ligand_atoms, contacts)
    severe_count = sum(1 for contact in contacts if contact.contact_type == "severe_clash")
    close_count = sum(1 for contact in contacts if contact.distance < CLOSE_CUTOFF_A)
    return severe_count, close_count, contacts


def main() -> int:
    severe_count, close_count, _contacts = run_precheck()
    print(f"severe_clashes_lt_{SEVERE_CUTOFF_A:.1f}A: {severe_count}")
    print(f"close_contacts_lt_{CLOSE_CUTOFF_A:.1f}A: {close_count}")
    print(f"clash CSV: {CSV_OUT}")
    print(f"clash Markdown: {MD_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
