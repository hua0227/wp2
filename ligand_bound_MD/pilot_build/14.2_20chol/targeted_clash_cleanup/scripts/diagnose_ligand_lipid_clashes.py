from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any


PILOT_ROOT = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol")
CLEANUP_ROOT = PILOT_ROOT / "targeted_clash_cleanup"
REPORTS_DIR = CLEANUP_ROOT / "reports"

BOUND_PSF = PILOT_ROOT / "outputs" / "TRKB_20chol_L002_14.2_bound.psf"
BOUND_PDB = PILOT_ROOT / "outputs" / "TRKB_20chol_L002_14.2_bound.pdb"
PREVIOUS_MINIMIZED_PDB = PILOT_ROOT / "clash_resolution" / "outputs" / "TRKB_20chol_L002_14.2_minimized.pdb"

CSV_OUT = REPORTS_DIR / "detailed_ligand_lipid_clashes.csv"
MD_OUT = REPORTS_DIR / "detailed_ligand_lipid_clashes.md"

SEVERE_CUTOFF_A = 1.0
CLOSE_CUTOFF_A = 1.5
TARGET_LIGAND_ATOMS = {"H9", "H11"}
TARGET_LIPID_SEGID = "MEMB"
TARGET_LIPID_RESID = "94"
TARGET_LIPID_RESNAME = "POPC"

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
    index: int
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
    source_label: str
    distance_a: float
    cutoff_class: str
    atom_pair_class: str
    ligand_atom: PdbAtom
    other_atom: PdbAtom
    other_category: str
    other_residue_key: tuple[str, str, str]
    is_popc94: bool
    involves_h9_h11: bool


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


def parse_pdb_atom_line(line: str, index: int) -> PdbAtom | None:
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
    return PdbAtom(index, serial, name, resname, segid, resid, element, x, y, z)


def read_pdb_atoms(path: Path) -> list[PdbAtom]:
    atoms: list[PdbAtom] = []
    atom_index = 0
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        atom = parse_pdb_atom_line(line, atom_index)
        atom_index += 1
        if atom is not None:
            atoms.append(atom)
    return atoms


def is_ligand(atom: PdbAtom) -> bool:
    return atom.segid == "LIG" or atom.resname == "L002"


def is_hydrogen(atom: PdbAtom) -> bool:
    return atom.element.upper() == "H"


def atom_pair_class(ligand_atom: PdbAtom, other_atom: PdbAtom) -> str:
    ligand_part = "H" if is_hydrogen(ligand_atom) else "heavy"
    other_part = "H" if is_hydrogen(other_atom) else "heavy"
    return f"{ligand_part}-{other_part}"


def classify_atom(atom: PdbAtom) -> str:
    resname = atom.resname.upper()
    segid = atom.segid.upper()
    if resname in PROTEIN_RESNAMES or segid.startswith("PRO"):
        return "protein"
    if resname in CHOLESTEROL_RESNAMES:
        return "cholesterol"
    if resname in LIPID_RESNAMES or segid == "MEMB":
        return "lipid"
    if resname in WATER_RESNAMES or segid.startswith("WT"):
        return "water"
    if resname in ION_RESNAMES or segid.startswith("ION"):
        return "ion"
    return "other/unknown"


def residue_key(atom: PdbAtom) -> tuple[str, str, str]:
    return (atom.segid, atom.resid, atom.resname)


def is_popc94(atom: PdbAtom) -> bool:
    return (
        atom.segid == TARGET_LIPID_SEGID
        and atom.resid == TARGET_LIPID_RESID
        and atom.resname == TARGET_LIPID_RESNAME
    )


def distance_a(a: PdbAtom, b: PdbAtom) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def find_contacts(source_label: str, atoms: list[PdbAtom]) -> list[Contact]:
    ligand_atoms = [atom for atom in atoms if is_ligand(atom)]
    other_atoms = [atom for atom in atoms if not is_ligand(atom)]
    contacts: list[Contact] = []
    close_sq = CLOSE_CUTOFF_A * CLOSE_CUTOFF_A
    severe_sq = SEVERE_CUTOFF_A * SEVERE_CUTOFF_A
    for ligand_atom in ligand_atoms:
        for other_atom in other_atoms:
            dx = ligand_atom.x - other_atom.x
            dy = ligand_atom.y - other_atom.y
            dz = ligand_atom.z - other_atom.z
            d2 = dx * dx + dy * dy + dz * dz
            if d2 >= close_sq:
                continue
            cutoff_class = "severe_lt_1.0A" if d2 < severe_sq else "close_1.0_to_1.5A"
            contacts.append(
                Contact(
                    source_label=source_label,
                    distance_a=math.sqrt(d2),
                    cutoff_class=cutoff_class,
                    atom_pair_class=atom_pair_class(ligand_atom, other_atom),
                    ligand_atom=ligand_atom,
                    other_atom=other_atom,
                    other_category=classify_atom(other_atom),
                    other_residue_key=residue_key(other_atom),
                    is_popc94=is_popc94(other_atom),
                    involves_h9_h11=ligand_atom.name in TARGET_LIGAND_ATOMS,
                )
            )
    contacts.sort(key=lambda item: (item.source_label, item.distance_a))
    return contacts


def contact_rows(contacts: list[Contact]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for contact in contacts:
        ligand = contact.ligand_atom
        other = contact.other_atom
        rows.append(
            {
                "source": contact.source_label,
                "cutoff_class": contact.cutoff_class,
                "distance_A": f"{contact.distance_a:.4f}",
                "atom_pair_class": contact.atom_pair_class,
                "other_category": contact.other_category,
                "is_MEMB_94_POPC": int(contact.is_popc94),
                "involves_LIG_H9_or_H11": int(contact.involves_h9_h11),
                "ligand_index": ligand.index,
                "ligand_serial": ligand.serial,
                "ligand_atom_name": ligand.name,
                "ligand_element": ligand.element,
                "ligand_resname": ligand.resname,
                "ligand_segid": ligand.segid,
                "ligand_resid": ligand.resid,
                "other_index": other.index,
                "other_serial": other.serial,
                "other_atom_name": other.name,
                "other_element": other.element,
                "other_resname": other.resname,
                "other_segid": other.segid,
                "other_resid": other.resid,
            }
        )
    return rows


def write_csv(path: Path, contacts: list[Contact]) -> None:
    fieldnames = [
        "source",
        "cutoff_class",
        "distance_A",
        "atom_pair_class",
        "other_category",
        "is_MEMB_94_POPC",
        "involves_LIG_H9_or_H11",
        "ligand_index",
        "ligand_serial",
        "ligand_atom_name",
        "ligand_element",
        "ligand_resname",
        "ligand_segid",
        "ligand_resid",
        "other_index",
        "other_serial",
        "other_atom_name",
        "other_element",
        "other_resname",
        "other_segid",
        "other_resid",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(contact_rows(contacts))


def count_contacts(contacts: list[Contact], source_label: str | None = None, severe_only: bool = False) -> int:
    return sum(
        1
        for contact in contacts
        if (source_label is None or contact.source_label == source_label)
        and (not severe_only or contact.cutoff_class == "severe_lt_1.0A")
    )


def source_labels(contacts: list[Contact]) -> list[str]:
    return sorted({contact.source_label for contact in contacts})


def popc94_atoms(atoms: list[PdbAtom]) -> list[PdbAtom]:
    return [atom for atom in atoms if is_popc94(atom)]


def ligand_atoms_clashing_with_popc94(contacts: list[Contact], source_label: str, severe_only: bool = True) -> set[str]:
    names: set[str] = set()
    for contact in contacts:
        if contact.source_label != source_label or not contact.is_popc94:
            continue
        if severe_only and contact.cutoff_class != "severe_lt_1.0A":
            continue
        names.add(contact.ligand_atom.name)
    return names


def residue_summary(contacts: list[Contact], source_label: str, severe_only: bool = True) -> list[tuple[tuple[str, str, str], int]]:
    counts: Counter[tuple[str, str, str]] = Counter()
    for contact in contacts:
        if contact.source_label != source_label:
            continue
        if severe_only and contact.cutoff_class != "severe_lt_1.0A":
            continue
        counts[contact.other_residue_key] += 1
    return counts.most_common()


def is_popc94_dominant(contacts: list[Contact], source_label: str) -> bool:
    severe = [contact for contact in contacts if contact.source_label == source_label and contact.cutoff_class == "severe_lt_1.0A"]
    if not severe:
        return False
    by_residue = residue_summary(contacts, source_label, severe_only=True)
    if not by_residue:
        return False
    top_key, top_count = by_residue[0]
    popc94_count = sum(1 for contact in severe if contact.is_popc94)
    return top_key == (TARGET_LIPID_SEGID, TARGET_LIPID_RESID, TARGET_LIPID_RESNAME) and popc94_count == top_count


def h9_h11_lines(contacts: list[Contact], source_label: str) -> list[str]:
    rows = [
        contact
        for contact in contacts
        if contact.source_label == source_label and contact.involves_h9_h11
    ]
    lines = [
        "| Cutoff | Distance A | Ligand | Other atom | Other residue | Pair class | POPC94 |",
        "|---|---:|---|---|---|---|---:|",
    ]
    for contact in rows[:40]:
        other = contact.other_atom
        lines.append(
            f"| {contact.cutoff_class} | {contact.distance_a:.4f} | {contact.ligand_atom.name} | "
            f"{other.name} {other.serial} | {other.resname} {other.segid}:{other.resid} | "
            f"{contact.atom_pair_class} | {int(contact.is_popc94)} |"
        )
    if not rows:
        lines.append("| none |  |  |  |  |  | 0 |")
    return lines


def write_markdown(path: Path, atoms_by_source: dict[str, list[PdbAtom]], contacts: list[Contact]) -> None:
    lines = [
        "# Detailed Ligand-Lipid Clash Diagnosis",
        "",
        "## Scope",
        "",
        "- Pilot: ligand_id 14.2 / resname L002 / system 20chol.",
        f"- Bound PSF: `{BOUND_PSF}`",
        f"- Bound PDB: `{BOUND_PDB}`",
        f"- Previous minimized PDB checked: `{PREVIOUS_MINIMIZED_PDB}`",
        "- Contacts include all ligand atoms, including hydrogens, against all non-ligand atoms.",
        f"- Severe cutoff: < {SEVERE_CUTOFF_A:.1f} A.",
        f"- Close-contact cutoff: < {CLOSE_CUTOFF_A:.1f} A, with severe contacts included in the total close count.",
        "",
        "## Summary",
        "",
        "| Source | Severe <1.0 A | Close <1.5 A | POPC94 severe | POPC94 close | Ligand atoms severe with POPC94 | H9/H11 severe with POPC94 | MEMB:94 POPC atoms | POPC94 dominant severe source |",
        "|---|---:|---:|---:|---:|---|---:|---:|---:|",
    ]
    for source in source_labels(contacts):
        severe = [c for c in contacts if c.source_label == source and c.cutoff_class == "severe_lt_1.0A"]
        close = [c for c in contacts if c.source_label == source]
        popc_severe = [c for c in severe if c.is_popc94]
        popc_close = [c for c in close if c.is_popc94]
        ligand_names = sorted(ligand_atoms_clashing_with_popc94(contacts, source, severe_only=True))
        h9_h11_popc_severe = [c for c in popc_severe if c.involves_h9_h11]
        lines.append(
            f"| {source} | {len(severe)} | {len(close)} | {len(popc_severe)} | {len(popc_close)} | "
            f"{', '.join(ligand_names) if ligand_names else 'none'} | {len(h9_h11_popc_severe)} | "
            f"{len(popc94_atoms(atoms_by_source[source]))} | {int(is_popc94_dominant(contacts, source))} |"
        )

    lines.extend(
        [
            "",
            "## Atom Pair Classes",
            "",
            "| Source | Cutoff | heavy-heavy | heavy-H | H-heavy | H-H |",
            "|---|---|---:|---:|---:|---:|",
        ]
    )
    for source in source_labels(contacts):
        for cutoff in ["severe_lt_1.0A", "close_1.0_to_1.5A"]:
            subset = [c for c in contacts if c.source_label == source and c.cutoff_class == cutoff]
            counts = Counter(c.atom_pair_class for c in subset)
            lines.append(
                f"| {source} | {cutoff} | {counts.get('heavy-heavy', 0)} | {counts.get('heavy-H', 0)} | "
                f"{counts.get('H-heavy', 0)} | {counts.get('H-H', 0)} |"
            )

    lines.extend(["", "## Top Severe-Clashing Residues", ""])
    for source in source_labels(contacts):
        lines.extend(
            [
                f"### {source}",
                "",
                "| Rank | Residue | Severe contact count | POPC94 |",
                "|---:|---|---:|---:|",
            ]
        )
        for rank, (key, count) in enumerate(residue_summary(contacts, source, severe_only=True)[:10], start=1):
            segid, resid, resname = key
            lines.append(
                f"| {rank} | {resname} {segid}:{resid} | {count} | "
                f"{int(key == (TARGET_LIPID_SEGID, TARGET_LIPID_RESID, TARGET_LIPID_RESNAME))} |"
            )
        if not residue_summary(contacts, source, severe_only=True):
            lines.append("|  | none | 0 | 0 |")
        lines.append("")

    lines.extend(["## LIG H9/H11 Contacts", ""])
    for source in source_labels(contacts):
        lines.extend([f"### {source}", "", *h9_h11_lines(contacts, source), ""])

    lines.extend(
        [
            "## Closest Contacts",
            "",
            "| Source | Cutoff | Distance A | Pair class | Ligand atom | Other atom | Other residue | Category | POPC94 |",
            "|---|---|---:|---|---|---|---|---|---:|",
        ]
    )
    for contact in sorted(contacts, key=lambda item: item.distance_a)[:80]:
        other = contact.other_atom
        lines.append(
            f"| {contact.source_label} | {contact.cutoff_class} | {contact.distance_a:.4f} | {contact.atom_pair_class} | "
            f"{contact.ligand_atom.name} {contact.ligand_atom.segid}:{contact.ligand_atom.resid} | "
            f"{other.name} {other.serial} | {other.resname} {other.segid}:{other.resid} | "
            f"{contact.other_category} | {int(contact.is_popc94)} |"
        )

    lines.extend(
        [
            "",
            "## Deletion Decision Input",
            "",
            "Branch B is permitted to remove only the whole POPC molecule MEMB:94 in a pilot-local cleaned system. "
            "The deletion decision uses the bound Branch B starting coordinates plus the prior failed-MD force diagnostic that implicated LIG H9/H11 and MEMB:94 POPC. "
            "The previous minimized coordinates are also listed because they show whether simple geometric contacts remained after the earlier minimization.",
        ]
    )
    bound_dominant = is_popc94_dominant(contacts, "bound_initial")
    minimized_dominant = is_popc94_dominant(contacts, "previous_minimized")
    lines.append(f"- MEMB:94 POPC dominant severe source in bound_initial: {'YES' if bound_dominant else 'NO'}")
    lines.append(f"- MEMB:94 POPC dominant severe source in previous_minimized: {'YES' if minimized_dominant else 'NO'}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_sources() -> dict[str, list[PdbAtom]]:
    sources = {"bound_initial": read_pdb_atoms(BOUND_PDB)}
    if PREVIOUS_MINIMIZED_PDB.exists():
        sources["previous_minimized"] = read_pdb_atoms(PREVIOUS_MINIMIZED_PDB)
    return sources


def main() -> int:
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    atoms_by_source = load_sources()
    all_contacts: list[Contact] = []
    for label, atoms in atoms_by_source.items():
        all_contacts.extend(find_contacts(label, atoms))
    write_csv(CSV_OUT, all_contacts)
    write_markdown(MD_OUT, atoms_by_source, all_contacts)

    print("Detailed ligand-lipid clash diagnosis")
    for source in source_labels(all_contacts):
        severe = [c for c in all_contacts if c.source_label == source and c.cutoff_class == "severe_lt_1.0A"]
        close = [c for c in all_contacts if c.source_label == source]
        popc_severe = [c for c in severe if c.is_popc94]
        popc_close = [c for c in close if c.is_popc94]
        h9_h11_popc = [c for c in popc_severe if c.involves_h9_h11]
        print(f"{source}: severe_lt_1.0A={len(severe)} close_lt_1.5A={len(close)}")
        print(f"{source}: MEMB94_POPC severe={len(popc_severe)} close={len(popc_close)}")
        print(f"{source}: LIG_H9_H11 severe_with_MEMB94_POPC={len(h9_h11_popc)}")
        print(f"{source}: MEMB94_POPC_dominant={is_popc94_dominant(all_contacts, source)}")
    print(f"CSV: {CSV_OUT}")
    print(f"Markdown: {MD_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
