from __future__ import annotations

from collections import Counter

from common_40chol import (
    BOUND_PDB,
    CLEANED_PDB,
    CLEANED_PSF,
    REPORTS_DIR,
    SCRIPTS_DIR,
    choose_cleanup_targets,
    contact_rows,
    copy_bound_to_cleaned,
    ensure_dirs,
    find_contacts,
    make_selection,
    read_pdb_atoms,
    write_csv,
)


CSV_OUT = REPORTS_DIR / "clash_diagnosis_40chol.csv"
MD_OUT = REPORTS_DIR / "clash_diagnosis_40chol.md"
CLEANUP_TCL = SCRIPTS_DIR / "make_cleaned_bound_system_40chol.tcl"


def category_table(contacts) -> list[str]:
    counts = Counter(contact.other_category for contact in contacts)
    lines = ["| Category | Count |", "|---|---:|"]
    for category in ["protein", "lipid", "cholesterol", "water", "ion", "unknown"]:
        lines.append(f"| {category} | {counts.get(category, 0)} |")
    return lines


def write_cleanup_tcl(targets) -> None:
    lines = [
        "# Pilot-local cleaned system builder for 14.2 / L002 / 40chol.",
        "# Deletes only dominant severe-clashing whole water/ion/lipid/cholesterol molecules selected by classify_40chol_ligand_clashes.py.",
        "# Does not modify original apo files or ligand .str.",
        f"set input_psf \"{str(BOUND_PDB.with_suffix('.psf')).replace(chr(92), '/')}\"",
        f"set input_pdb \"{str(BOUND_PDB).replace(chr(92), '/')}\"",
        f"set output_psf \"{str(CLEANED_PSF).replace(chr(92), '/')}\"",
        f"set output_pdb \"{str(CLEANED_PDB).replace(chr(92), '/')}\"",
        "",
    ]
    if not targets:
        lines.extend(
            [
                'puts "No targeted cleanup deletion selected. Copying bound PSF/PDB to cleaned outputs."',
                "file copy -force $input_psf $output_psf",
                "file copy -force $input_pdb $output_pdb",
                'puts "Wrote cleaned PSF: $output_psf"',
                'puts "Wrote cleaned PDB: $output_pdb"',
                "quit",
                "",
            ]
        )
    else:
        selection = make_selection(targets)
        lines.extend(
            [
                "mol new $input_psf type psf waitfor all",
                "mol addfile $input_pdb type pdb waitfor all",
                f"set del [atomselect top {{{selection}}}]",
                'puts "Deleting selected whole molecules/residues:"',
                "puts [$del get {segid resid resname name index}]",
                f"set keep [atomselect top {{not ({selection})}}]",
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
    CLEANUP_TCL.write_text("\n".join(lines), encoding="utf-8")


def write_markdown(path, atoms, contacts, targets, cleanup_reason, molecule_counts) -> None:
    severe = [contact for contact in contacts if contact.contact_type == "severe_clash"]
    close = [contact for contact in contacts if contact.distance < 1.5]
    ligand_heavy = sum(1 for atom in atoms if (atom.segid == "LIG" or atom.resname == "L002") and atom.element.upper() != "H")
    other_heavy = sum(1 for atom in atoms if not (atom.segid == "LIG" or atom.resname == "L002") and atom.element.upper() != "H")
    lines = [
        "# 40chol Ligand Clash Diagnosis",
        "",
        "## Scope",
        "",
        "Pilot: 14.2 / L002 / 40chol ligand-bound system. Ligand heavy atoms were checked against non-ligand heavy atoms.",
        "",
        "## Summary",
        "",
        f"- Bound PDB: `{BOUND_PDB}`",
        f"- Ligand heavy atoms checked: {ligand_heavy}",
        f"- Non-ligand heavy atoms checked: {other_heavy}",
        f"- Severe clashes < 1.0 A: {len(severe)}",
        f"- Close contacts < 1.5 A: {len(close)}",
        "- Severe clashes are included in the close-contact count.",
        "",
        "## Severe Clash Categories",
        "",
        *category_table(severe),
        "",
        "## Close Contact Categories",
        "",
        *category_table(close),
        "",
        "## Targeted Cleanup Decision",
        "",
        f"- Decision: {cleanup_reason}",
    ]
    if targets:
        lines.append("- Selected whole molecules/residues:")
        for category, segid, resid, resname in targets:
            lines.append(f"  - {category}: {resname} {segid}:{resid}")
    else:
        lines.append("- No molecule was deleted; cleaned outputs are pilot-local copies of the bound PSF/PDB.")
    if molecule_counts:
        lines.extend(["", "## Severe Clash Molecule Counts", "", "| Category | Molecule | Severe count |", "|---|---|---:|"])
        for (category, segid, resid, resname), count in molecule_counts.most_common(15):
            lines.append(f"| {category} | {resname} {segid}:{resid} | {count} |")
    lines.extend(["", "## Closest Contacts", "", "| Type | Distance A | Category | Ligand atom | Other atom | Other residue |", "|---|---:|---|---|---|---|"])
    for contact in contacts[:50]:
        lig = contact.ligand_atom
        other = contact.other_atom
        lines.append(
            f"| {contact.contact_type} | {contact.distance:.3f} | {contact.other_category} | "
            f"{lig.name} {lig.segid}:{lig.resid} | {other.name} {other.serial} | {other.resname} {other.segid}:{other.resid} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    ensure_dirs()
    atoms = read_pdb_atoms(BOUND_PDB)
    contacts = find_contacts(atoms)
    targets, cleanup_reason, molecule_counts = choose_cleanup_targets(contacts)
    write_csv(CSV_OUT, contact_rows(contacts))
    write_markdown(MD_OUT, atoms, contacts, targets, cleanup_reason, molecule_counts)
    write_cleanup_tcl(targets)
    if not targets:
        copy_bound_to_cleaned()
    severe_count = sum(1 for contact in contacts if contact.contact_type == "severe_clash")
    close_count = sum(1 for contact in contacts if contact.distance < 1.5)
    print(f"severe_clashes_lt_1.0A: {severe_count}")
    print(f"close_contacts_lt_1.5A: {close_count}")
    print(f"targeted_cleanup_decision: {cleanup_reason}")
    print(f"selected_deletions: {len(targets)}")
    for category, segid, resid, resname in targets:
        print(f"- {category} {resname} {segid}:{resid}")
    print(f"classification CSV: {CSV_OUT}")
    print(f"classification Markdown: {MD_OUT}")
    print(f"cleaned-system Tcl: {CLEANUP_TCL}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
