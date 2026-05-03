from __future__ import annotations

import math
from collections import defaultdict

import parmed as pmd

from residual_common import (
    CLEANED_PSF,
    MEMB49_KEY,
    ORIGINAL_APO_PDB,
    REPORTS_DIR,
    STAGE4_PDB,
    atom_category,
    coords_have_nan,
    distance_a,
    ensure_dirs,
    is_heavy,
    is_ligand,
    is_memb49,
    read_gro_box_lengths_a,
    read_pdb_atoms,
    write_csv,
)


GEOM_CSV = REPORTS_DIR / "MEMB49_geometry_diagnostics.csv"
CONTACTS_CSV = REPORTS_DIR / "MEMB49_contacts.csv"
SUMMARY_JSON = REPORTS_DIR / "MEMB49_geometry_summary.json"


def coord(atom) -> tuple[float, float, float]:
    return (atom.x, atom.y, atom.z)


def psf_memb49_bonds() -> list[tuple[int, str, int, str]]:
    structure = pmd.load_file(str(CLEANED_PSF))
    memb49_atom_indices = {
        atom.idx for atom in structure.atoms if getattr(atom.residue, "segid", "") == "MEMB" and str(atom.residue.number) == "49" and atom.residue.name == "POPC"
    }
    bonds: list[tuple[int, str, int, str]] = []
    for bond in structure.bonds:
        if bond.atom1.idx in memb49_atom_indices or bond.atom2.idx in memb49_atom_indices:
            bonds.append((bond.atom1.idx, bond.atom1.name, bond.atom2.idx, bond.atom2.name))
    return bonds


def mapped_heavy_rmsd(stage_atoms, apo_atoms, box) -> tuple[float, int, str]:
    stage_by_name = {atom.name: atom for atom in stage_atoms if is_memb49(atom) and is_heavy(atom)}
    apo_by_name = {atom.name: atom for atom in apo_atoms if is_memb49(atom) and is_heavy(atom)}
    common = sorted(set(stage_by_name) & set(apo_by_name))
    if not common:
        return float("nan"), 0, "no_common_atoms"
    # Translate apo to the closest image relative to stage by centroid.
    stage_centroid = tuple(sum(getattr(stage_by_name[name], axis) for name in common) / len(common) for axis in ["x", "y", "z"])
    apo_centroid = tuple(sum(getattr(apo_by_name[name], axis) for name in common) / len(common) for axis in ["x", "y", "z"])
    translation = []
    for s, a, length in zip(stage_centroid, apo_centroid, box):
        delta = s - a
        delta -= round(delta / length) * length
        translation.append(delta)
    total = 0.0
    for name in common:
        s = coord(stage_by_name[name])
        a0 = coord(apo_by_name[name])
        a = (a0[0] + translation[0], a0[1] + translation[1], a0[2] + translation[2])
        total += distance_a(s, a) ** 2
    return math.sqrt(total / len(common)), len(common), f"{translation[0]:.4f},{translation[1]:.4f},{translation[2]:.4f}"


def main() -> int:
    ensure_dirs()
    atoms = read_pdb_atoms(STAGE4_PDB)
    apo_atoms = read_pdb_atoms(ORIGINAL_APO_PDB)
    box = read_gro_box_lengths_a()
    coords = [coord(atom) for atom in atoms]
    memb49_indices = [atom.index for atom in atoms if is_memb49(atom)]
    ligand_indices = [atom.index for atom in atoms if is_ligand(atom)]
    non_ligand_non_memb49_indices = [atom.index for atom in atoms if not is_ligand(atom) and not is_memb49(atom)]
    memb49_atoms = [atoms[index] for index in memb49_indices]
    abnormal_coord = coords_have_nan([coords[index] for index in memb49_indices])

    bond_rows: list[dict[str, object]] = []
    abnormal_bonds: list[dict[str, object]] = []
    for a1, n1, a2, n2 in psf_memb49_bonds():
        if a1 >= len(atoms) or a2 >= len(atoms):
            continue
        direct = distance_a(coords[a1], coords[a2])
        pbc = distance_a(coords[a1], coords[a2], box)
        flags: list[str] = []
        if pbc < 0.7:
            flags.append("short_lt_0.7A")
        if pbc > 2.0:
            flags.append("long_gt_2.0A")
        if direct > 2.0 and pbc <= 2.0:
            flags.append("wrapped_direct_long")
        row = {
            "record_type": "bond",
            "atom1_index": a1,
            "atom1_name": n1,
            "atom2_index": a2,
            "atom2_name": n2,
            "direct_length_A": f"{direct:.6g}",
            "minimum_image_length_A": f"{pbc:.6g}",
            "flags": ";".join(flags) if flags else "ok",
        }
        bond_rows.append(row)
        if any(flag in flags for flag in ["short_lt_0.7A", "long_gt_2.0A"]):
            abnormal_bonds.append(row)

    ligand_min = min(distance_a(coords[i], coords[j], box) for i in memb49_indices for j in ligand_indices)
    contacts: list[dict[str, object]] = []
    for i in memb49_indices:
        atom_i = atoms[i]
        for j in non_ligand_non_memb49_indices:
            atom_j = atoms[j]
            dist = distance_a(coords[i], coords[j], box)
            contacts.append(
                {
                    "distance_A": dist,
                    "memb49_atom_index": i,
                    "memb49_atom_name": atom_i.name,
                    "other_index": j,
                    "other_atom_name": atom_j.name,
                    "other_resname": atom_j.resname,
                    "other_segid": atom_j.segid,
                    "other_resid": atom_j.resid,
                    "other_category": atom_category(atom_j),
                    "other_is_ligand": int(is_ligand(atom_j)),
                }
            )
    contacts.sort(key=lambda row: row["distance_A"])
    contact_rows = [{**row, "rank": rank, "distance_A": f"{row['distance_A']:.6g}"} for rank, row in enumerate(contacts[:20], start=1)]

    heavy_rmsd, mapped_count, apo_translation = mapped_heavy_rmsd(atoms, apo_atoms, box)
    memb49_centroid = tuple(sum(getattr(atom, axis) for atom in memb49_atoms) / len(memb49_atoms) for axis in ["x", "y", "z"])
    box_margin = min(
        memb49_centroid[0] % box[0],
        box[0] - (memb49_centroid[0] % box[0]),
        memb49_centroid[1] % box[1],
        box[1] - (memb49_centroid[1] % box[1]),
        memb49_centroid[2] % box[2],
        box[2] - (memb49_centroid[2] % box[2]),
    )
    near_boundary = box_margin < 5.0

    summary_rows = [
        {"record_type": "summary", "metric": "MEMB49_atom_count", "value": len(memb49_indices), "flags": "ok"},
        {"record_type": "summary", "metric": "MEMB49_bond_count_from_psf", "value": len(bond_rows), "flags": "ok"},
        {"record_type": "summary", "metric": "abnormal_short_or_long_bond_count", "value": len(abnormal_bonds), "flags": "abnormal" if abnormal_bonds else "ok"},
        {"record_type": "summary", "metric": "MEMB49_ligand_min_distance_A", "value": f"{ligand_min:.6g}", "flags": "clash" if ligand_min < 1.5 else "ok"},
        {"record_type": "summary", "metric": "MEMB49_vs_original_apo_heavy_RMSD_A", "value": f"{heavy_rmsd:.6g}", "flags": "ok" if math.isfinite(heavy_rmsd) else "not_available"},
        {"record_type": "summary", "metric": "MEMB49_vs_original_apo_mapped_heavy_atom_count", "value": mapped_count, "flags": "ok" if mapped_count else "not_available"},
        {"record_type": "summary", "metric": "apo_to_stage4_translation_A", "value": apo_translation, "flags": "ok"},
        {"record_type": "summary", "metric": "MEMB49_centroid_A", "value": f"{memb49_centroid[0]:.4f},{memb49_centroid[1]:.4f},{memb49_centroid[2]:.4f}", "flags": "ok"},
        {"record_type": "summary", "metric": "MEMB49_min_box_boundary_margin_A", "value": f"{box_margin:.6g}", "flags": "near_boundary" if near_boundary else "ok"},
        {"record_type": "summary", "metric": "MEMB49_has_nan_coords", "value": abnormal_coord, "flags": "abnormal" if abnormal_coord else "ok"},
    ]
    geometry_rows = summary_rows + bond_rows
    geometry_fieldnames = [
        "record_type",
        "metric",
        "value",
        "atom1_index",
        "atom1_name",
        "atom2_index",
        "atom2_name",
        "direct_length_A",
        "minimum_image_length_A",
        "flags",
    ]
    write_csv(GEOM_CSV, geometry_rows, geometry_fieldnames)
    write_csv(CONTACTS_CSV, contact_rows)
    from residual_common import write_json

    write_json(
        SUMMARY_JSON,
        {
            "memb49_atom_count": len(memb49_indices),
            "bond_count": len(bond_rows),
            "abnormal_bond_count": len(abnormal_bonds),
            "ligand_min_distance_A": ligand_min,
            "apo_heavy_rmsd_A": heavy_rmsd,
            "apo_mapped_heavy_atom_count": mapped_count,
            "box_margin_A": box_margin,
            "near_boundary": near_boundary,
            "has_nan_coords": abnormal_coord,
            "closest_contacts": contact_rows,
            "abnormal_bonds": abnormal_bonds,
        },
    )
    print(f"MEMB49_abnormal_bonds: {len(abnormal_bonds)}")
    print(f"MEMB49_ligand_min_distance_A: {ligand_min:.6g}")
    print(f"MEMB49_apo_heavy_RMSD_A: {heavy_rmsd:.6g}")
    print(f"MEMB49_min_box_boundary_margin_A: {box_margin:.6g}")
    print(f"MEMB49_near_box_boundary: {near_boundary}")
    print(f"Geometry CSV: {GEOM_CSV}")
    print(f"Contacts CSV: {CONTACTS_CSV}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
