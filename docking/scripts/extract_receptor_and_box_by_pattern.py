from pathlib import Path
import csv
import mdtraj as md
import numpy as np

systems = {
    "20chol": {
        "full": Path(r"C:\TRKB_WP2\docking\receptors\20chol\TRKB_20chol_full.pdb"),
        "protein": Path(r"C:\TRKB_WP2\docking\receptors\20chol\TRKB_20chol_protein.pdb"),
        "box": Path(r"C:\TRKB_WP2\docking\configs\TRKB_20chol_box.txt"),
        "report": Path(r"C:\TRKB_WP2\docking\receptors\20chol\TRKB_20chol_target_residue_mapping.csv"),
    },
    "40chol": {
        "full": Path(r"C:\TRKB_WP2\docking\receptors\40chol\TRKB_40chol_full.pdb"),
        "protein": Path(r"C:\TRKB_WP2\docking\receptors\40chol\TRKB_40chol_protein.pdb"),
        "box": Path(r"C:\TRKB_WP2\docking\configs\TRKB_40chol_box.txt"),
        "report": Path(r"C:\TRKB_WP2\docking\receptors\40chol\TRKB_40chol_target_residue_mapping.csv"),
    },
}

for name, paths in systems.items():
    print(f"\n=== Processing {name} ===")

    traj = md.load(str(paths["full"]))

    protein_atoms = traj.topology.select("protein")
    if len(protein_atoms) == 0:
        raise RuntimeError(f"No protein atoms found in {paths['full']}")

    protein = traj.atom_slice(protein_atoms)
    protein.save_pdb(str(paths["protein"]))
    print("Saved protein receptor:", paths["protein"])
    print("Protein atoms:", len(protein_atoms))

    chains = {}
    for res in protein.topology.residues:
        if res.is_protein:
            chains.setdefault(res.chain.index, []).append(res)

    target_residues = []
    mapping_rows = []

    for chain_id, residues in chains.items():
        for i, res in enumerate(residues):
            if i + 7 >= len(residues):
                continue

            r0 = residues[i]
            r4 = residues[i + 4]
            r7 = residues[i + 7]

            # Target pattern: Y433, V437, S440
            if r0.name == "TYR" and r4.name == "VAL" and r7.name == "SER":
                target_residues.extend([r0, r4, r7])

                mapping_rows.append([name, chain_id, "Y433_like", i, r0.resSeq, r0.name, r0.index])
                mapping_rows.append([name, chain_id, "V437_like", i + 4, r4.resSeq, r4.name, r4.index])
                mapping_rows.append([name, chain_id, "S440_like", i + 7, r7.resSeq, r7.name, r7.index])

                print("FOUND target-like residues:")
                print(f"  Chain {chain_id}: Y433-like resSeq={r0.resSeq}, V437-like resSeq={r4.resSeq}, S440-like resSeq={r7.resSeq}")

    if len(target_residues) == 0:
        raise RuntimeError(
            f"Could not find Y433/V437/S440-like pattern in {name}. "
            f"Inspect protein sequence manually."
        )

    with open(paths["report"], "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "system",
            "chain_index",
            "target_label",
            "index_in_chain",
            "pdb_resSeq",
            "residue_name",
            "mdtraj_residue_index",
        ])
        writer.writerows(mapping_rows)

    print("Saved mapping report:", paths["report"])

    target_atom_indices = []
    for res in target_residues:
        for atom in res.atoms:
            target_atom_indices.append(atom.index)

    xyz_nm = protein.xyz[0, target_atom_indices, :]
    center_nm = xyz_nm.mean(axis=0)
    center_angstrom = center_nm * 10.0

    # Initial search box; later可根据可视化微调
    size_x = 22.0
    size_y = 22.0
    size_z = 24.0

    with open(paths["box"], "w") as f:
        f.write(f"center_x = {center_angstrom[0]:.3f}\n")
        f.write(f"center_y = {center_angstrom[1]:.3f}\n")
        f.write(f"center_z = {center_angstrom[2]:.3f}\n")
        f.write(f"size_x = {size_x:.1f}\n")
        f.write(f"size_y = {size_y:.1f}\n")
        f.write(f"size_z = {size_z:.1f}\n")

    print("Saved docking box:", paths["box"])
    print(f"Docking center Angstrom: {center_angstrom[0]:.3f}, {center_angstrom[1]:.3f}, {center_angstrom[2]:.3f}")