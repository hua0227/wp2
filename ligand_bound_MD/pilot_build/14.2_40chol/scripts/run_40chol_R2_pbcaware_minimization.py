from __future__ import annotations

import math

import openmm
from openmm import LocalEnergyMinimizer, Platform, VerletIntegrator
from openmm.app import HBonds
from openmm.unit import kilojoules_per_mole, nanometer, picoseconds

from common_40chol import (
    CLEANED_PDB,
    CLEANED_PSF,
    MAX_MINIMIZATION_ITERATIONS,
    MINIMIZED_PDB,
    REPORTS_DIR,
    build_system,
    choose_opencl_platform,
    coords_have_nan,
    ensure_dirs,
    force_rows,
    pbc_restraint_force,
    positions_to_angstrom,
    r2_restraint_map,
    read_pdb_atoms,
    retention_and_clashes,
    write_csv,
    write_pdb,
)


SUMMARY_CSV = REPORTS_DIR / "R2_minimization_40chol_summary.csv"
SUMMARY_MD = REPORTS_DIR / "R2_minimization_40chol_summary.md"


def describe_force_atom(row: dict[str, object]) -> str:
    return f"{row['atom_name']} {row['resname']} {row['segid']}:{row['resid']} index={row['index']}"


def main() -> int:
    ensure_dirs()
    print("40chol R2 PBC-aware restrained minimization")
    print("OpenCL only. No MD will be run by this script.")
    print(f"Input cleaned PSF: {CLEANED_PSF}")
    print(f"Input cleaned PDB: {CLEANED_PDB}")
    if not CLEANED_PSF.exists() or not CLEANED_PDB.exists():
        print("FAILED: cleaned PSF/PDB missing")
        return 2

    atoms = read_pdb_atoms(CLEANED_PDB)
    psf, pdb, system = build_system(CLEANED_PSF, CLEANED_PDB, constraints=HBonds, rigid_water=True)
    restraint_map = r2_restraint_map(atoms, protein_k=500.0, ligand_k=25.0, lipid_k=5.0)
    restraint_force = pbc_restraint_force(restraint_map, pdb.positions)
    restraint_force.setForceGroup(system.getNumForces())
    system.addForce(restraint_force)
    platform = choose_opencl_platform(Platform)
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = openmm.Context(system, integrator, platform)
    context.setPositions(pdb.positions)
    ref_coords = positions_to_angstrom(pdb.positions)

    initial_state = context.getState(getEnergy=True, getForces=True)
    initial_energy = initial_state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    print(f"Selected OpenMM platform: {platform.getName()}")
    print(f"R2 restraint atoms: {len(restraint_map)}")
    print(f"Initial restrained energy (kJ/mol): {initial_energy:.12g}")

    LocalEnergyMinimizer.minimize(
        context,
        tolerance=1.0 * kilojoules_per_mole / nanometer,
        maxIterations=MAX_MINIMIZATION_ITERATIONS,
    )
    state = context.getState(getEnergy=True, getForces=True, getPositions=True)
    final_energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    coords = positions_to_angstrom(state.getPositions(asNumpy=False))
    rows, force_nan = force_rows(atoms, state.getForces(asNumpy=False))
    top_force = rows[0]
    retention, clashes = retention_and_clashes(atoms, ref_coords, coords)
    coord_nan = coords_have_nan(coords)
    energy_finite = math.isfinite(final_energy)
    key_dist = retention["key_distances_A"]
    passed_gate = (
        energy_finite
        and not force_nan
        and not coord_nan
        and clashes["severe_lt_1A"] == 0
        and retention["near_chains"] == 2
        and all(value <= 6.0 for value in key_dist.values())
    )

    write_pdb(
        CLEANED_PDB,
        coords,
        MINIMIZED_PDB,
        [
            "REMARK 40chol R2 PBC-aware restrained minimization",
            "REMARK no MD was run",
            f"REMARK final energy {final_energy:.12g} kJ/mol",
        ],
    )

    summary = {
        "initial_energy_kJ_mol": f"{initial_energy:.12g}",
        "final_energy_kJ_mol": f"{final_energy:.12g}",
        "energy_finite": int(energy_finite),
        "force_nan": int(force_nan),
        "coordinate_nan": int(coord_nan),
        "max_force_kJ_mol_nm": f"{top_force['force_magnitude_kJ_mol_nm']:.12g}",
        "top_force_atom": describe_force_atom(top_force),
        "ligand_min_distance_TYR13_A": f"{key_dist['TYR13']:.6g}",
        "ligand_min_distance_VAL17_A": f"{key_dist['VAL17']:.6g}",
        "ligand_min_distance_SER20_A": f"{key_dist['SER20']:.6g}",
        "ligand_key_min_distance_A": f"{retention['key_min_distance_A']:.6g}",
        "ligand_near_chains_count": retention["near_chains"],
        "ligand_rmsd_vs_cleaned_start_A": f"{retention['ligand_rmsd_A']:.6g}",
        "severe_lt_1A": clashes["severe_lt_1A"],
        "close_lt_1p5A": clashes["close_lt_1p5A"],
        "restraint_atom_count": len(restraint_map),
        "probe_gate_passed": int(passed_gate),
        "output_pdb": str(MINIMIZED_PDB),
    }
    write_csv(SUMMARY_CSV, [summary])
    lines = [
        "# 40chol R2 Minimization Summary",
        "",
        "## Scope",
        "",
        "14.2 / L002 / 40chol non-production pilot. OpenCL was used. No MD was run in this minimization step.",
        "",
        "## Results",
        "",
        f"- Cleaned PSF: `{CLEANED_PSF}`",
        f"- Cleaned PDB: `{CLEANED_PDB}`",
        f"- Minimized PDB: `{MINIMIZED_PDB}`",
        f"- Initial restrained energy: {summary['initial_energy_kJ_mol']} kJ/mol",
        f"- Final restrained energy: {summary['final_energy_kJ_mol']} kJ/mol",
        f"- Energy finite: {bool(energy_finite)}",
        f"- Force NaN: {bool(force_nan)}",
        f"- Coordinate NaN: {bool(coord_nan)}",
        f"- Max force: {summary['max_force_kJ_mol_nm']} kJ/mol/nm",
        f"- Top force atom: {summary['top_force_atom']}",
        f"- Ligand min distance to TYR13: {summary['ligand_min_distance_TYR13_A']} A",
        f"- Ligand min distance to VAL17: {summary['ligand_min_distance_VAL17_A']} A",
        f"- Ligand min distance to SER20: {summary['ligand_min_distance_SER20_A']} A",
        f"- Ligand near chains count: {summary['ligand_near_chains_count']}",
        f"- Severe clashes < 1.0 A after minimization: {summary['severe_lt_1A']}",
        f"- Close contacts < 1.5 A after minimization: {summary['close_lt_1p5A']}",
        f"- Very-short probe gate passed: {bool(passed_gate)}",
        "",
        "R2 restraints used periodicdistance anchors: protein backbone medium, ligand heavy atoms weak, lipid/cholesterol heavy atoms extremely weak, water/ions unrestrained.",
    ]
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Final restrained energy (kJ/mol): {final_energy:.12g}")
    print(f"Max force (kJ/mol/nm): {top_force['force_magnitude_kJ_mol_nm']:.12g}")
    print(f"Top force atom: {summary['top_force_atom']}")
    print(f"Ligand min distance TYR13/VAL17/SER20 (A): {key_dist['TYR13']:.6g}/{key_dist['VAL17']:.6g}/{key_dist['SER20']:.6g}")
    print(f"Ligand near chains count: {retention['near_chains']}")
    print(f"Severe/close after minimization: {clashes['severe_lt_1A']}/{clashes['close_lt_1p5A']}")
    print(f"Very-short probe gate passed: {passed_gate}")
    print(f"Summary CSV: {SUMMARY_CSV}")
    print(f"Summary Markdown: {SUMMARY_MD}")
    del context
    del integrator
    return 0 if energy_finite and not force_nan and not coord_nan else 1


if __name__ == "__main__":
    raise SystemExit(main())
