from __future__ import annotations

import math
from pathlib import Path

import parmed as pmd

from branchB_common import (
    BRANCH_B_MINIMIZED_PDB,
    CLEANED_PDB,
    CLEANED_PSF,
    GRO_BOX,
    OUTPUTS_DIR,
    REPORTS_DIR,
    build_charmm_system,
    choose_opencl_platform,
    direct_rmsd_a,
    ensure_dirs,
    is_heavy,
    is_ligand,
    is_water_or_ion,
    positions_to_angstrom,
    read_gro_box_lengths_a,
    read_pdb_atoms,
    write_csv,
    write_pdb_with_coords,
)


CSV_OUT = REPORTS_DIR / "coordinate_handoff_check.csv"
MD_OUT = REPORTS_DIR / "coordinate_handoff_check.md"
ROUNDTRIP_PDB = OUTPUTS_DIR / "branchB_minimized_roundtrip.pdb"


def psf_atoms() -> list[dict[str, str]]:
    structure = pmd.load_file(str(CLEANED_PSF))
    rows: list[dict[str, str]] = []
    for atom in structure.atoms:
        rows.append(
            {
                "name": atom.name,
                "resname": atom.residue.name,
                "resid": str(atom.residue.number),
                "segid": getattr(atom.residue, "segid", ""),
                "element": atom.element_name if getattr(atom, "element_name", "") else "",
            }
        )
    return rows


def compare_psf_pdb_order(psf_rows: list[dict[str, str]], pdb_path: Path) -> tuple[bool, list[str], list[str]]:
    atoms = read_pdb_atoms(pdb_path)
    mismatches: list[str] = []
    warnings: list[str] = []
    if len(psf_rows) != len(atoms):
        mismatches.append(f"atom count mismatch: psf={len(psf_rows)} pdb={len(atoms)}")
        return False, mismatches, warnings
    for index, (psf_atom, pdb_atom) in enumerate(zip(psf_rows, atoms)):
        core_expected = (psf_atom["name"], psf_atom["resname"], psf_atom["segid"])
        core_observed = (pdb_atom.name, pdb_atom.resname, pdb_atom.segid)
        if core_expected != core_observed:
            mismatches.append(f"index {index}: psf_core={core_expected} pdb_core={core_observed}")
            if len(mismatches) >= 20:
                break
        elif psf_atom["resid"] != pdb_atom.resid:
            if is_water_or_ion(pdb_atom):
                if len(warnings) < 20:
                    warnings.append(
                        f"index {index}: water/ion resid differs psf={psf_atom['resid']} pdb={pdb_atom.resid}"
                    )
            else:
                mismatches.append(
                    f"index {index}: non-water resid differs psf={(psf_atom['name'], psf_atom['resname'], psf_atom['resid'], psf_atom['segid'])} "
                    f"pdb={(pdb_atom.name, pdb_atom.resname, pdb_atom.resid, pdb_atom.segid)}"
                )
                if len(mismatches) >= 20:
                    break
    return not mismatches, mismatches, warnings


def compare_ligand_names(psf_rows: list[dict[str, str]], pdb_path: Path) -> tuple[bool, str]:
    psf_names = [
        row["name"]
        for row in psf_rows
        if row["segid"] == "LIG" or row["resname"] == "L002"
    ]
    pdb_names = [atom.name for atom in read_pdb_atoms(pdb_path) if is_ligand(atom)]
    return psf_names == pdb_names, ";".join([",".join(psf_names), ",".join(pdb_names)])


def coords_have_nan(coords: list[tuple[float, float, float]]) -> bool:
    return any(not math.isfinite(value) for xyz in coords for value in xyz)


def box_check() -> dict[str, object]:
    import openmm
    from openmm import Platform, VerletIntegrator
    from openmm.app import HBonds
    from openmm.unit import nanometer, picoseconds

    psf, pdb, _params, system = build_charmm_system(CLEANED_PSF, BRANCH_B_MINIMIZED_PDB, constraints=HBonds, rigid_water=True)
    platform = choose_opencl_platform(Platform)
    integrator = VerletIntegrator(0.001 * picoseconds)
    context = openmm.Context(system, integrator, platform)
    context.setPositions(pdb.positions)
    state = context.getState(getPositions=True)
    coords_a = positions_to_angstrom(state.getPositions(asNumpy=False))
    write_pdb_with_coords(
        BRANCH_B_MINIMIZED_PDB,
        coords_a,
        ROUNDTRIP_PDB,
        [
            "REMARK Branch B minimized coordinate handoff roundtrip",
            "REMARK no dynamics was run",
        ],
    )
    vectors = state.getPeriodicBoxVectors()
    lengths_nm = (
        vectors[0].value_in_unit(nanometer).x,
        vectors[1].value_in_unit(nanometer).y,
        vectors[2].value_in_unit(nanometer).z,
    )
    expected_a = read_gro_box_lengths_a(GRO_BOX)
    observed_a = tuple(value * 10.0 for value in lengths_nm)
    max_box_delta_a = max(abs(a - b) for a, b in zip(expected_a, observed_a))
    del context
    del integrator
    return {
        "platform": platform.getName(),
        "roundtrip_coords_a": coords_a,
        "expected_box_A": expected_a,
        "observed_box_A": observed_a,
        "max_box_delta_A": max_box_delta_a,
        "periodic_box_ok": max_box_delta_a < 1.0e-4,
    }


def main() -> int:
    ensure_dirs()
    psf_rows = psf_atoms()
    cleaned_atoms = read_pdb_atoms(CLEANED_PDB)
    minimized_atoms = read_pdb_atoms(BRANCH_B_MINIMIZED_PDB)
    cleaned_order_ok, cleaned_mismatches, cleaned_warnings = compare_psf_pdb_order(psf_rows, CLEANED_PDB)
    minimized_order_ok, minimized_mismatches, minimized_warnings = compare_psf_pdb_order(psf_rows, BRANCH_B_MINIMIZED_PDB)
    ligand_cleaned_ok, ligand_cleaned_names = compare_ligand_names(psf_rows, CLEANED_PDB)
    ligand_minimized_ok, ligand_minimized_names = compare_ligand_names(psf_rows, BRANCH_B_MINIMIZED_PDB)
    context_result = box_check()
    roundtrip_atoms = read_pdb_atoms(ROUNDTRIP_PDB)
    input_coords = [(atom.x, atom.y, atom.z) for atom in minimized_atoms]
    roundtrip_coords = [(atom.x, atom.y, atom.z) for atom in roundtrip_atoms]
    heavy_indices = [index for index, atom in enumerate(minimized_atoms) if is_heavy(atom)]
    ligand_heavy_indices = [index for index, atom in enumerate(minimized_atoms) if is_ligand(atom) and is_heavy(atom)]
    heavy_rmsd = direct_rmsd_a(heavy_indices, input_coords, roundtrip_coords)
    ligand_heavy_rmsd = direct_rmsd_a(ligand_heavy_indices, input_coords, roundtrip_coords)
    input_nan = coords_have_nan(input_coords)
    roundtrip_nan = coords_have_nan(roundtrip_coords)
    handoff_pass = (
        len(psf_rows) == len(cleaned_atoms) == len(minimized_atoms) == len(roundtrip_atoms)
        and cleaned_order_ok
        and minimized_order_ok
        and ligand_cleaned_ok
        and ligand_minimized_ok
        and not input_nan
        and not roundtrip_nan
        and bool(context_result["periodic_box_ok"])
        and heavy_rmsd < 1.0e-3
        and ligand_heavy_rmsd < 1.0e-3
    )
    rows = [
        {"check": "psf_atom_count", "value": len(psf_rows), "pass": 1},
        {"check": "cleaned_pdb_atom_count", "value": len(cleaned_atoms), "pass": int(len(cleaned_atoms) == len(psf_rows))},
        {"check": "minimized_pdb_atom_count", "value": len(minimized_atoms), "pass": int(len(minimized_atoms) == len(psf_rows))},
        {"check": "roundtrip_pdb_atom_count", "value": len(roundtrip_atoms), "pass": int(len(roundtrip_atoms) == len(psf_rows))},
        {"check": "cleaned_psf_pdb_order_consistent", "value": cleaned_mismatches[0] if cleaned_mismatches else "yes", "pass": int(cleaned_order_ok)},
        {"check": "minimized_psf_pdb_order_consistent", "value": minimized_mismatches[0] if minimized_mismatches else "yes", "pass": int(minimized_order_ok)},
        {"check": "cleaned_water_ion_resid_warnings", "value": cleaned_warnings[0] if cleaned_warnings else "none", "pass": 1},
        {"check": "minimized_water_ion_resid_warnings", "value": minimized_warnings[0] if minimized_warnings else "none", "pass": 1},
        {"check": "cleaned_ligand_atom_names_consistent", "value": ligand_cleaned_names, "pass": int(ligand_cleaned_ok)},
        {"check": "minimized_ligand_atom_names_consistent", "value": ligand_minimized_names, "pass": int(ligand_minimized_ok)},
        {"check": "input_minimized_has_nan_coords", "value": input_nan, "pass": int(not input_nan)},
        {"check": "roundtrip_has_nan_coords", "value": roundtrip_nan, "pass": int(not roundtrip_nan)},
        {"check": "roundtrip_heavy_atom_rmsd_A", "value": f"{heavy_rmsd:.8g}", "pass": int(heavy_rmsd < 1.0e-3)},
        {"check": "roundtrip_ligand_heavy_atom_rmsd_A", "value": f"{ligand_heavy_rmsd:.8g}", "pass": int(ligand_heavy_rmsd < 1.0e-3)},
        {"check": "periodic_box_expected_A", "value": context_result["expected_box_A"], "pass": 1},
        {"check": "periodic_box_observed_A", "value": context_result["observed_box_A"], "pass": int(context_result["periodic_box_ok"])},
        {"check": "periodic_box_max_delta_A", "value": f"{float(context_result['max_box_delta_A']):.8g}", "pass": int(context_result["periodic_box_ok"])},
        {"check": "openmm_platform", "value": context_result["platform"], "pass": int(context_result["platform"] == "OpenCL")},
        {"check": "coordinate_handoff_pass", "value": handoff_pass, "pass": int(handoff_pass)},
    ]
    write_csv(CSV_OUT, rows, ["check", "value", "pass"])
    md_lines = [
        "# Branch B Coordinate Handoff Check",
        "",
        "No dynamics was run. The minimized coordinates were loaded into an OpenMM Context and immediately written back out.",
        "",
        f"- Cleaned PSF atoms: {len(psf_rows)}",
        f"- Cleaned PDB atoms: {len(cleaned_atoms)}",
        f"- Minimized PDB atoms: {len(minimized_atoms)}",
        f"- Roundtrip PDB: `{ROUNDTRIP_PDB}`",
        f"- Heavy-atom roundtrip RMSD: {heavy_rmsd:.8g} A",
        f"- Ligand heavy-atom roundtrip RMSD: {ligand_heavy_rmsd:.8g} A",
        f"- NaN coordinates in input/roundtrip: {input_nan} / {roundtrip_nan}",
        f"- Periodic box max delta vs GRO: {float(context_result['max_box_delta_A']):.8g} A",
        f"- Coordinate handoff pass: {'YES' if handoff_pass else 'NO'}",
        f"- Water/ion residue-number warnings: {len(cleaned_warnings)} cleaned, {len(minimized_warnings)} minimized. These do not affect atom order because name/resname/segid order is consistent.",
        "",
        "| Check | Value | Pass |",
        "|---|---|---:|",
    ]
    for row in rows:
        md_lines.append(f"| {row['check']} | {row['value']} | {row['pass']} |")
    MD_OUT.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    print(f"coordinate_handoff_pass: {handoff_pass}")
    print(f"psf_atom_count: {len(psf_rows)}")
    print(f"cleaned_pdb_atom_count: {len(cleaned_atoms)}")
    print(f"minimized_pdb_atom_count: {len(minimized_atoms)}")
    print(f"roundtrip_heavy_atom_rmsd_A: {heavy_rmsd:.8g}")
    print(f"roundtrip_ligand_heavy_atom_rmsd_A: {ligand_heavy_rmsd:.8g}")
    print(f"periodic_box_ok: {context_result['periodic_box_ok']}")
    print(f"CSV: {CSV_OUT}")
    print(f"Markdown: {MD_OUT}")
    return 0 if handoff_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
