from __future__ import annotations

import math
import sys
from pathlib import Path


PILOT_ROOT = Path(r"C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol")
LOGS_DIR = PILOT_ROOT / "logs"
OUTPUTS_DIR = PILOT_ROOT / "outputs"

BOUND_PSF = OUTPUTS_DIR / "TRKB_20chol_L002_14.2_bound.psf"
BOUND_PDB = OUTPUTS_DIR / "TRKB_20chol_L002_14.2_bound.pdb"
GRO_BOX = Path(r"C:\TRKB_WP2\TRKB_20CHOL\gromacs\step5_input.gro")

TOPPAR_DIR = Path(r"C:\TRKB_WP2\TRKB_20CHOL\toppar")
LIGAND_STR = Path(r"C:\TRKB_WP2\ligand_bound_MD\cgenff_parameterization\returned_str\L002_14.2.str")
PARAMETER_FILES = [
    TOPPAR_DIR / "top_all36_prot.rtf",
    TOPPAR_DIR / "top_all36_lipid.rtf",
    TOPPAR_DIR / "top_all36_cgenff.rtf",
    TOPPAR_DIR / "par_all36m_prot.prm",
    TOPPAR_DIR / "par_all36_lipid.prm",
    TOPPAR_DIR / "par_all36_cgenff.prm",
    TOPPAR_DIR / "toppar_water_ions.str",
    TOPPAR_DIR / "toppar_all36_lipid_cholesterol.str",
    LIGAND_STR,
]


def choose_platform(platform_cls, preferred: str = "OpenCL"):
    try:
        return platform_cls.getPlatformByName(preferred)
    except Exception:
        for index in range(platform_cls.getNumPlatforms()):
            platform = platform_cls.getPlatform(index)
            if platform.getName() != "CUDA":
                return platform
    raise RuntimeError("No non-CUDA OpenMM platform is available")


def read_gro_box_vectors(path: Path):
    from openmm import Vec3
    from openmm.unit import nanometer

    values = [float(value) for value in path.read_text(encoding="utf-8", errors="replace").splitlines()[-1].split()]
    if len(values) == 3:
        ax, by, cz = values
        return (
            Vec3(ax, 0.0, 0.0) * nanometer,
            Vec3(0.0, by, 0.0) * nanometer,
            Vec3(0.0, 0.0, cz) * nanometer,
        )
    if len(values) == 9:
        v1x, v2y, v3z, v1y, v1z, v2x, v2z, v3x, v3y = values
        return (
            Vec3(v1x, v1y, v1z) * nanometer,
            Vec3(v2x, v2y, v2z) * nanometer,
            Vec3(v3x, v3y, v3z) * nanometer,
        )
    raise ValueError(f"Unsupported GRO box vector line with {len(values)} values: {values}")


def main() -> int:
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    print("OpenMM ligand-bound read test for pilot 14.2 / L002 / 20chol")
    print("No minimization and no MD will be run.")
    print(f"Input PSF: {BOUND_PSF}")
    print(f"Input PDB: {BOUND_PDB}")
    print(f"GRO box: {GRO_BOX}")

    if not BOUND_PSF.exists() or not BOUND_PDB.exists():
        print("SKIP: ligand-bound PSF/PDB are missing; psfgen build must succeed before read test.")
        return 2

    try:
        import openmm
        from openmm import Platform, VerletIntegrator
        from openmm.app import CharmmParameterSet, CharmmPsfFile, NoCutoff, PDBFile, PME
        from openmm.unit import kilojoules_per_mole, nanometer, picoseconds
    except Exception as exc:
        print(f"FAILED: OpenMM import failed: {exc!r}")
        return 2

    print("Available OpenMM platforms:")
    for index in range(Platform.getNumPlatforms()):
        print(f"- {Platform.getPlatform(index).getName()}")

    missing = [str(path) for path in [BOUND_PSF, BOUND_PDB, GRO_BOX, *PARAMETER_FILES] if not path.exists()]
    if missing:
        print("FAILED: missing required files:")
        for path in missing:
            print(f"- {path}")
        return 2

    try:
        psf = CharmmPsfFile(str(BOUND_PSF))
        pdb = PDBFile(str(BOUND_PDB))
        box_vectors = read_gro_box_vectors(GRO_BOX)
        if hasattr(psf, "setBoxVectors"):
            psf.setBoxVectors(*box_vectors)
        else:
            a = box_vectors[0].value_in_unit(nanometer)[0]
            b = box_vectors[1].value_in_unit(nanometer)[1]
            c = box_vectors[2].value_in_unit(nanometer)[2]
            psf.setBox(a * nanometer, b * nanometer, c * nanometer)

        params = CharmmParameterSet(*[str(path) for path in PARAMETER_FILES])
        print("CharmmParameterSet loaded.")
        try:
            system = psf.createSystem(
                params,
                nonbondedMethod=PME,
                nonbondedCutoff=1.2 * nanometer,
                switchDistance=1.0 * nanometer,
                constraints=None,
            )
        except Exception as pme_exc:
            print(f"PME createSystem failed: {pme_exc!r}")
            print("Trying NoCutoff only to isolate parameter/topology readability.")
            system = psf.createSystem(params, nonbondedMethod=NoCutoff, constraints=None)

        selected_platform = choose_platform(Platform, preferred="OpenCL")
        print(f"Selected OpenMM platform: {selected_platform.getName()}")
        integrator = VerletIntegrator(0.001 * picoseconds)
        context = openmm.Context(system, integrator, selected_platform)
        context.setPositions(pdb.positions)
        state = context.getState(getEnergy=True)
        energy = state.getPotentialEnergy()
        energy_kj_mol = energy.value_in_unit(kilojoules_per_mole)
        print(f"Initial potential energy (kJ/mol): {energy_kj_mol}")
        print(f"Initial potential energy finite: {math.isfinite(energy_kj_mol)}")
        del context
        del integrator
        return 0 if math.isfinite(energy_kj_mol) else 1
    except Exception as exc:
        print(f"FAILED: OpenMM read test failed: {exc!r}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
