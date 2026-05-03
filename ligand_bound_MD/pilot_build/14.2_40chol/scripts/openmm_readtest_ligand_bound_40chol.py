from __future__ import annotations

import math

import openmm
from openmm import Platform, VerletIntegrator
from openmm.app import CharmmParameterSet, CharmmPsfFile, NoCutoff, PDBFile, PME
from openmm.unit import kilojoules_per_mole, nanometer, picoseconds

from common_40chol import BOUND_PDB, BOUND_PSF, GRO_BOX, PARAMETER_FILES, choose_opencl_platform, ensure_dirs, set_psf_box_from_gro


def main() -> int:
    ensure_dirs()
    print("OpenMM ligand-bound read test for pilot 14.2 / L002 / 40chol")
    print("OpenCL only. No minimization and no MD will be run.")
    print(f"Input PSF: {BOUND_PSF}")
    print(f"Input PDB: {BOUND_PDB}")
    print(f"GRO box: {GRO_BOX}")
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
        set_psf_box_from_gro(psf)
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
        except Exception as exc:
            print(f"PME createSystem failed: {exc!r}")
            print("Trying NoCutoff only to isolate topology/parameter readability.")
            system = psf.createSystem(params, nonbondedMethod=NoCutoff, constraints=None)
        platform = choose_opencl_platform(Platform)
        print(f"Selected OpenMM platform: {platform.getName()}")
        integrator = VerletIntegrator(0.001 * picoseconds)
        context = openmm.Context(system, integrator, platform)
        context.setPositions(pdb.positions)
        state = context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
        print(f"Initial potential energy (kJ/mol): {energy}")
        print(f"Initial potential energy finite: {math.isfinite(energy)}")
        del context
        del integrator
        return 0 if math.isfinite(energy) else 1
    except Exception as exc:
        print(f"FAILED: OpenMM read test failed: {exc!r}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
