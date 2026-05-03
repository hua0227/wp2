from pathlib import Path

from openmm import unit, Platform, LangevinMiddleIntegrator
from openmm.app import (
    CharmmPsfFile,
    CharmmCrdFile,
    CharmmParameterSet,
    Simulation,
    PME,
    HBonds,
    PDBFile,
    StateDataReporter,
)

root = Path(".")

psf_file = root / "step5_assembly.psf"
crd_file = root / "step5_assembly.crd"
gro_file = root / "gromacs" / "step5_input.gro"
toppar_dir = root / "toppar"

output_dir = root / "openmm_test_output"
output_dir.mkdir(exist_ok=True)


def read_gro_box_lengths(gro_path):
    lines = [line.strip() for line in gro_path.read_text().splitlines() if line.strip()]
    last = lines[-1].split()
    values = [float(x) for x in last]
    if len(values) < 3:
        raise ValueError("Could not read box lengths from .gro file.")
    return (
        values[0] * unit.nanometer,
        values[1] * unit.nanometer,
        values[2] * unit.nanometer,
    )


print("[1] Loading structure...")
psf = CharmmPsfFile(str(psf_file))
crd = CharmmCrdFile(str(crd_file))

a, b, c = read_gro_box_lengths(gro_file)
psf.setBox(a, b, c)
print("Box:", a, b, c)
print("Atoms:", len(crd.positions))

print("[2] Loading parameters...")
param_files = [
    toppar_dir / "top_all36_prot.rtf",
    toppar_dir / "top_all36_lipid.rtf",
    toppar_dir / "top_all36_cgenff.rtf",

    toppar_dir / "par_all36m_prot.prm",
    toppar_dir / "par_all36_lipid.prm",
    toppar_dir / "par_all36_cgenff.prm",

    toppar_dir / "toppar_water_ions.str",
    toppar_dir / "toppar_all36_lipid_cholesterol.str",
]
params = CharmmParameterSet(*[str(p) for p in param_files])

print("[3] Creating system...")
system = psf.createSystem(
    params,
    nonbondedMethod=PME,
    nonbondedCutoff=1.2 * unit.nanometer,
    constraints=HBonds,
)

integrator = LangevinMiddleIntegrator(
    310 * unit.kelvin,
    1.0 / unit.picosecond,
    0.002 * unit.picoseconds,
)

platform_names = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]

if "OpenCL" in platform_names:
    platform = Platform.getPlatformByName("OpenCL")
    properties = {"Precision": "mixed"}
    print("[4] Using OpenCL platform.")
    simulation = Simulation(psf.topology, system, integrator, platform, properties)
elif "CPU" in platform_names:
    platform = Platform.getPlatformByName("CPU")
    print("[4] Using CPU platform.")
    simulation = Simulation(psf.topology, system, integrator, platform)
else:
    platform = Platform.getPlatformByName("Reference")
    print("[4] Using Reference platform.")
    simulation = Simulation(psf.topology, system, integrator, platform)

simulation.context.setPositions(crd.positions)

print("[5] Checking initial energy...")
state = simulation.context.getState(getEnergy=True)
print("Initial potential energy:", state.getPotentialEnergy())

simulation.reporters.append(
    StateDataReporter(
        str(output_dir / "minimize_log.txt"),
        1,
        step=True,
        potentialEnergy=True,
        temperature=True,
        progress=False,
        speed=True,
        totalSteps=1,
    )
)

print("[6] Minimizing energy...")
simulation.minimizeEnergy(maxIterations=500)

state = simulation.context.getState(getEnergy=True, getPositions=True)
print("Final potential energy:", state.getPotentialEnergy())

with open(output_dir / "minimized.pdb", "w") as f:
    PDBFile.writeFile(psf.topology, state.getPositions(), f)

print("[7] Done.")
print("Output written to:", output_dir)
print("Check file:", output_dir / "minimized.pdb")