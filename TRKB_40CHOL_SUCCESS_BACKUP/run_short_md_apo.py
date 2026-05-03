from pathlib import Path

from openmm import unit, Platform, LangevinMiddleIntegrator, MonteCarloBarostat
from openmm.app import (
    CharmmPsfFile,
    CharmmCrdFile,
    CharmmParameterSet,
    Simulation,
    PME,
    HBonds,
    PDBFile,
    StateDataReporter,
    DCDReporter,
    CheckpointReporter,
)

root = Path(".")

psf_file = root / "step5_assembly.psf"
crd_file = root / "step5_assembly.crd"
gro_file = root / "gromacs" / "step5_input.gro"
toppar_dir = root / "toppar"

minimized_pdb = root / "openmm_test_output" / "minimized.pdb"

output_dir = root / "openmm_short_md_output"
output_dir.mkdir(exist_ok=True)


def read_gro_box_lengths(gro_path):
    lines = [line.strip() for line in gro_path.read_text().splitlines() if line.strip()]
    last = lines[-1].split()
    values = [float(x) for x in last]
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

# 短测试可以先加各向同性压力控制；正式膜蛋白建议后面改成更贴近 CHARMM-GUI 的膜系统平衡方案
system.addForce(MonteCarloBarostat(1.0 * unit.atmosphere, 310 * unit.kelvin, 50))

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

if minimized_pdb.exists():
    print("[5] Loading minimized positions...")
    pdb = PDBFile(str(minimized_pdb))
    simulation.context.setPositions(pdb.positions)
else:
    print("[5] Minimized PDB not found. Using original CRD positions.")
    simulation.context.setPositions(crd.positions)

simulation.context.setVelocitiesToTemperature(310 * unit.kelvin)

simulation.reporters.append(
    StateDataReporter(
        str(output_dir / "short_md_log.txt"),
        100,
        step=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        speed=True,
        separator="\t",
    )
)

simulation.reporters.append(
    DCDReporter(str(output_dir / "short_md.dcd"), 100)
)

simulation.reporters.append(
    CheckpointReporter(str(output_dir / "short_md.chk"), 500)
)

print("[6] Running 5000 steps short MD...")
simulation.step(5000)

state = simulation.context.getState(getPositions=True, getEnergy=True)
print("Final potential energy:", state.getPotentialEnergy())

with open(output_dir / "short_md_final.pdb", "w") as f:
    PDBFile.writeFile(psf.topology, state.getPositions(), f)

print("[7] Done.")
print("Outputs written to:", output_dir)