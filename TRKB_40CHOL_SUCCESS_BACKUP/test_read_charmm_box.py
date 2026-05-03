from pathlib import Path

from openmm import unit, Platform
from openmm.app import (
    CharmmPsfFile,
    CharmmCrdFile,
    CharmmParameterSet,
    PME,
    HBonds,
)

root = Path(".")

psf_file = root / "step5_assembly.psf"
crd_file = root / "step5_assembly.crd"
gro_file = root / "gromacs" / "step5_input.gro"
toppar_dir = root / "toppar"

print("[1] Checking files...")
for f in [psf_file, crd_file, gro_file, toppar_dir]:
    print(f, "OK" if f.exists() else "MISSING")

if not psf_file.exists():
    raise FileNotFoundError("Missing step5_assembly.psf")
if not crd_file.exists():
    raise FileNotFoundError("Missing step5_assembly.crd")
if not gro_file.exists():
    raise FileNotFoundError("Missing gromacs/step5_input.gro")
if not toppar_dir.exists():
    raise FileNotFoundError("Missing toppar directory")


def read_gro_box_lengths(gro_path):
    """
    Read the last line of a GROMACS .gro file.
    For a rectangular CHARMM-GUI membrane box, the first three values are x, y, z in nm.
    """
    lines = [line.strip() for line in gro_path.read_text().splitlines() if line.strip()]
    last = lines[-1].split()
    values = [float(x) for x in last]

    if len(values) < 3:
        raise ValueError("Could not read box lengths from the last line of .gro file.")

    a = values[0] * unit.nanometer
    b = values[1] * unit.nanometer
    c = values[2] * unit.nanometer
    return a, b, c


print("\n[2] Loading PSF and CRD...")
psf = CharmmPsfFile(str(psf_file))
crd = CharmmCrdFile(str(crd_file))

print("Atoms in PSF:", len(list(psf.topology.atoms())))
print("Positions in CRD:", len(crd.positions))

print("\n[3] Reading box from GROMACS .gro...")
a, b, c = read_gro_box_lengths(gro_file)
print("Box lengths:", a, b, c)

psf.setBox(a, b, c)
print("Periodic box set successfully.")

print("\n[4] Loading CHARMM parameters...")

# 优先只读取本体系需要的核心 CHARMM36 参数，避免加载过多无关 str 文件产生大量重复参数 warning
param_candidates = [
    toppar_dir / "top_all36_prot.rtf",
    toppar_dir / "top_all36_lipid.rtf",
    toppar_dir / "top_all36_cgenff.rtf",

    toppar_dir / "par_all36m_prot.prm",
    toppar_dir / "par_all36_lipid.prm",
    toppar_dir / "par_all36_cgenff.prm",

    toppar_dir / "toppar_water_ions.str",
    toppar_dir / "toppar_all36_lipid_cholesterol.str",
]

param_files = [str(p) for p in param_candidates if p.exists()]

print("Parameter files used:")
for p in param_files:
    print(" ", p)

if len(param_files) == 0:
    raise RuntimeError("No CHARMM parameter files found.")

params = CharmmParameterSet(*param_files)

print("\n[5] Creating OpenMM System...")
system = psf.createSystem(
    params,
    nonbondedMethod=PME,
    nonbondedCutoff=1.2 * unit.nanometer,
    constraints=HBonds,
)

print("Number of particles:", system.getNumParticles())

print("\n[6] Available OpenMM platforms:")
for i in range(Platform.getNumPlatforms()):
    print("-", Platform.getPlatform(i).getName())

print("\nSUCCESS: CHARMM PSF/CRD + GROMACS box can be read by OpenMM.")