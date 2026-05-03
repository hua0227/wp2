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
toppar_dir = root / "toppar"

print("[1] Checking files...")
for f in [psf_file, crd_file, toppar_dir]:
    print(f, "OK" if f.exists() else "MISSING")

if not psf_file.exists():
    raise FileNotFoundError("Missing step5_assembly.psf")
if not crd_file.exists():
    raise FileNotFoundError("Missing step5_assembly.crd")
if not toppar_dir.exists():
    raise FileNotFoundError("Missing toppar directory")

param_files = []
for pattern in ["*.rtf", "*.top", "*.prm", "*.par", "*.str", "*.inp"]:
    param_files.extend(sorted(toppar_dir.rglob(pattern)))

print("\n[2] Parameter files found:")
for f in param_files:
    print(" ", f)

if len(param_files) == 0:
    raise RuntimeError("No CHARMM parameter files found in toppar.")

print("\n[3] Loading PSF and CRD...")
psf = CharmmPsfFile(str(psf_file))
crd = CharmmCrdFile(str(crd_file))

print("Atoms in PSF:", len(list(psf.topology.atoms())))
print("Positions in CRD:", len(crd.positions))

print("\n[4] Loading CHARMM parameters...")
params = CharmmParameterSet(*[str(f) for f in param_files])

print("\n[5] Setting periodic box...")
# 尝试从 CRD 读取盒子信息；如果失败，再手动处理
try:
    psf.setBox(*crd.boxLengths)
    print("Box lengths:", crd.boxLengths)
except Exception as e:
    print("Could not set box from CRD:", e)
    print("If this fails later, we will use the GROMACS .gro box instead.")

print("\n[6] Creating OpenMM System...")
system = psf.createSystem(
    params,
    nonbondedMethod=PME,
    nonbondedCutoff=1.2 * unit.nanometer,
    constraints=HBonds,
)

print("Number of particles:", system.getNumParticles())

print("\n[7] Available OpenMM platforms:")
for i in range(Platform.getNumPlatforms()):
    print("-", Platform.getPlatform(i).getName())

print("\nSUCCESS: CHARMM files can be read by OpenMM.")