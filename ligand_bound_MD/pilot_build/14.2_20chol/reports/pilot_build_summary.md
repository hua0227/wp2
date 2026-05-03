# Ligand-Bound Topology Pilot Build Summary

## Purpose

This pilot continues the ligand-bound topology build for one case only: ligand 14.2 / resname L002 / 20chol. The goal was to run the VMD/psfgen build, then run an OpenMM read test that creates a Context and evaluates the initial potential energy without minimization or MD.

## Inputs

- VMD executable: `D:\software\VMD2\vmd.exe`
- OpenMM Python: `C:\Users\14566\miniconda3\envs\trkb_openmm\python.exe`
- Apo PSF: `C:\TRKB_WP2\TRKB_20CHOL\step5_assembly.psf`
- Apo short-MD PDB: `C:\TRKB_WP2\TRKB_20CHOL\openmm_short_md_output\short_md_final.pdb`
- GRO box source: `C:\TRKB_WP2\TRKB_20CHOL\gromacs\step5_input.gro`
- Ligand stream file: `C:\TRKB_WP2\ligand_bound_MD\cgenff_parameterization\returned_str\L002_14.2.str`
- Ligand psfgen-ready PDB: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\outputs\L002_14.2_20chol_psfgen_ready.pdb`

## Coordinate Staging

Direct `readpsf` with `short_md_final.pdb` fails because the short-MD PDB uses PDB chain/residue naming that does not match the CHARMM-GUI PSF segment/residue identity. To keep the post-short-MD coordinates without altering the original apo files, a pilot-local coordinate PDB was generated:

`C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\outputs\TRKB_20chol_shortmd_psfgen_ready.pdb`

This file uses atom identity/segid/resid columns from `step5_assembly.pdb` and xyz coordinates from `short_md_final.pdb`.

## Clash Precheck

- Severe clashes, ligand heavy to apo heavy < 1.0 A: 9
- Close contacts, ligand heavy to apo heavy < 1.5 A: 17
- Severe clash categories: lipid 9; protein 0; water 0; ion 0; cholesterol 0; unknown 0
- Close contact categories: lipid 16; cholesterol 1; protein 0; water 0; ion 0; unknown 0

Detailed reports:
- `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\reports\clash_precheck.csv`
- `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\reports\clash_precheck.md`

## VMD / psfgen Result

- VMD found: YES
- VMD path: `D:\software\VMD2\vmd.exe`
- psfgen script: `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\scripts\build_ligand_bound_psfgen.tcl`
- psfgen command exit code: 0
- Ligand-bound PSF generated: YES
- Ligand-bound PDB generated: YES
- Coordinates guessed by psfgen: 0 atoms

Outputs:
- `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\outputs\TRKB_20chol_L002_14.2_bound.psf`
- `C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\outputs\TRKB_20chol_L002_14.2_bound.pdb`

psfgen log:
`C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\logs\psfgen_build.log`

The psfgen log contains ERROR/WARNING lines. The most relevant ones are parser limitations while reading broad CHARMM stream/topology content, including unrecognized `set/if` statements and unsupported `LONEPAIR` / some internal-coordinate statements in stream files. Despite these messages, L002 was recognized, the LIG segment was built, ligand coordinates were assigned, 0 atoms needed guessed coordinates, and PSF/PDB files were written.

## OpenMM Read Test

- OpenMM read test executed: YES
- OpenMM read test success: YES
- OpenMM version environment: `trkb_openmm`
- Available platforms printed: Reference, CPU, OpenCL
- Selected platform: OpenCL
- CUDA used: NO
- `CharmmParameterSet` loaded: YES
- `createSystem` succeeded: YES
- Context creation succeeded: YES
- Missing parameters: no missing-parameter error observed
- Initial potential energy finite: YES
- Initial potential energy: `19539250927026.24 kJ/mol`

OpenMM log:
`C:\TRKB_WP2\ligand_bound_MD\pilot_build\14.2_20chol\logs\openmm_readtest.log`

The OpenMM log contains non-fatal warnings about replacing duplicate dihedral definitions. It does not show missing parameter errors. The energy is finite but very large, consistent with the pre-existing ligand-lipid close contacts/clashes and the fact that this stage intentionally does not relax the system.

## Explicit Non-Actions

- No MD was run.
- No energy minimization was run.
- The original apo PSF/PDB/GRO files were not modified.
- `C:\TRKB_WP2\TRKB_20CHOL\step5_assembly.psf` was not changed.
- `C:\TRKB_WP2\ligand_bound_MD\cgenff_parameterization\returned_str\L002_14.2.str` was not changed.
- Docking and receptor generation were not rerun.
- No Top10 batch processing was performed.

## Next Suggestions

1. Inspect the 9 severe ligand-lipid clashes and 17 close contacts before any downstream relaxation.
2. Consider reducing psfgen stream-parser noise by using only topology files required for this specific apo composition, while keeping parameter loading complete in OpenMM.
3. Use the generated PSF/PDB only as a topology/read-test pilot output. Any simulation-preparation step should explicitly handle clash resolution and relaxation in a separate stage.
