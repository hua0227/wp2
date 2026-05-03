# Targeted Clash Cleanup Summary

## Scope

- ligand_id: 14.2
- resname: L002
- system: 20chol
- This work created only pilot-local cleaned files under `targeted_clash_cleanup`.
- No original apo files were modified.
- The ligand `.str` file was not modified.
- Docking and receptor generation were not rerun.

## Original Failure Reason

- The short restrained MD pilot stopped at step 1 by safety rule.
- Step-1 temperature: 1797.49 K.
- Step-1 potential energy: 4.16529e+08 kJ/mol.
- No NaN was present in the final safe output.
- The ligand remained near TYR13/VAL17/SER20 and near two chains at the safety stop.
- Prior force diagnostics implicated LIG H9/H11 and nearby MEMB:94 POPC atoms.

## Detailed Clash Diagnosis

- Bound starting coordinates: 27 severe contacts <1.0 A and 81 close contacts <1.5 A.
- MEMB:94 POPC in bound coordinates: 22 severe and 63 close contacts.
- MEMB:94 POPC dominant severe source in bound coordinates: YES.
- LIG H9/H11 severe contacts with MEMB:94 POPC in bound coordinates: 1.
- Previous minimized coordinates: 0 severe contacts <1.0 A and 1 close contacts <1.5 A.
- Atom-pair classes for bound severe contacts: {'heavy-heavy': 9, 'H-H': 6, 'H-heavy': 5, 'heavy-H': 7}.

## Branch A: No Lipid Deletion

- Input: original bound PSF/PDB.
- POPC deleted: NO.
- Initial energy after constraints: 1.95686678524e+13 kJ/mol.
- Final energy: 1.95686678272e+13 kJ/mol.
- Energy finite: YES.
- Remaining severe contacts <1.0 A: 27.
- Remaining close contacts <1.5 A: 81.
- LIG H9/H11 POPC severe/close contacts: 1 / 3.
- Ligand min distance to TYR13/VAL17/SER20: 3.02229 A.
- Ligand near chains count: 2.
- Ligand RMSD relative to pre-minimized pose: 0.109691 A.
- Result: retained ligand position, but did not reduce the dominant POPC94 clash.

## Branch B: Remove MEMB:94 POPC

- Input: pilot-local VMD-cleaned PSF/PDB.
- POPC deleted: YES, exactly the whole MEMB:94 POPC molecule.
- Remaining MEMB:94 POPC atoms in cleaned PDB: 0.
- Initial energy after constraints: 2805010204.39 kJ/mol.
- Final energy: 212021.732669 kJ/mol.
- Energy finite: YES.
- Remaining severe contacts <1.0 A: 0.
- Remaining close contacts <1.5 A: 0.
- LIG H9/H11 POPC severe/close contacts: 0 / 0.
- Ligand min distance to TYR13/VAL17/SER20: 3.24676 A.
- Ligand near chains count: 2.
- Ligand RMSD relative to pre-minimized pose: 1.31314 A.
- Periodic minimum-image distances were used for branch retention checks.
- Result: removed the severe/close ligand contacts while retaining the ligand near the pocket and both chains.

## Optional Ultra-Conservative 1 ps Check

- Branch checked: Branch B only.
- Conditions: OpenCL, 0.1 fs timestep, 50 K, strong restraints, no trajectory, not production MD.
- Completed 1 ps: NO.
- Last finite step: 10 / 10000.
- Nonfinite energy detected: NO.
- Stop reason: temperature_above_200K.
- Final severe/close contacts at last finite step: 0 / 0.
- Ligand key min distance at last finite step: 3.24779 A.
- Ligand near chains at last finite step: 2.

## Recommendation

- Branch B is the correct structural cleanup branch, but it should not proceed to a longer or production MD run yet because the optional 1 ps check still triggered a temperature safety stop.
- Recommended next step: continue with Branch B only and perform additional constrained/round-trip relaxation or exact-coordinate handoff before another short restrained MD attempt.
- Branch A should be rejected because POPC94 clashes remain.
- No production MD was run.
