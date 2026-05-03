# Batch Gate Status Report

This report tracks the gate-based ligand-bound non-production pilot pipeline for 2.3 / L006, 6.2 / L008, and 17.2 / L010 across 20chol and 40chol.

- Pipeline order: preflight -> psfgen build -> OpenMM read test -> clash diagnosis -> targeted cleanup -> R2 minimization -> very-short probe -> 10 ps -> 50 ps -> 100 ps
- Any ligand-system stops at the first failing or unsafe stage.
- Other ligand-systems may continue independently.

## Final Status Table

| ligand_id | resname | system | final_stage | completed_10ps | completed_50ps | completed_100ps | stop_reason | error |
|---|---|---|---|---:|---:|---:|---|---|
| 2.3 | L006 | 20chol | completed_100ps | True | True | True | none |  |
| 2.3 | L006 | 40chol | completed_100ps | True | True | True | none |  |
| 6.2 | L008 | 20chol | completed_100ps | True | True | True | none |  |
| 6.2 | L008 | 40chol | completed_100ps | True | True | True | none |  |
| 17.2 | L010 | 20chol | completed_100ps | True | True | True | none |  |
| 17.2 | L010 | 40chol | completed_100ps | True | True | True | none |  |
