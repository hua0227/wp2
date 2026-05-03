# CGenFF Penalty Review

This review marks risk based on CGenFF/ParamChem penalty values only. No ligand is automatically discarded.
High penalties do not necessarily mean a ligand cannot be used, but they do require manual inspection and possibly parameter optimization before MD.

## Penalty Summary

| Ligand | Resname | Param | Charge | Max | Level | Manual review | Notes |
|---|---|---:|---:|---:|---|---:|---|
| 8.1 | L001 | 54.8 | 22.945 | 54.8 | high | 1 | high CGenFF analogy penalty; inspect or optimize parameters before MD |
| 14.2 | L002 | 19.5 | 22.916 | 22.916 | moderate | 0 | basic validation recommended |
| 19.1 | L003 | 19.5 | 23.991 | 23.991 | moderate | 0 | basic validation recommended |
| 9.2 | L004 | 11.6 | 21.149 | 21.149 | moderate | 0 | basic validation recommended |
| 12.2 | L005 | 54.8 | 24.018 | 54.8 | high | 1 | high CGenFF analogy penalty; inspect or optimize parameters before MD |
| 2.3 | L006 | 11 | 19.747 | 19.747 | moderate | 0 | basic validation recommended |
| 4.3 | L007 | 54.8 | 31.016 | 54.8 | high | 1 | high CGenFF analogy penalty; inspect or optimize parameters before MD |
| 6.2 | L008 | 21 | 24.632 | 24.632 | moderate | 0 | basic validation recommended |
| 8.3 | L009 | 54.8 | 22.945 | 54.8 | high | 1 | high CGenFF analogy penalty; inspect or optimize parameters before MD |
| 17.2 | L010 | 20.6 | 28.851 | 28.851 | moderate | 0 | basic validation recommended |

## Preliminary OpenMM Read-Test Candidates

- 14.2 (L002): moderate, max penalty 22.916
- 19.1 (L003): moderate, max penalty 23.991
- 9.2 (L004): moderate, max penalty 21.149
- 2.3 (L006): moderate, max penalty 19.747
- 6.2 (L008): moderate, max penalty 24.632
- 17.2 (L010): moderate, max penalty 28.851

## Ligands Requiring Caution

- 8.1 (L001): high, max penalty 54.8
- 12.2 (L005): high, max penalty 54.8
- 4.3 (L007): high, max penalty 54.8
- 8.3 (L009): high, max penalty 54.8

## Interpretation Notes

- `low`: max penalty below 10.
- `moderate`: max penalty from 10 to below 50.
- `high`: max penalty from 50 to below 100.
- `very_high`: max penalty 100 or higher.
- `unknown`: penalty could not be parsed or `.str` was missing.
- High penalty is a risk marker, not an automatic exclusion rule.
