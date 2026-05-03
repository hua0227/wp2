# CGenFF / ParamChem Upload Package

Submit the `.mol2` files in `upload_mol2` to CGenFF/ParamChem.

Each ligand needs a returned `.str` file. Recommended returned file names:

| Rank | Ligand ID | Residue name | Upload mol2 | Expected returned str |
|---:|---|---|---|---|
| 1 | 8.1 | L001 | `L001_8.1.mol2` | `L001_8.1.str` |
| 2 | 14.2 | L002 | `L002_14.2.mol2` | `L002_14.2.str` |
| 3 | 19.1 | L003 | `L003_19.1.mol2` | `L003_19.1.str` |
| 4 | 9.2 | L004 | `L004_9.2.mol2` | `L004_9.2.str` |
| 5 | 12.2 | L005 | `L005_12.2.mol2` | `L005_12.2.str` |
| 6 | 2.3 | L006 | `L006_2.3.mol2` | `L006_2.3.str` |
| 7 | 4.3 | L007 | `L007_4.3.mol2` | `L007_4.3.str` |
| 8 | 6.2 | L008 | `L008_6.2.mol2` | `L008_6.2.str` |
| 9 | 8.3 | L009 | `L009_8.3.mol2` | `L009_8.3.str` |
| 10 | 17.2 | L010 | `L010_17.2.mol2` | `L010_17.2.str` |

Place returned `.str` files in:

`C:\TRKB_WP2\ligand_bound_MD\cgenff_parameterization\returned_str`

Do not manually fabricate `.str` files.

After ParamChem/CGenFF returns the `.str` files, inspect CGenFF penalties carefully.
Pay special attention to high-penalty bond, angle, and dihedral parameters before any MD setup.
