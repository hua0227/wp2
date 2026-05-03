# Ligand-Bound MD Preflight Summary

## Purpose

This preflight checks ligand parameter file readability, staged input consistency, atom counts, CGenFF naming for docked poses, and visual-only complex previews before any ligand-bound MD setup.

No MD, minimization, receptor regeneration, docking rerun, `.str` modification, PSF generation, or original apo system modification was performed.

Complex preview PDB files are only for PyMOL visualization. They are not simulation-ready topologies and must not be used directly as proof that OpenMM ligand-bound MD setup is complete.

## CharmmParameterSet Load Check

- Successful ligand `.str` loads: 10
- Failed ligand `.str` loads: NONE

## Atom Count Consistency

- Ligands with atom counts consistent in both systems: NONE

## Named Pose Generation

- Ligands with named poses generated for both systems: NONE

## Complex Preview Generation

- Ligands with complex previews generated for both systems: NONE

## High-Penalty Ligands Requiring Caution

- 8.1, 12.2, 4.3, 8.3

High penalty does not automatically exclude a ligand, but it flags parameters for careful manual review or optimization before MD.

## Recommended First-Round OpenMM Read-Test Ligands

- 14.2, 19.1, 2.3, 6.2, 17.2

These recommendations prioritize moderate-penalty ligands requested for first-round parameter/read testing and successful `.str` loading.
Because atom-count matching failed for the current docked poses, complex-level OpenMM read tests still require a corrected ligand pose/topology alignment first.
