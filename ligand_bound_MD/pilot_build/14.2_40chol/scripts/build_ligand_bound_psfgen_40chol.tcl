# Pilot VMD/psfgen build script for ligand 14.2 / L002 / 40chol.
# This script appends a standalone ligand segment to the existing 40chol apo PSF/PDB.
# It does not modify original apo files, modify ligand .str, run minimization, or run MD.

package require psfgen
resetpsf

set toppar_dir [file normalize "C:/TRKB_WP2/TRKB_40CHOL/toppar"]
set ligand_str [file normalize "C:/TRKB_WP2/ligand_bound_MD/cgenff_parameterization/returned_str/L002_14.2.str"]
set apo_psf [file normalize "C:/TRKB_WP2/TRKB_40CHOL/step5_assembly.psf"]
set apo_psfgen_pdb [file normalize "C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_40chol/outputs/TRKB_40chol_psfgen_ready_apo.pdb"]
set ligand_pdb [file normalize "C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_40chol/outputs/L002_14.2_40chol_psfgen_ready.pdb"]
set out_psf [file normalize "C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_40chol/outputs/TRKB_40chol_L002_14.2_bound.psf"]
set out_pdb [file normalize "C:/TRKB_WP2/ligand_bound_MD/pilot_build/14.2_40chol/outputs/TRKB_40chol_L002_14.2_bound.pdb"]

topology [file join $toppar_dir "top_all36_prot.rtf"]
topology [file join $toppar_dir "top_all36_lipid.rtf"]
topology [file join $toppar_dir "top_all36_cgenff.rtf"]
topology [file join $toppar_dir "toppar_water_ions.str"]
topology [file join $toppar_dir "toppar_all36_lipid_cholesterol.str"]
topology $ligand_str

readpsf $apo_psf pdb $apo_psfgen_pdb

segment LIG {
    first none
    last none
    residue 1 L002
}
coordpdb $ligand_pdb LIG

guesscoord
writepsf $out_psf
writepdb $out_pdb

puts "Wrote ligand-bound PSF: $out_psf"
puts "Wrote ligand-bound PDB: $out_pdb"
quit
