# Pilot VMD/psfgen build script for ligand 6.2 / L008 / 20chol.
package require psfgen
resetpsf

set toppar_dir [file normalize "C:/TRKB_WP2/TRKB_20CHOL/toppar" ]
set ligand_str [file normalize "C:/TRKB_WP2/ligand_bound_MD/cgenff_parameterization/returned_str/L008_6.2.str" ]
set apo_psf [file normalize "C:/TRKB_WP2/TRKB_20CHOL/step5_assembly.psf" ]
set apo_psfgen_pdb [file normalize "C:/TRKB_WP2/ligand_bound_MD/pilot_build/6.2_20chol/outputs/TRKB_20chol_psfgen_ready_apo.pdb" ]
set ligand_pdb [file normalize "C:/TRKB_WP2/ligand_bound_MD/pilot_build/6.2_20chol/outputs/L008_6.2_20chol_psfgen_ready.pdb" ]
set out_psf [file normalize "C:/TRKB_WP2/ligand_bound_MD/pilot_build/6.2_20chol/outputs/TRKB_20chol_L008_6.2_bound.psf" ]
set out_pdb [file normalize "C:/TRKB_WP2/ligand_bound_MD/pilot_build/6.2_20chol/outputs/TRKB_20chol_L008_6.2_bound.pdb" ]

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
    residue 1 L008
}
coordpdb $ligand_pdb LIG

guesscoord
writepsf $out_psf
writepdb $out_pdb
puts "Wrote ligand-bound PSF: $out_psf"
puts "Wrote ligand-bound PDB: $out_pdb"
quit
