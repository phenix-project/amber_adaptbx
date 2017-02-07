phenix.refine  \
   4phenix_$1.pdb $1.mtz \
   c_beta_restraints=False discard_psi_phi=False \
   strategy=individual_sites+individual_adp \
   refinement.main.number_of_macro_cycles=10 \
   topology_file_name=4amber_$1.prmtop \
   amber.coordinate_file_name=4amber_$1.rst7 \
   use_amber=True \
   prefix=amber serial=1 \
   write_geo=False --overwrite cdl=True | tee amber1.log 
