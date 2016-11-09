#!/bin/sh

folder=../../files/2igd/
basename='2igd' 

phenix.AmberPrep $folder/$basename.pdb
phenix.refine 4phenix_$basename.pdb\
              $folder/$basename.mtz\
              amber.topology_file_name=4amber_$basename.prmtop\
              amber.coordinate_file_name=4amber_$basename.rst7\
              use_amber=True \
              wxc_scale=0.025 \
              --overwrite \
              amber.bellymask=':WAT' \
              refinement.main.number_of_macro_cycles=1
