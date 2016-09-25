#!/bin/sh

root_name=../../files/vAla3/vAla3
phenix.refine $root_name.pdb\
              $root_name.cif\
              $root_name.mtz\
              topology_file_name=$root_name.prmtop\
              coordinate_file_name=$root_name.rst7\
              use_amber=True \
              wxc_scale=0.025 \
              --overwrite \
              refinement.main.number_of_macro_cycles=1
