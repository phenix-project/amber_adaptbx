from __future__ import division
from libtbx.utils import Sorry
try:
  import sander
  import sanderles
except:
  raise Sorry('Unable to import "sander". Check that $AMBERHOME is set correctly to the Amber directory.')
from amber_adaptbx.utils import (
    expand_coord_to_unit_cell,
    get_amber_structs,
    bond_rmsd,
    bond_rmsZ,
    angle_rmsZ,
    angle_rmsd,
    is_prmtop_LES,
    collapse_grad_to_asu,
    check_file,
    print_sites_cart,
)
from amber_adaptbx.geometry import geometry_manager
from amber_adaptbx.energy import energies, SanderStruct

master_phil_str = """
  use_amber = False
    .type = bool
    .help = Use Amber for all the gradients in refinement
    .short_caption = Enable Amber
  topology_file_name = None
    .type = path
    .help = A topology file needed by Amber. Can be generated using phenix.AmberPrep.
    .style = bold input_file
  coordinate_file_name = None
    .type = path
    .help = A coordinate file needed by Amber. Can be generated using phenix.AmberPrep.
    .style = bold input_file
  order_mapping_file_name = None
    .type = path
    .help = A pickled file that maps amber's coordinates to phenix's coordinates. If given, reuse this map.
    .style = bold input_file
  wxc_factor = .1
    .type = float
    .style = hidden
  automatic_wxc_scale = False
    .type = bool
    .style = hidden
    .help = Use the ratio of the restraints gradident norm and the Amber \
            gradient norm to set wxc_scale
  restraint_wt = 0.
    .type = float
    .style = hidden
  restraintmask = ''
    .type = str
    .style = hidden
  reference_file_name = ''
    .type = str
    .style = hidden
  bellymask = ''
    .type = str
    .style = hidden
    .help = If given, turn on belly in sander
  netcdf_trajectory_file_name = ''
    .type = str
    .style = hidden
    .help = If given, turn on writing netcdf trajectory
  print_amber_energies = False
    .type = bool
    .help = Print details of Amber energies during refinement
"""
