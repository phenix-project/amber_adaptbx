from __future__ import division
from libtbx.utils import Sorry
try:
  import sander
  import sanderles
except ImportError as e:
  sander = None
  sanderles = None

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
  order_file_name = None
    .type = path
    .help = A file that maps amber atom numbers to phenix atom numbers.
    .style = bold input_file
  wxc_factor = 1.
    .type = float
    .style = hidden
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
  qmmask = ''
    .type = str
    .style = hidden
    .help = If given, turn on QM/MM with the given mask
  qmcharge = 0
    .type = int
    .style = hidden
    .help = Charge of the QM/MM region
  netcdf_trajectory_file_name = ''
    .type = str
    .style = hidden
    .help = If given, turn on writing netcdf trajectory
  print_amber_energies = False
    .type = bool
    .help = Print details of Amber energies during refinement
"""
obsoleted = '''
  automatic_wxc_scale = False
    .type = bool
    .style = hidden
    .help = Use the ratio of the restraints gradient norm and the Amber \
            gradient norm to set wxc_scale
'''
