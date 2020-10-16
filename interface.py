import sys
import numpy as np
import amber_adaptbx
from libtbx.utils import Sorry

# TODO: create writer in another place? where?
from parmed.amber.netcdffiles import NetCDFTraj

def get_amber_struct_object(params):
  amber_params = params.amber
  ridingH = True
  if getattr(params, "hydrogens", False):  # ensemble refinement does not have this
    if params.hydrogens.refine in ['riding', 'Auto']:
      ridingH = True
    elif params.hydrogens.refine in ['individual']:
      ridingH = False
    else:
      raise Sorry("Hydrogens.refine parameter '%s' unknown!"
                  % params.hydrogens.refine)
  amber_structs = amber_adaptbx.SanderStruct(
      parm_file_name=amber_params.topology_file_name,
      rst_file_name=amber_params.coordinate_file_name,
      ridingH=ridingH,
  )

  if amber_params.bellymask:
    try:
      amber_structs.inp.ibelly = 1
      amber_structs.inp.bellymask = amber_params.bellymask
    except AttributeError:
      raise AttributeError(
          'Setting bellymask for pysander does not work with AmberTools <= 16')

  if amber_params.restraint_wt > 0.:
    try:
      amber_structs.inp.ntr = 1
      amber_structs.inp.restraint_wt = amber_params.restraint_wt
      amber_structs.inp.restraintmask = amber_params.restraintmask
      if amber_params.reference_file_name:
        amber_structs.inp.refc = amber_params.reference_file_name
      else:
        amber_structs.inp.refc = amber_params.coordinate_file_name
    except AttributeError as e:
      raise AttributeError(
          'Setting amber restraint for pysander does not work with AmberTools <= 16')

  amber_structs.qm_inp=amber_structs.sander_engine.QmInputOptions()
  if amber_params.qmmask:
     amber_structs.inp.ifqnt=1
     amber_structs.qm_inp.qm_theory = 'PM6'
     amber_structs.qm_inp.qmmask = amber_params.qmmask
     amber_structs.qm_inp.qmcut = 8.0
     amber_structs.qm_inp.qmcharge= amber_params.qmcharge
     amber_structs.qm_inp.diag_routine=1
  else:
     amber_structs.inp.ifqnt=0

  if not amber_structs.sander_engine.is_setup():
    amber_structs.sander_engine.setup(
      amber_structs.parm,
      amber_structs.parm.coordinates,
      amber_structs.parm.box,
      amber_structs.inp,
      amber_structs.qm_inp
      )

  if amber_params.order_file_name is not None:
    order_map_file_name = amber_params.order_file_name
    amber_structs.order_map_file_name = order_map_file_name
    mapped_arr = np.loadtxt(order_map_file_name, dtype='i4').transpose()
    amber_structs.order_converter = dict(a2p=mapped_arr[0], p2a=mapped_arr[1])

  if amber_params.netcdf_trajectory_file_name:
    n_atoms = len(amber_structs.parm.atoms)
    amber_structs.writer = NetCDFTraj.open_new(amber_params.netcdf_trajectory_file_name,
             n_atoms, box=True, crds=True, frcs=False)
  return amber_structs, amber_structs.sander_engine
