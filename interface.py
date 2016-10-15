import amber_adaptbx
import numpy as np
from libtbx.utils import Sorry


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

  amber_structs.sander_engine.setup(
      amber_structs.parm,
      amber_structs.rst.coordinates,
      amber_structs.rst.box,
      amber_structs.inp
  )

  if amber_params.order_mapping_file_name is not None:
    order_map_file_name = amber_params.order_mapping_file_name
    amber_structs.order_map_file_name = order_map_file_name
    mapped_arr = np.loadtxt(order_map_file_name, dtype='i4').transpose()
    amber_structs.order_converter = dict(a2p=mapped_arr[0], p2a=mapped_arr[1])

  return amber_structs, amber_structs.sander_engine
