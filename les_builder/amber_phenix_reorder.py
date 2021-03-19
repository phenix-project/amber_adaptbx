#!/usr/bin/env phenix.python

from __future__ import print_function
from iotbx import pdb
import numpy as np
from collections import OrderedDict

__all__ = ['initialize_order_converter',
           'get_indices_convert_dict_from_array',
           'reorder_coords_phenix_to_amber',
           'reorder_force_amber_to_phenix']

"""
# Simple example to show how to map phenix's coordinates to amber's coordinates
# Note: If you are using IPython, you can copy whole code block and run them

>>> import numpy as np
>>> from scitbx.array_family import flex
>>> phenix_coords = flex.vec3_double([[0., 2., 3.], 
...                                   [4., 5., 6.],
...                                   [10., 12., 13.],
...                                   [7., 8., 9.]])
>>> amber_coords = np.array([[10., 12., 13.],
...                          [7., 8., 9],
...                          [4., 5., 6.],
...                          [0., 2., 3.]])

>>> # Aim: creat coordinates map (amber->phenix, phenix->amber)
>>> from amber_adaptbx.amber_phenix_reorder import get_indices_convert_dict_from_array
>>> order_converter = get_indices_convert_dict_from_array(phenix_coords, amber_coords)
>>> order_converter
{'a2p': array([3, 2, 0, 1]), 'p2a': array([2, 3, 1, 0])}

>>> # Test
>>> # convert phenix_coords to amber_coordinates
>>> from amber_adaptbx.amber_phenix_reorder import reorder_coords_phenix_to_amber
>>> reorder_coords_phenix_to_amber(phenix_coords, order_converter['p2a'])
array([[ 10.,  12.,  13.],
       [  7.,   8.,   9.],
       [  4.,   5.,   6.],
       [  0.,   2.,   3.]])
# which is equal to origial amber_coords
# so on for force conversion
"""

def initialize_order_converter(self):
  '''make order_converter for geometry_manager class

  Parameters
  ----------
  self : instance of geometry_manager class
  '''
  geometry_manager = self.__class__
  # only compute the order_converter with the original sites_cart
  # subsequent geometry_manager objects will share the same `order_converter`
  # TODO: refactor
  if self.amber_structs.order_map_file_name is not None and self.amber_structs.order_converter is not None:
    geometry_manager.order_converter = self.amber_structs.order_converter
  else:
    #
    # this functionality has been moved to AmberPrep and should be removed
    # at least for the Amber Geom. Restraints Manager.
    #
    # compute reorder map based on original ASU pdb and unitcell rst7 (and parm7) files
    asu_n_atoms = len(self.sites_cart)
    n_models = int(len(self.amber_structs.parm.coordinates) / asu_n_atoms)
    geometry_manager.order_converter = get_indices_convert_dict_from_array(self.sites_cart,
            self.amber_structs.parm.coordinates[:asu_n_atoms])

    # asu
    asu_p2a_indices = geometry_manager.order_converter['p2a']
    asu_a2p_indices = geometry_manager.order_converter['a2p']

    # unitcell
    uc_p2a_indices = []
    uc_a2p_indices = []

    for index in range(n_models):
      # need to increase atom indices
      uc_p2a_indices.extend((asu_p2a_indices + index * asu_n_atoms).tolist())
      uc_a2p_indices.extend((asu_a2p_indices + index * asu_n_atoms).tolist())

    # extend array for unitcell
    geometry_manager.order_converter['p2a'] = np.array(uc_p2a_indices)
    geometry_manager.order_converter['a2p'] = np.array(uc_a2p_indices)

    # save to disk for debugging
    order_converter = geometry_manager.order_converter
    saved_arr = np.array([order_converter['a2p'], order_converter['p2a']], dtype='i4')
    # 1st column: amber -> phenix
    # 2nd column: phenix -> amber
    filename = 'amber_phenix_atom_order_map.txt'
    np.savetxt(filename, saved_arr.transpose(), fmt='%5d')
    # end debugging

def reorder_force_amber_to_phenix(frc, new_indices):
  """Convert amber force to phenix's order
  
  Parameters
  ----------
  frc : 1D-array force, calculated by amber
  new_indices : 1D array
      used for remapping from amber order to phenix order

  Returns
  -------
  reordered numpy 1D-array forces for phenix
  """
  frc = np.asarray(frc)
  n_atoms = int(frc.shape[0]/3)
  new_frc = frc.reshape(n_atoms, 3)[new_indices]
  return new_frc.flatten()

def reorder_coords_phenix_to_amber(coords, new_indices):
  """Convert phenix site_cart to amber format
  
  Paramters
  ---------
  coords: flex.vec3_double
  new_indices : 1D numpy array

  Returns
  -------
  reordered numpy 2D-array coordinates for sander
  """
  coords_2d = np.asarray(coords)
  return coords_2d[new_indices]

def get_indices_convert_dict(fn):
  """return a dict(a2p=arr0, p2a=arr1).

  Key "a2p" means converting amber order to phenix order
  Key "p2a" means converting phenix order to amber order

  Parameters
  ----------
  fn : str, pdb filename
  """
  pdb_inp = pdb.input(file_name=fn)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  
  newids = OrderedDict((atom.id_str(), idx) for (idx, atom) in enumerate(pdb_hierarchy.atoms()))
  oldids= OrderedDict((atom.id_str(), idx) for (idx, atom) in enumerate(pdb_inp.atoms()))
  
  return {'p2a': np.array([newids[atom.id_str()] for atom in pdb_inp.atoms()]),
          'a2p': np.array([oldids[atom.id_str()] for atom in pdb_hierarchy.atoms()])}



if __name__ == '__main__':
  import parmed as pmd
  from numpy.testing import assert_almost_equal as aa_eq

  root = './test_4lzt_les/'
  fn, tn, rst7 = [root + filename for filename in ["4lztabH.pdb", "4lztab.parm7", "4lztabH.rst7"]]

  pdb_input = pdb.input(fn)
  hi = pdb_input.construct_hierarchy()
  parm = pmd.load_file(tn, rst7)
  
  p2a_indices = get_indices_convert_dict(fn)['p2a']
  a2p_indices = get_indices_convert_dict(fn)['a2p']
  
  a0 = reorder_coords_phenix_to_amber(hi.atoms().extract_xyz(), p2a_indices)
  print(a0.shape)
  aa_eq(a0, parm.coordinates.flatten())
  print("OK")
