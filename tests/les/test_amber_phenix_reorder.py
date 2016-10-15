# run with py.test
import os
import numpy as np
import pytest
from numpy.testing import assert_almost_equal as aa_eq
from numpy.testing import assert_equal as eq

from iotbx import pdb
from scitbx.array_family import flex

import parmed as pmd
from amber_adaptbx.les_builder.amber_phenix_reorder import (reorder_coords_phenix_to_amber,
        get_indices_convert_dict, get_indices_convert_dict_from_array)
from amber_adaptbx import expand_coord_to_unit_cell, SanderStruct
from amber_adaptbx.tests.utils import get_fn

def test_normal_pdb():
  fn = get_fn('4lzt/4lzt.not_les.pdb')
  tn = get_fn('4lzt/4lzt.parm7')
  rst7 = get_fn('4lzt/4lzt.rst7')

  pdb_input = pdb.input(fn)
  symm = pdb_input.crystal_symmetry()
  pdb_hc = pdb_input.construct_hierarchy()

  phenix_coords_uc = expand_coord_to_unit_cell(pdb_hc.atoms().extract_xyz(), symm)

  indices_dict = get_indices_convert_dict(fn)
  p2a_indices = indices_dict['p2a']
  a2p_indices = indices_dict['a2p']
  parm = pmd.load_file(tn, rst7)

  a0 = reorder_coords_phenix_to_amber(phenix_coords_uc, p2a_indices)

  # make sure phenix and amber read the same coordinates
  aa_eq(pdb_input.atoms().extract_xyz(), parm.coordinates, decimal=2)
  aa_eq(a0, parm.coordinates, decimal=2)

  # make sure the a2p_indices and p2a_indices are identical for non-LES
  indices_dict = get_indices_convert_dict_from_array(pdb_input.atoms().extract_xyz(), parm.coordinates)
  p2a_indices = indices_dict['p2a']
  a2p_indices = indices_dict['a2p']
  eq(a2p_indices, p2a_indices)

def test_converter():
  fn = get_fn('4lzt_pawel/4lztabH.pdb')
  tn = get_fn('4lzt_pawel/4lztab.parm7')
  rst7 = get_fn('4lzt_pawel/4lztabH.rst7')
  
  pdb_input = pdb.input(fn)
  parm = pmd.load_file(tn, rst7)

  symm = pdb_input.crystal_symmetry()
  pdb_hc = pdb_input.construct_hierarchy()
  
  phenix_coords_uc = expand_coord_to_unit_cell(pdb_hc.atoms().extract_xyz(), symm)

  indices_dict = get_indices_convert_dict(fn)
  p2a_indices = indices_dict['p2a']
  a2p_indices = indices_dict['a2p']
  # make sure two indices are not equal
  with pytest.raises(AssertionError):
      aa_eq(p2a_indices, a2p_indices)

  a0 = reorder_coords_phenix_to_amber(phenix_coords_uc, p2a_indices)
  aa_eq(pdb_input.atoms().extract_xyz(), parm.coordinates)
  aa_eq(a0, parm.coordinates)

def test_converter_by_array():
  fn = get_fn('4lzt_pawel/4lztabH.pdb')
  tn = get_fn('4lzt_pawel/4lztab.parm7')
  rst7 = get_fn('4lzt_pawel/4lztabH.rst7')
  
  pdb_input = pdb.input(fn)
  parm = pmd.load_file(tn, rst7)

  symm = pdb_input.crystal_symmetry()

  # reorder
  pdb_hc = pdb_input.construct_hierarchy()
  
  phenix_coords_uc = expand_coord_to_unit_cell(pdb_hc.atoms().extract_xyz(), symm)

  indices_dict = get_indices_convert_dict_from_array(phenix_coords_uc, parm.coordinates)
  p2a_indices = indices_dict['p2a']
  a2p_indices = indices_dict['a2p']

  a0 = reorder_coords_phenix_to_amber(phenix_coords_uc, p2a_indices)

  # make sure the original pdb and rst7 have the same coordinates
  aa_eq(pdb_input.atoms().extract_xyz(), parm.coordinates, decimal=2)

  # make sure reordered coordinates is equal to original coordinates
  aa_eq(a0, parm.coordinates, decimal=2)

  # raise KeyError if provide wrong coordinates
  def make_sure_key_error_being_raised():
    phenix_coords = flex.vec3_double([[0., 2., 3.2],  # should be [0., 2., 3.]
                                      [4., 5., 6.],
                                      [10., 12., 13.],
                                      [7., 8., 9.]])
    amber_coords = np.array([[10., 12., 13.],
                             [7., 8., 9],
                             [4., 5., 6.],
                             [0., 2., 3.]])
    get_indices_convert_dict_from_array(phenix_coords, amber_coords)

  with pytest.raises(KeyError):
      make_sure_key_error_being_raised()

def test_sander_struct():
  fn = get_fn('4lzt_pawel/4lztabH.pdb')
  tn = get_fn('4lzt_pawel/4lztab.parm7')
  rst7 = get_fn('4lzt_pawel/4lztabH.rst7')
  s_struct = SanderStruct(tn, rst7)
  parm = pmd.load_file(tn, rst7)
  aa_eq(parm.coordinates, s_struct.parm.coordinates)
