import os
import subprocess
import pytest
import libtbx.load_env
from amber_adaptbx.tests.utils import tempfolder, get_fn
import parmed as pmd

def test_ensure_amber_and_phenix_have_the_same_atom_order_in_pdb_output():
  command_refine_template = [
    'phenix.refine',
    '4phenix_2igd.pdb',
    get_fn('2igd/2igd.mtz'),
    'topology_file_name=4amber_2igd.prmtop',
    'coordinate_file_name=4amber_2igd.rst7',
    'wxc_scale=0.025',
    '--overwrite',
  ]

  command_build = [
    'phenix.AmberPrep',
    get_fn('2igd/2igd.pdb'),
  ]

  command_refine_use_amber = command_refine_template[:]
  command_refine_use_amber.append('use_amber=True')
  command_refine_use_phenix = command_refine_template[:]
  command_refine_use_phenix.append('use_amber=False')
  print(command_refine_use_amber)
  print(command_refine_use_phenix)

  pdb_fn = '4phenix_2igd_refine_001.pdb'
  with tempfolder():
    subprocess.check_call(command_build)
    subprocess.check_call(command_refine_use_amber)
    amber_parm = pmd.load_file(pdb_fn)
  with tempfolder():
    subprocess.check_call(command_build)
    subprocess.check_call(command_refine_use_phenix)
    phenix_parm = pmd.load_file(pdb_fn)
  for atom_amber, atom_phenix in zip(amber_parm.atoms, phenix_parm.atoms):
      assert atom_amber.name == atom_phenix.name
