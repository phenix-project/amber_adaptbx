import os
import subprocess
import pytest
import libtbx.load_env
from amber_adaptbx.tests.utils import tempfolder, get_fn
import parmed as pmd

def test_ensure_amber_and_phenix_have_the_same_atom_order_in_pdb_output():
  original_pdb = get_fn('2igd/2igd.pdb')
  command_build = [
    'phenix.AmberPrep',
    original_pdb
  ]

  with tempfolder():
    fn = '4phenix_2igd.pdb'
    subprocess.check_call(command_build + ['LES=False'])
    non_les_parm = pmd.load_file(fn)
  with tempfolder():
    subprocess.check_call(command_build + ['LES=True'])
    fn = '4phenix_2igd.LES.pdb'
    les_parm = pmd.load_file(fn)

  for atom_non_les_parm, atom_les_parm in zip(non_les_parm.residues[0], les_parm.residues[0]):
    assert atom_non_les_parm.name == atom_les_parm.name
