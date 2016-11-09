import os
import subprocess
import pytest
import libtbx.load_env
from amber_adaptbx.tests.utils import tempfolder, get_fn
import parmed as pmd

def test_keeping_original_resnum():
  code = '1gdu'
  original_pdb = get_fn('{}/{}.pdb'.format(code, code))
  command_build = [
    'phenix.AmberPrep',
    original_pdb
  ]

  # LES=False
  with tempfolder():
    subprocess.check_call(command_build)
    parm = pmd.load_file('4phenix_{}.pdb'.format(code))
    assert parm.residues[0].number == 16
    assert parm.residues[-1].number == 2201

  # LES=True
  with tempfolder():
    subprocess.check_call(command_build + ['LES=True'])
    parm = pmd.load_file('4phenix_{}.LES.pdb'.format(code))
    assert parm.residues[0].number == 16
    assert parm.residues[-1].number == 2201
