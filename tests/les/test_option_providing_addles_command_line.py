import os
import subprocess
import pytest
import libtbx.load_env
from amber_adaptbx.tests.utils import tempfolder, get_fn
import parmed as pmd

def test_adding_addles_option_to_command_line():
  addles_input = get_fn('2igd/addles_test/addles.in') 
  addles_input_wrong_header  = get_fn('2igd/addles_test/addles.wrong_header.in') 
  pdb_fn = get_fn('2igd/2igd.pdb')
  prmtop = '4amber_2igd.prmtop'

  command_build_addles_input_cmd = [
          'phenix.AmberPrep',
          pdb_fn,
          'LES=True',
          'addles_input={}'.format(addles_input)
  ]
  with tempfolder():
    subprocess.check_call(command_build_addles_input_cmd)
    parm = pmd.load_file(prmtop)
    assert len(parm.atoms) == 5148

  command_build_normal = [
          'phenix.AmberPrep',
          pdb_fn,
          'LES=True',
  ]
  with tempfolder():
    subprocess.check_call(command_build_normal)
    parm = pmd.load_file(prmtop)
    assert len(parm.atoms) == 5240

  command_build_wrong_header = [
          'phenix.AmberPrep',
          pdb_fn,
          'LES=True',
          'addles_input={}'.format(addles_input_wrong_header)
  ]
  with tempfolder():
    with pytest.raises(subprocess.CalledProcessError):
      subprocess.check_call(command_build_wrong_header)

  command_build_non_existing_file = [
          'phenix.AmberPrep',
          pdb_fn,
          'LES=True',
          'addles_input={}'.format('some_weird_random_name_you_are_not_lucky.in')
  ]
  with tempfolder():
    with pytest.raises(subprocess.CalledProcessError):
      subprocess.check_call(command_build_non_existing_file)
