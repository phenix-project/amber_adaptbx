import os
import subprocess
import pytest
import libtbx.load_env
from amber_adaptbx.tests.utils import tempfolder, get_fn

@pytest.mark.parametrize('use_amber', [True])
def test_non_LES_refinement_vAla3(use_amber):
  command_refine = [
    'phenix.refine',
    get_fn('vAla3/vAla3.pdb'),
    get_fn('vAla3/vAla3.cif'),
    get_fn('vAla3/vAla3.mtz'),
    'topology_file_name={}'.format(get_fn('vAla3/vAla3.prmtop')),
    'coordinate_file_name={}'.format(get_fn('vAla3/vAla3.rst7')),
    'use_amber={}'.format(use_amber),
    'wxc_scale=0.025',
    '--overwrite',
    'refinement.main.number_of_macro_cycles=1'
  ]

  expected_lines = """
Start R-work = 0.0150, R-free = 0.0149
Final R-work = 0.0061, R-free = 0.0060
"""

  with tempfolder():
    output = subprocess.check_output(command_refine)
    assert expected_lines in output 
    print(output[:-100])
