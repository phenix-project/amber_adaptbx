import os
import subprocess
import pytest
import libtbx.load_env
from numpy.testing import assert_almost_equal as aa_eq
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

  expected_r = {
          # 'Start R-work': 0.0154, 'Start R-free': 0.0154,
          'Final R-work': 0.0066, 'Final R-free': 0.0065
  }

  with tempfolder():
    output = subprocess.check_output(command_refine)
    for line in output.split('\n'):
        if 'Final R-work' in line:
            break
    final_r_work = float(line.split()[3].strip(','))
    final_r_free = float(line.split()[6])
    aa_eq([final_r_work,], [expected_r['Final R-work'],], decimal=3)
    aa_eq([final_r_free,], [expected_r['Final R-free'],], decimal=3)
