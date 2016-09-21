#!/usr/bin/env phenix.python
'''Testing for both LES and non-LES
'''
import os
import subprocess
import pytest
import libtbx.load_env
from amber_adaptbx.tests.utils import (tempfolder,
        get_prmtop_and_rst7_and_pdb_filenames_from_pdb,
)
from amber_adaptbx.tests.config import PDB_COLLECTION

@pytest.mark.medium
@pytest.mark.parametrize('pdb_file', PDB_COLLECTION)
@pytest.mark.parametrize('LES', [True, False])
@pytest.mark.parametrize('use_amber', [True, False])
def test_geometry_minimization_command_line(pdb_file, LES, use_amber):
  """ Test two commands with all combinations (>= 2 pdb files, LES=True/False, use_amber=True/False)
  """
  with tempfolder():
    prmtop_file, rst7_file, new_pdb = get_prmtop_and_rst7_and_pdb_filenames_from_pdb(pdb_file, LES=LES)
    command_build = [
            'phenix.AmberPrep',
            pdb_file,
            'LES={}'.format(LES)
    ]
    output = subprocess.check_output(command_build)
    command_minimization = [
            'phenix.geometry_minimization',
            new_pdb,
            'use_amber={}'.format(use_amber),
            'amber.topology_file_name={}'.format(prmtop_file),
            'amber.coordinate_file_name={}'.format(rst7_file),
            'max_iterations=1'
    ]
    print('geometry_minimization command: ', ' '.join(command_minimization))
    output = subprocess.check_output(command_minimization)
    assert 'Write PDB file' in output
    assert '_minimized.pdb' in output
