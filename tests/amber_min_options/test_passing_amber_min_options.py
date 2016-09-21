import os
import subprocess
import pytest
from amber_adaptbx.tests.utils import get_fn, tempfolder, assert_file_has_line
from amber_adaptbx.tests.config import PDB_COLLECTION

@pytest.mark.medium
@pytest.mark.parametrize('pdb_file', PDB_COLLECTION)
@pytest.mark.parametrize('task', ['amber_h', 'amber_all'])
@pytest.mark.parametrize('amber_min_options', [
    'maxcyc=10',
    'maxcyc=10, restraint_wt=10.',
])
def test_passing_amber_min_options(pdb_file, task, amber_min_options):
  command_min = [
          'phenix.AmberPrep',
          pdb_file,
          'minimise={}'.format(task),
          'amber_min_options="{}"'.format(amber_min_options)
  ]

  print('\n-->' + ' '.join(command_min))
  with tempfolder():
    base = os.path.basename(pdb_file).split('.')[0]
    mdout_file_name = base + '_' + task + '.out'
    subprocess.check_call(command_min)
    line = 'maxcyc  =      10, ncyc    =     200, ntmin   =       1'
    assert_file_has_line(mdout_file_name, line)
    if amber_min_options == 'maxcyc=10, restraint_wt=10.':
      assert_file_has_line(mdout_file_name, 'restraint_wt =  10.00000')

@pytest.mark.slow
@pytest.mark.parametrize('pdb_file', PDB_COLLECTION[:1])
@pytest.mark.parametrize('task', ['amber_h', 'amber_all'])
def test_default_min(pdb_file, task):
  command_min = [
          'phenix.AmberPrep',
          pdb_file,
          'minimise={}'.format(task)
  ]

  with tempfolder():
    base = os.path.basename(pdb_file).split('.')[0]
    mdout_file_name = base + '_' + task + '.out'
    print('--> mdout_file_name = {}'.format(mdout_file_name))
    subprocess.check_call(command_min)
    if task == 'amber_h':
      maxcyc = 1000
    else:
      maxcyc = 500
    if task == 'amber_h':
      line = 'maxcyc  =    1000, ncyc    =     200, ntmin   =       1'
    elif task == 'amber_all':
      line = 'maxcyc  =     500, ncyc    =     200, ntmin   =       1'
    assert_file_has_line(mdout_file_name, line)
