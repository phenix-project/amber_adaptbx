import os
import subprocess
import pytest
from numpy.testing import assert_almost_equal as aa_eq
import libtbx.load_env
from amber_adaptbx.tests.utils import get_fn, tempfolder
from glob import glob

def equal_files(fn1, fn2):
  # create equal_files here so pytest will raise difference if getting error
  with open(fn1) as fh1:
    with open(fn2) as fh2:
      for line0, line1 in zip(fh1.readlines(), fh2.readlines()):
        assert line0 == line1

@pytest.mark.parametrize('code', ['1gdu_tiny'])
def test_bfactor_occupancy_by_comparing_pdb_files(code):
 pdb_fn = get_fn('fake/{}.pdb'.format(code))

 non_les_pdb = '4phenix_{}.pdb'.format(code)
 les_pdb = '4phenix_{}.LES.pdb'.format(code)
 expected_pdb_non_les = get_fn('fake/{}'.format(non_les_pdb))
 expected_pdb_les = get_fn('fake/{}'.format(les_pdb))

 command_build = [
         'phenix.AmberPrep',
         pdb_fn,
         'LES=True'
 ]

 with tempfolder():
   subprocess.check_call(command_build)
   equal_files(non_les_pdb, expected_pdb_non_les)
   equal_files(les_pdb, expected_pdb_les)
