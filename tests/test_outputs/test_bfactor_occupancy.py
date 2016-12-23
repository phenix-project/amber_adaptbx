import os
import subprocess
from glob import glob
import pytest
from numpy.testing import assert_almost_equal as aa_eq
import libtbx.load_env
from amber_adaptbx.tests.utils import get_fn, tempfolder
import parmed as pmd

def equal_bfactor_and_occupancy(fn1, fn2):
    parm1 = pmd.load_file(fn1)
    parm2 = pmd.load_file(fn2)

    for atom1, atom2 in zip(parm1.atoms, parm2.atoms):
        assert atom1.bfactor == atom2.bfactor
        assert atom1.occupancy == atom2.occupancy

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
   equal_bfactor_and_occupancy(non_les_pdb, expected_pdb_non_les)
   equal_bfactor_and_occupancy(les_pdb, expected_pdb_les)
