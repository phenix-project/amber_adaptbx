import os
import subprocess
from glob import glob
import pytest
from numpy.testing import assert_almost_equal as aa_eq
import libtbx.load_env
from amber_adaptbx.tests.utils import get_fn, tempfolder
import parmed as pmd

def equal_bfactor_and_occupancy(parm1, parm2):
    for atom1, atom2 in zip(parm1.atoms, parm2.atoms):
        assert atom1.bfactor == atom2.bfactor
        assert atom1.occupancy == atom2.occupancy

@pytest.mark.parametrize('code', ['1gdu_tiny'])
def test_bfactor_occupancy_by_comparing_pdb_files(code):
 pdb_fn = get_fn('fake/{}.pdb'.format(code))

 # 3-name rules, there is no ".LES.' in final pdb anymore.
 final_pdb_fn = '4phenix_{}.pdb'.format(code)
 les_pdb = '4phenix_{}.LES.pdb'.format(code)

 parm_expected_pdb_non_les = pmd.load_file(get_fn('fake/{}'.format(final_pdb_fn)))
 parm_expected_pdb_les = pmd.load_file(get_fn('fake/{}'.format(les_pdb)))

 command_build = [
         'phenix.AmberPrep',
         pdb_fn,
         'use_amber_unitcell=False',
 ]

 # 'LES=False'
 with tempfolder():
   subprocess.check_call(command_build + ['LES=False',])
   parm_non_les = pmd.load_file(final_pdb_fn)
   subprocess.check_call(command_build + ['LES=True',])
   parm_les = pmd.load_file(final_pdb_fn)

 equal_bfactor_and_occupancy(parm_non_les, parm_expected_pdb_non_les)
 equal_bfactor_and_occupancy(parm_les, parm_expected_pdb_les)
