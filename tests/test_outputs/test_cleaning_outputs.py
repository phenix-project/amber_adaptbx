import os
import subprocess
import pytest
from numpy.testing import assert_almost_equal as aa_eq
import libtbx.load_env
from amber_adaptbx.tests.utils import get_fn, tempfolder
from glob import glob


@pytest.mark.parametrize('code', ['2igd'])
def test_cleaning_unused_files(code):
    pdb_fn = get_fn('{code}/{code}.pdb'.format(code=code))

    command_non_les = ['phenix.AmberPrep', pdb_fn, 'LES=False']
    expected_files_non_les = set([
        '4phenix_{}.pdb'.format(code),
        '4amber_{}.prmtop'.format(code),
        '4amber_{}.rst7'.format(code),
    ])

    command_les = ['phenix.AmberPrep', pdb_fn, 'LES=True']

    expected_files_les = set()
    expected_files_les.update(expected_files_non_les)
    expected_files_les.update(
        set([
            '4phenix_{}.pdb'.format(code),
            '4amber_{}.prmtop'.format(code),
            '4amber_{}.rst7'.format(code),
            '4amber_{}.LES.pdb'.format(code),
            'addles.in',
        ]))

    with tempfolder():
        subprocess.check_call(command_les)
        file_set = set(glob('*'))
        assert sorted(expected_files_les) == sorted(file_set)

    with tempfolder():
        subprocess.check_call(command_non_les)
        file_set = set(glob('*'))
        assert sorted(expected_files_non_les) == sorted(file_set)
