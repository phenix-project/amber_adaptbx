#!/usr/bin/env phenix.python
import subprocess
import pytest
from numpy.testing import assert_almost_equal as aa_eq
from amber_adaptbx.tests.utils import (tempfolder, get_fn,
        get_energy_and_forces,
        get_prmtop_and_rst7_and_pdb_filenames_from_pdb,
)
import sander
import parmed as pmd

@pytest.mark.parametrize('pdb_file', [
    get_fn('2igd/2igd_simplified.pdb')
])
def test_build_from_pdb_that_does_not_have_remark_290(pdb_file):
    command_build = [
            'phenix.AmberPrep',
            pdb_file
    ]

    with tempfolder():
        # just make sure no error
        subprocess.check_call(command_build)
