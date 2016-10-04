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

@pytest.mark.parametrize('pdb_file, expected_energy', [
    (get_fn('2igd/2igd_simplified.pdb'), -1664.0907),
])
def test_build_from_pdb_that_does_not_have_remark_290(pdb_file, expected_energy):
    command_build = [
            'phenix.AmberPrep',
            pdb_file
    ]
    with tempfolder():
        subprocess.check_call(command_build)
        prmtop_file, rst7_file, _ = get_prmtop_and_rst7_and_pdb_filenames_from_pdb(pdb_file, LES=False)
        box = pmd.load_file(rst7_file).box
        with sander.setup(prmtop_file, rst7_file, box=box, mm_options=sander.pme_input()):
            ene, _ = sander.energy_forces()
        aa_eq([ene.tot], [expected_energy], decimal=4)
