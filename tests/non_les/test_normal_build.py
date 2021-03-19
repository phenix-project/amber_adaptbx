#!/usr/bin/env phenix.python
from __future__ import print_function
import os
import subprocess
import pytest
import libtbx.load_env
from amber_adaptbx.tests.utils import (
    tempfolder,
    get_fn,
    assert_energy_and_forces,
    run_sander_minimization,
    get_prmtop_and_rst7_and_pdb_filenames_from_pdb, )


@pytest.mark.slow
@pytest.mark.minimization
@pytest.mark.parametrize('pdb_file', [
    get_fn('2igd/2igd.pdb'),
    get_fn('4lzt/4lzt_no_BHOH.pdb'),
    get_fn('4lzt/4lzt.pdb'),
])
def test_geometry_minimization(pdb_file):
    """
  """
    command_build = 'phenix.AmberPrep {}'.format(pdb_file)
    with tempfolder():
        prmtop_file, rst7_file, new_pdb = get_prmtop_and_rst7_and_pdb_filenames_from_pdb(
            pdb_file)
        output = subprocess.check_output(command_build, shell=True)
        command_list = [
            'phenix.geometry_minimization',
            new_pdb,
            'use_amber=True',
            'amber.topology_file_name={}'.format(prmtop_file),
            'amber.coordinate_file_name={}'.format(rst7_file),
        ]
        print(('geometry_minimization command: ', ' '.join(command_list)))
        output = subprocess.check_output(command_list)
        assert 'Write PDB file' in output
        assert '_minimized.pdb' in output


@pytest.mark.slow
@pytest.mark.no_assertion
@pytest.mark.parametrize('pdb_file', [
    get_fn('2igd/2igd.pdb'),
    get_fn('4lzt/4lzt_no_BHOH.pdb'),
    get_fn('4lzt/4lzt.pdb'),
])
def test_geometry_minimization_command_line_with_amber_h_option(pdb_file):
    """ ensure there is no error, there is no assertion """
    command = 'phenix.AmberPrep {} minimise=amber_h'.format(pdb_file)
    with tempfolder():
        subprocess.check_output(command.split())


@pytest.mark.slow
@pytest.mark.no_assertion
@pytest.mark.phenix
@pytest.mark.parametrize('pdb_file, mtz_file', [
    (get_fn('4lzt/4lzt.pdb'), get_fn('4lzt/4lzt.mtz')),
    (get_fn('4lzt/4lzt_no_BHOH.pdb'), get_fn('4lzt/4lzt.mtz')),
    (get_fn('2igd/2igd.pdb'), get_fn('2igd/2igd.mtz')),
])
def test_refinement_phenix(pdb_file, mtz_file):
    """ ensure there is no error, there is no assertion """
    print(('pdb_file', pdb_file))
    prmtop_file, rst7_file, new_pdb_file = get_prmtop_and_rst7_and_pdb_filenames_from_pdb(
        pdb_file)

    command_build = [
        'phenix.AmberPrep', '{pdb_file}'.format(pdb_file=pdb_file)
    ]
    command_refine = [
        'phenix.refine',
        '{pdb_file}'.format(pdb_file=new_pdb_file),
        '{mtz_file}'.format(mtz_file=mtz_file),
        'amber.topology_file_name={}'.format(prmtop_file),
        'amber.coordinate_file_name={}'.format(rst7_file),
        'refinement.main.number_of_macro_cycles=2',
    ]
    if '4lzt' in pdb_file:
        command_refine.append(
            'refinement.input.xray_data.r_free_flags.generate=True')
    with tempfolder():
        subprocess.check_output(command_build)
        print(' '.join(command_refine))
        # subprocess.check_output(command_refine)
        subprocess.call(command_refine)
