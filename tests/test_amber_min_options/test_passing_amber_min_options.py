import os
import subprocess
import time
import pytest
from amber_adaptbx.tests.utils import (
    get_fn,
    tempfolder,
    assert_file_has_line,
    get_prmtop_and_rst7_and_pdb_filenames_from_pdb,
    get_minimized_pdb_filename,
    get_minimized_rst7_filename, )
from amber_adaptbx.tests.config import PDB_COLLECTION


@pytest.mark.parametrize('pdb_file', PDB_COLLECTION)
@pytest.mark.parametrize('minimization_type', ['amber_h', 'amber_all'])
@pytest.mark.parametrize('minimization_options', [
    'maxcyc=10',
    'maxcyc=10, restraint_wt=10.',
])
def test_passing_minimization_options(pdb_file, minimization_type,
                                      minimization_options):
    command_min = [
        'phenix.AmberPrep', pdb_file, 'minimise={}'.format(minimization_type),
        'clean=False', 'minimization_options="{}"'.format(minimization_options)
    ]

    print('\n-->' + ' '.join(command_min))
    with tempfolder():
        base = os.path.basename(pdb_file).split('.')[0]
        mdout_file_name = base + '.min.out'
        subprocess.check_call(command_min)
        line = 'maxcyc  =      10, ncyc    =     200, ntmin   =       1'
        assert_file_has_line(mdout_file_name, line)
        if minimization_options == 'maxcyc=10, restraint_wt=10.':
            assert_file_has_line(mdout_file_name, 'restraint_wt =  10.00000')


@pytest.mark.parametrize('pdb_file', PDB_COLLECTION[:1])
@pytest.mark.parametrize('minimization_type', ['amber_h', 'amber_all'])
def test_default_min(pdb_file, minimization_type):
    command_min = [
        'phenix.AmberPrep', pdb_file, 'clean=False',
        'minimise={}'.format(minimization_type)
    ]

    with tempfolder():
        base = os.path.basename(pdb_file).split('.')[0]
        mdout_file_name = base + '.min.out'
        print('--> mdout_file_name = {}'.format(mdout_file_name))
        process = subprocess.Popen(command_min)
        # we don't need to have minimization finished
        # only care about the input part of mdout
        time.sleep(15)
        process.terminate()
        if minimization_type == 'amber_h':
            maxcyc = 1000
        else:
            maxcyc = 500
        if minimization_type == 'amber_h':
            line = 'maxcyc  =    1000, ncyc    =     200, ntmin   =       1'
        elif minimization_type == 'amber_all':
            line = 'maxcyc  =      50, ncyc    =     200, ntmin   =       1'
        assert_file_has_line(mdout_file_name, line)


@pytest.mark.parametrize('pdb_file', [get_fn('2igd/2igd.pdb')])
@pytest.mark.parametrize('minimization_type', ['amber_all', 'amber_h'])
@pytest.mark.parametrize('LES', [True, False])
@pytest.mark.parametrize('mtz_file', [get_fn('2igd/2igd.mtz')])
@pytest.mark.slow
def test_run_refinement_after_minimization_that_used_sander(
        pdb_file, minimization_type, LES, mtz_file):
    # TODO: assertion?
    if minimization_type in ['amber_h', 'amber_all']:
        minimization_options = "maxcyc=100"
    else:
        raise ValueError(
            'wrong minimization_type: {}'.format(minimization_type))
    prmtop_file, _, _ = get_prmtop_and_rst7_and_pdb_filenames_from_pdb(
        pdb_file, LES=LES)
    minimized_rst7_file = get_minimized_rst7_filename(
        pdb_file, minimization_type=minimization_type, LES=LES)
    minimized_pdb_file = get_minimized_pdb_filename(
        pdb_file, minimization_type=minimization_type, LES=LES)
    command_min = [
        'phenix.AmberPrep', pdb_file, 'minimise={}'.format(minimization_type),
        'LES={}'.format(LES), 'clean=False',
        'minimization_options="{}"'.format(minimization_options)
    ]
    command_refinement = [
        'phenix.refine',
        minimized_pdb_file,
        mtz_file,
        'amber.topology_file_name = {}'.format(prmtop_file),
        'amber.coordinate_file_name = {}'.format(minimized_rst7_file),
        'refinement.main.number_of_macro_cycles=2',
    ]
    print('\n-->' + ' '.join(command_min))
    print('\n-->' + ' '.join(command_refinement))
    with tempfolder():
        output_min = subprocess.check_output(command_min)
        output_refinement = subprocess.check_output(command_refinement)
