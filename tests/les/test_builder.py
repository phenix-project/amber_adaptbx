#!/usr/bin/env phenix.python
import os
import subprocess
import pytest
import libtbx.load_env
from amber_adaptbx.tests.utils import (tempfolder, get_fn,
        assert_energy_and_forces,
        run_sander_minimization,
        get_prmtop_and_rst7_and_pdb_filenames_from_pdb,
        get_minimized_pdb_filename,
        get_minimized_rst7_filename,
)
from amber_adaptbx.tests.config import (PDB_COLLECTION, saved_2igd_prmtop_file,
        saved_2igd_rst7_file,
)

@pytest.mark.parametrize('prmtop_file, rst7_file', [(saved_2igd_prmtop_file, saved_2igd_rst7_file)])
@pytest.mark.saved
def test_run_sander_LES_minimization_from_saved_files(prmtop_file, rst7_file):
  """ ensure we can run minimization """
  with tempfolder():
    output = run_sander_minimization(prmtop_file=prmtop_file,
              rst7_file=rst7_file,
              maxcyc=10)
    assert 'FINAL RESULTS' in output

@pytest.mark.parametrize('pdb_file', PDB_COLLECTION)
def test_run_sander_LES_minimization_from_LES_build(pdb_file):
  """ ensure we can run minimization, maxcyc=10
  Tests for 2igd.pdb, 4lzt_no_BHOH.pdb
  """
  with tempfolder():
    command = 'phenix.AmberPrep {pdb_file} LES=True'.format(pdb_file=pdb_file)
    # run building LES prmtop
    subprocess.check_call(command.split())
    prmtop_file, rst7_file,_ = get_prmtop_and_rst7_and_pdb_filenames_from_pdb(pdb_file, LES=True)
    output = run_sander_minimization(prmtop_file=prmtop_file,
              rst7_file=rst7_file,
              maxcyc=10)
    assert 'FINAL RESULTS' in output, 'minimization must be finished'

@pytest.mark.parametrize('pdb_file', [
    get_fn('2g3i/2g3i.pdb'),
])
def test_minus_0_coordinates_use_amber_all_for_minimise(pdb_file):
  """ ensure there is no error if having -0.0 coordinates """
  command_build = [
          'phenix.AmberPrep',
          pdb_file,
          'LES=True',
          'minimise=amber_all',
          "minimization_options='maxcyc=10'"
  ]
  amber_minimimized_pdb = get_minimized_pdb_filename(pdb_file, LES=True, minimization_type='amber_all')
  prmtop, _, _ = get_prmtop_and_rst7_and_pdb_filenames_from_pdb(pdb_file, LES=True)
  minimized_rst7 = get_minimized_rst7_filename(pdb_file, minimization_type='amber_all', LES=True)

  command_geom = [
    'phenix.geometry_minimization',
    amber_minimimized_pdb,
    'amber.topology_file_name',
    prmtop,
    'amber.coordinate_file_name',
    minimized_rst7

  ]
  with tempfolder():
    subprocess.check_call(command_build)
    # use geometry_minimization to trigger computing reorder map
    subprocess.check_call(command_geom)
    
@pytest.mark.parametrize('pdb_file', PDB_COLLECTION)
def test_command_line_build(pdb_file):
  command = 'phenix.AmberPrep {} LES=True'.format(pdb_file)
  prmtop_file, rst7_file, _ = get_prmtop_and_rst7_and_pdb_filenames_from_pdb(pdb_file)
  with tempfolder():
    subprocess.check_output(command.split())
    # ensure not error
    output = run_sander_minimization(prmtop_file=prmtop_file, rst7_file=rst7_file)
    assert 'FINAL RESULTS' in output

@pytest.mark.parametrize('pdb_file', PDB_COLLECTION)
def test_geometry_minimization_from_AmberPrep_with_amber_h_option(pdb_file):
  """ ensure there is no error, there is no assertion """
  command = "phenix.AmberPrep {} LES=True minimise=amber_h minimization_options='maxcyc=2'".format(pdb_file)
  with tempfolder():
    subprocess.check_output(command.split())

@pytest.mark.parametrize('pdb_file', [
    get_fn('4lzt/4lzt_no_BHOH.pdb'),
    get_fn('2igd/2igd.pdb')
])
def test_geometry_minimization_from_AmberPrep_with_amber_all_option(pdb_file):
  """ ensure there is no error, there is no assertion """
  command = "phenix.AmberPrep {} LES=True minimise=amber_all minimization_options='maxcyc=2'".format(pdb_file)
  with tempfolder():
    prmtop_file, rst7_file, new_pdb_file = get_prmtop_and_rst7_and_pdb_filenames_from_pdb(pdb_file)
    subprocess.check_output(command.split())

@pytest.mark.parametrize('pdb_file', [
    get_fn('2igd/2igd.pdb'),
])
def test_geometry_minimization_from_AmberPrep_with_amber_all_option_with_assertion_and_larger_maxcyc(pdb_file):
  """ ensure there is no error, there is no assertion """
  command = [
          "phenix.AmberPrep",
          pdb_file,
          "LES=True", "minimise=amber_all",
          "minimization_options='maxcyc=2'",
          "clean=off",
  ]
  with tempfolder():
    prmtop_file, rst7_file, new_pdb_file = get_prmtop_and_rst7_and_pdb_filenames_from_pdb(pdb_file)
    subprocess.check_output(command)
    assert(os.path.exists(get_minimized_pdb_filename(pdb_file, minimization_type='amber_all', LES=True)))
