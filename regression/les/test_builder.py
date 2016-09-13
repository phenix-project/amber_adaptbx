#!/usr/bin/env phenix.python
import os
import subprocess
import pytest
import libtbx.load_env
from amber_adaptbx.regression.utils import (tempfolder, get_fn,
        assert_energy_and_forces,
        run_sander_minimization,
)

# TODO: add more pdb code
asu_pdb = get_fn('2igd.pdb')
saved_rst7_file = get_fn('2igdab.LES.rst7')
saved_prmtop_file = get_fn('2igdab.LES.prmtop')
saved_pdb_file = get_fn('2igdab_4phenix.LES.pdb')

pdb_code = '2igd'
prmtop_file = '4amber_{code}.LES.prmtop'.format(code=pdb_code)
rst7_file = '4amber_{code}.LES.rst7'.format(code=pdb_code)

def test_run_sander_LES_min():
  """ ensure we can run minimization """
  with tempfolder():
    output = run_sander_minimization(prmtop_file=saved_prmtop_file,
              rst7_file=saved_rst7_file)
    assert 'FINAL RESULTS' in output
    
def test_command_line_build():
  command = 'phenix.AmberPrep {} LES=True'.format(asu_pdb)
  with tempfolder():
    subprocess.check_output(command.split())
    # ensure not error
    output = run_sander_minimization(prmtop_file=prmtop_file, rst7_file=rst7_file)
    assert 'FINAL RESULTS' in output
    assert_energy_and_forces(prmtop_file, rst7_file,
            saved_prmtop_file, saved_rst7_file)

@pytest.mark.minimization
def test_geometry_minimization():
  """ Test two commands

  Examples
  --------
  phenix.AmberPrep 2igd.pdb
  phenix.geometry_minimization 4phenix_2igd.LES.pdb use_amber=True \
      amber.topology_file_name=4amber_2igd.LES.prmtop \
      amber.coordinate_file_name=4amber_2igd.LES.rst7 \
      max_iterations=2
  """
  command = 'phenix.AmberPrep {} LES=True | grep phenix.geometry_minimization'.format(asu_pdb)
  with tempfolder():
    output = subprocess.check_output(command, shell=True)
    # e.g output = phenix.geometry_minimization 4phenix_2igd.LES.pdb use_amber=True \ 
    # amber.topology_file_name=4amber_2igd.LES.prmtop amber.coordinate_file_name=4amber_2igd.LES.rst7

    # run very short minimization
    minization_command = output.strip() + ' max_iterations=2'
    minimization_output = subprocess.check_output(minization_command.split())

@pytest.mark.slowtest
@pytest.mark.no_assertion
def test_command_line_minimization_amber_h():
  """ ensure there is no error, there is no assertion """
  command = 'phenix.AmberPrep {} LES=True minimise=amber_h'.format(asu_pdb)
  with tempfolder():
    subprocess.check_output(command.split())

@pytest.mark.slowtest
@pytest.mark.no_assertion
def test_command_line_minimization_phenix_all():
  """ ensure there is no error, there is no assertion """
  command = 'phenix.AmberPrep {} LES=True minimise=phenix_all'.format(asu_pdb)
  with tempfolder():
    subprocess.check_output(command.split())

@pytest.mark.slowtest
@pytest.mark.refinement
@pytest.mark.pdb_2igd
@pytest.mark.no_assertion
def test_refine():
  """ just ensure no error. There is no assertion here """
  command_list = [
          'phenix.refine',
          saved_pdb_file,
          get_fn('2igd.mtz'),
          'topology_file_name={}'.format(saved_prmtop_file),
          'amber.coordinate_file_name={}'.format(saved_rst7_file),
          'use_amber=True',
          'strategy=individual_sites+individual_adp+occupancies',
          'write_geo=False',
          'refinement.main.number_of_macro_cycles=1',
  ]
  with tempfolder():
    output = subprocess.check_output(command_list)
