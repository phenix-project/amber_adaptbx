import os
import subprocess
import pytest
import libtbx.load_env
from amber_adaptbx.tests.utils import (tempfolder, get_fn,
        assert_energy_and_forces,
        run_sander_minimization,
        get_prmtop_and_rst7_and_pdb_filenames_from_pdb,
)

from amber_adaptbx.tests.config import (PDB_COLLECTION, PDB_MTZ_COLLECTION,
        MTZ_COLLECTION,
        saved_2igd_prmtop_file,
        saved_2igd_rst7_file,
        saved_2igd_pdb_file,
        saved_2igd_mtz_file
)


@pytest.mark.slow
@pytest.mark.parametrize('prmtop_file, rst7_file, pdb_file, mtz_file, use_amber', [
   (saved_2igd_prmtop_file, saved_2igd_rst7_file, saved_2igd_pdb_file, saved_2igd_mtz_file, True),
   (saved_2igd_prmtop_file, saved_2igd_rst7_file, saved_2igd_pdb_file, saved_2igd_mtz_file, False),
])
def test_refinement_from_saved_2igd_files(prmtop_file, rst7_file, pdb_file, mtz_file, use_amber):
  """ just ensure no error. There is no assertion here """
  command_list = [
          'phenix.refine',
          pdb_file,
          mtz_file,
          'topology_file_name={}'.format(prmtop_file),
          'amber.coordinate_file_name={}'.format(rst7_file),
          'use_amber={}'.format(use_amber),
          'strategy=individual_sites+individual_adp+occupancies',
          'write_geo=False',
          'refinement.main.number_of_macro_cycles=1',
  ]
  with tempfolder():
    subprocess.call(command_list)

@pytest.mark.slow
@pytest.mark.parametrize('pdb_file, mtz_file', PDB_MTZ_COLLECTION)
def test_refinement_after_running_AmberPrep_use_amber_True(pdb_file, mtz_file):
  """ just ensure no error. There is no assertion here """
  command_build = [
          'phenix.AmberPrep',
          pdb_file,
          'LES=True'
  ]

  prmtop_file, rst7_file, new_pdb_file = get_prmtop_and_rst7_and_pdb_filenames_from_pdb(pdb_file)
  use_amber = True

  command_refine = [
          'phenix.refine',
          new_pdb_file,
          mtz_file,
          'topology_file_name={}'.format(prmtop_file),
          'amber.coordinate_file_name={}'.format(rst7_file),
          'use_amber={}'.format(use_amber),
          'strategy=individual_sites+individual_adp+occupancies',
          'write_geo=False',
          'refinement.main.number_of_macro_cycles=1',
          'refinement.input.xray_data.r_free_flags.generate=True',
  ]
  with tempfolder():
    subprocess.check_call(command_build)
    subprocess.check_call(command_refine)

@pytest.mark.slow
@pytest.mark.parametrize('pdb_file, mtz_file', PDB_MTZ_COLLECTION)
def test_refinement_after_running_AmberPrep_use_amber_True_from_DAC_script(pdb_file, mtz_file):
  """ just ensure no error. There is no assertion here """
  command_build = [
          'phenix.AmberPrep',
          pdb_file,
          'LES=True'
  ]
  prmtop_file, rst7_file, new_pdb_file = get_prmtop_and_rst7_and_pdb_filenames_from_pdb(pdb_file, LES=True)
  use_amber = True
  command_refine = [
          'phenix.refine',
           new_pdb_file,
           mtz_file,
          'topology_file_name={}'.format(prmtop_file),
          'amber.coordinate_file_name={}'.format(rst7_file),
          'use_amber={}'.format(use_amber),
          'refinement.target_weights.optimize_xyz_weight=False',
          'strategy=individual_sites+individual_adp',
          'prefix=amber',
          'serial=1',
          'write_geo=False',
          'refinement.main.number_of_macro_cycles=10',
  ]
  # below command is not from Dave's script
  # Why adding? because got error: No array of R-free flags found.
  # phenix suggested: refinement.input.xray_data.r_free_flags.generate=True
  # TODO
  command_refine.append('refinement.input.xray_data.r_free_flags.generate=True')
  with tempfolder():
    subprocess.check_call(command_build)
    subprocess.check_call(command_refine)
