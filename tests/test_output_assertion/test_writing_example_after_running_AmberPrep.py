import subprocess
import pytest
from amber_adaptbx.tests.utils import (tempfolder, get_fn,
        assert_energy_and_forces,
        run_sander_minimization,
        get_prmtop_and_rst7_and_pdb_filenames_from_pdb,
        get_minimized_pdb_filename,
)
from amber_adaptbx.tests.config import (PDB_COLLECTION, saved_2igd_prmtop_file,
        saved_2igd_rst7_file,
)

def test_writing_example_after_running_AmberPrep_non_LES_without_minimization():
  pdb_file = get_fn('2igd/2igd.pdb')
  LES = False

  expected_line = """
==================================================
Done.  Three new files have been made:
      4phenix_2igd.pdb
      4amber_2igd.prmtop
      4amber_2igd.rst7
==================================================
"""
  # always create new command_build
  command_build = [
          'phenix.AmberPrep',
          pdb_file,
          'LES={}'.format(LES),
  ]
  with tempfolder():
    output_build = subprocess.check_output(command_build)
    assert expected_line in output_build

def test_writing_example_after_running_AmberPrep_LES_without_minimization():
  pdb_file = get_fn('2igd/2igd.pdb')
  LES = True

  expected_line = """
==================================================
Done.  Three new files have been made:
      4phenix_2igd.LES.pdb
      4amber_2igd.LES.prmtop
      4amber_2igd.LES.rst7
==================================================
"""
  # always create new command_build
  command_build = [
          'phenix.AmberPrep',
          pdb_file,
          'LES={}'.format(LES),
  ]
  with tempfolder():
    output_build = subprocess.check_output(command_build)
    assert expected_line in output_build

@pytest.mark.parametrize('minimization_type', ['phenix_all', 'amber_h', 'amber_all'])
@pytest.mark.medium
def test_writing_example_after_running_AmberPrep_LES_with_minimization(minimization_type):
  pdb_file = get_fn('2igd/2igd.pdb')
  LES = True

  expected_line_dict = {
          'amber_h': """
==================================================
Done.  Three new files have been made:
      4phenix_2igd.LES.min.amber_h.pdb
      4amber_2igd.LES.prmtop
      4amber_2igd.LES.min.amber_h.rst7
==================================================
""",
          'amber_all': """
==================================================
Done.  Three new files have been made:
      4phenix_2igd.LES.min.amber_all.pdb
      4amber_2igd.LES.prmtop
      4amber_2igd.LES.min.amber_all.rst7
==================================================
""",
          'phenix_all': """
==================================================
Done.  Three new files have been made:
      4phenix_2igd.LES.min.phenix_all.pdb
      4amber_2igd.LES.prmtop
      4amber_2igd.LES.rst7
==================================================
"""
  }
  # always create new command_build
  command_build = [
          'phenix.AmberPrep',
          pdb_file,
          'LES={}'.format(LES),
  ]
  with tempfolder():
    command_build.append('minimise={}'.format(minimization_type))
    if minimization_type in ['amber_all', 'amber_h']:
      command_build.append('minimization_options="maxcyc=2"')
    else:
      command_build.append('minimization_options="max_iterations=2"')
    output_build = subprocess.check_output(command_build)
    assert expected_line_dict[minimization_type] in output_build

@pytest.mark.parametrize('minimization_type', ['phenix_all', 'amber_h', 'amber_all'])
@pytest.mark.medium
def test_writing_example_after_running_AmberPrep_not_LES_with_minimization(minimization_type):
  pdb_file = get_fn('2igd/2igd.pdb')

  expected_line_dict = {
          'amber_h': """
==================================================
Done.  Three new files have been made:
      4phenix_2igd.min.amber_h.pdb
      4amber_2igd.prmtop
      4amber_2igd.min.amber_h.rst7
==================================================
""",
          'amber_all': """
==================================================
Done.  Three new files have been made:
      4phenix_2igd.min.amber_all.pdb
      4amber_2igd.prmtop
      4amber_2igd.min.amber_all.rst7
==================================================
""",
          'phenix_all': """
==================================================
Done.  Three new files have been made:
      4phenix_2igd.min.phenix_all.pdb
      4amber_2igd.prmtop
      4amber_2igd.rst7
==================================================
"""
  }
  command_build = [
          'phenix.AmberPrep',
          pdb_file,
  ]
  with tempfolder():
    command_build.append('minimise={}'.format(minimization_type))
    if minimization_type in ['amber_all', 'amber_h']:
      command_build.append('minimization_options="maxcyc=2"')
    else:
      command_build.append('minimization_options="max_iterations=2"')
    output_build = subprocess.check_output(command_build)
    assert expected_line_dict[minimization_type] in output_build
