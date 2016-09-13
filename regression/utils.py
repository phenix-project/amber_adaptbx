import os
from contextlib import contextmanager
import libtbx.load_env
import libtbx
import tempfile
from shutil import rmtree
import subprocess
from amber_adaptbx import is_prmtop_LES
import parmed as pmd
from numpy.testing import assert_almost_equal as aa_eq

MDIN_TEMPLATE = """test minimization
&cntrl
 imin=1, maxcyc={maxcyc}, ntpr=10, ntxo=1,
 /
 &ewald    verbose=1   /
"""

@contextmanager
def tempfolder():
  """run everything in temp folder
  """
  my_temp = tempfile.mkdtemp()
  cwd = os.getcwd()
  os.chdir(my_temp)
  yield
  os.chdir(cwd)
  rmtree(my_temp)

def get_fn(basename):
  # basename is relative to amber_adaptbx/regression/files/ folder
  # e.g get_fn('4amber_2igd.LES.save.rst7')
  return os.path.join(libtbx.env.dist_path("amber_adaptbx"),
                      'regression',
                      'files',
                      basename)

def get_energy_and_forces(prmtop_file, rst7_file):
  # work for both LES and non-LES topology
  import sander, sanderles
  md_engine = sanderles if is_prmtop_LES(prmtop_file) else sander 
  parm = pmd.load_file(prmtop_file, rst7_file)
  with md_engine.setup(prmtop_file, rst7_file, box=parm.box, mm_options=sanderles.pme_input()):
    ene, force = sanderles.energy_forces() 
  return ene, force

def assert_energy_and_forces(prmtop_file,
                             rst7_file,
                             saved_prmtop_file,
                             saved_rst7_file):
  """make sure to reproduce energy and force from saved prmtop
  """
  ene, forces = get_energy_and_forces(prmtop_file, rst7_file)
  saved_ene, saved_forces = get_energy_and_forces(saved_prmtop_file, saved_rst7_file)
  aa_eq(forces, saved_forces)
  aa_eq([ene.tot], [saved_ene.tot])

def run_sander_minimization(prmtop_file, rst7_file, maxcyc=2):
  """run sander minimization and return the content of mdout

  Notes: It's better to be used with tempfolder. prmtop_file can be either normal
  or LES topology

  Examples
  --------
  >>> with tempfolder():
  ...  mdout_content = run_sander_minimization('my.parm7', 'my.rst7', maxcyc=2)
  """
  command = 'sander -O -i mdin -p {0} -c {1}'.format(prmtop_file, rst7_file)
  if is_prmtop_LES(prmtop_file):
    command = command.replace('sander', 'sander.LES')
  with open('mdin', 'w') as mdin:
    mdin.write(MDIN_TEMPLATE.format(maxcyc=maxcyc))
  subprocess.check_call(command.split())
  with open('mdout') as fh:
    return fh.read()
