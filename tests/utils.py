import os
import math
from contextlib import contextmanager
import tempfile
from shutil import rmtree
import subprocess
from amber_adaptbx import is_prmtop_LES
import parmed as pmd
import numpy as np
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
  # basename is relative to amber_adaptbx/tests/files/ folder
  # e.g get_fn('4amber_2igd.LES.save.rst7')
  fn = os.path.join(os.path.dirname(__file__),
                      'files',
                      basename)
  assert os.path.exists(fn), 'File must exists {}'.format(fn)
  return fn

def rmsd(a1, a2):
  """rmsd for two array with the same shape.
  Arrays will be flattened

  Parameters
  ----------
  a1, a2: np.ndarray
  flatten : bool, default True
      if True: always flatten two input arrays
  """
  a1 = np.asarray(a1).flatten()
  a2 = np.asarray(a2).flatten()
  tmp = sum((a1- a2)**2)
  return math.sqrt(tmp / a1.shape[0])

def get_prmtop_and_rst7_and_pdb_filenames_from_pdb(pdb_file, LES=False):
  """e.g: get_prmtop_and_rst7_and_pdb_filenames_from_pdb('my_pdb.pdb') will return '4amber_my_pdb.prmtop'
  and '4amber_2igd.rst7'
  """
  basename = os.path.basename(pdb_file)
  root = basename.split('.')[0]
  if LES:
    return ('4amber_' + root + '.LES.prmtop',
            '4amber_' + root + '.LES.rst7',
            '4phenix_' + root + '.LES.pdb')
  else:
    return ('4amber_' + root + '.prmtop',
            '4amber_' + root + '.rst7',
            '4phenix_' + root + '.pdb')

def get_minimized_pdb_filename(pdb_file, minimization_type, LES=False):
  basename = os.path.basename(pdb_file)
  root = basename.split('.')[0]
  if LES:
    return '4phenix_' + root + '.LES.min.{}.pdb'.format(minimization_type)
  else:  
    return '4phenix_' + root + '.min.{}.pdb'.format(minimization_type)

def get_minimized_rst7_filename(pdb_file, minimization_type, LES=False):
  basename = os.path.basename(pdb_file)
  root = basename.split('.')[0]
  if LES:
    return '4amber_' + root + '.LES.min.{}.rst7'.format(minimization_type)
  else:  
    return '4amber_' + root + '.min.{}.rst7'.format(minimization_type)

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
  aa_eq(forces, saved_forces, decimal=4)
  aa_eq([ene.tot], [saved_ene.tot], decimal=4)

def assert_file_has_line(filename, line):
  with open(filename) as fh:
    assert line in fh.read(), '{} must has {}'.format(filename, line)

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

if __name__ == '__main__':
  assert(get_minimized_pdb_filename('2igd.pdb', minimization_type='amber_h') ==
         '4phenix_2igd.min.amber_h.pdb')
  assert(get_minimized_pdb_filename('2igd.pdb', LES=True, minimization_type='amber_h') ==
         '4phenix_2igd.LES.min.amber_h.pdb')
  assert(get_minimized_rst7_filename('2igd.pdb', minimization_type='amber_h') ==
         '4amber_2igd.min.amber_h.rst7')
  assert(get_minimized_rst7_filename('2igd.pdb', LES=True, minimization_type='amber_h') ==
         '4amber_2igd.LES.min.amber_h.rst7')
  assert (get_prmtop_and_rst7_and_pdb_filenames_from_pdb('fake_dir/my.pdb') ==
          ('4amber_my.prmtop', '4amber_my.rst7', '4phenix_my.pdb'))
  assert (get_prmtop_and_rst7_and_pdb_filenames_from_pdb('fake_dir/my.pdb', LES=True) ==
          ('4amber_my.LES.prmtop', '4amber_my.LES.rst7', '4phenix_my.LES.pdb'))
