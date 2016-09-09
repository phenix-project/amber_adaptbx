import os
from contextlib import contextmanager
import libtbx.load_env
import libtbx
import tempfile
from shutil import rmtree

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
  # basename is relative to amber_adaptbx folder
  # e.g get_fn('test_les/2igd/4amber_2igd.LES.save.rst7')
  return os.path.join(libtbx.env.dist_path("amber_adaptbx"),
       basename)
