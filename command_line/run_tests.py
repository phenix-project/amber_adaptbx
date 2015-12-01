import os, sys

from libtbx import easy_run
import libtbx.load_env

def run():
  print 'Running Amber tests'
  print '  minimisation...'
  cmd = "libtbx.python %s" % (os.path.join(
    libtbx.env.dist_path("amber_adaptbx"),
    "tst_minimization.py"))
  print cmd
  easy_run.call(cmd)
  print '  refinement...'
  cmd = "libtbx.python %s" % (os.path.join(
    libtbx.env.dist_path("amber_adaptbx"),
    "tst_refinement.py"))
  print cmd
  easy_run.call(cmd)

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
