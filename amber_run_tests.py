# LIBTBX_SET_DISPATCHER_NAME amber.run_tests
# LIBTBX_SET_DISPATCHER_NAME phenix.amber.run_tests

import os, sys
from libtbx import easy_run
import libtbx.load_env

def main():
  # require: pytest
  # Do not want to run in source folder?
  # mkdir $HOME/TMP (or anywhere)
  # cd $HOME/TMP
  # phenix.amber.run_tests
  try:
    import pytest
    print('pytest version = ', pytest.__version__)
  except ImportError:
    raise ImportError("amber_adaptbx requires pytest")
  test_folder = os.path.join(libtbx.env.dist_path("amber_adaptbx"),
                'tests')
  run_command = [
    'phenix.python', '-m', 'pytest',
    test_folder,
    '-v',
    '-k', '"not slow"',
    '-p', 'no:cacheprovider'
  ]
  easy_run.call(' '.join(run_command))

if __name__=="__main__":
  sys.exit(main())
