from __future__ import print_function
import os, sys
import libtbx.load_env

parent_dir = os.path.dirname(libtbx.env.dist_path("amber_adaptbx"))

def is_energy_outlier(residue_name):
  rd = repo_dir()
  if rd is None:
    return None
  outliers = os.path.join(rd, "outliers_min_energy.dat")
  if not os.path.exists(outliers): return False
  f = file(outliers, "rb")
  lines = f.readlines()
  f.close()
  for line in lines:
    if line.find(residue_name.upper()) != -1:
      return True
  return False

def repo_dir(verbose=False):
  env_dir = os.environ.get("AMBER_LIBRARY_DIR", None)
  if env_dir is not None:
    return env_dir
  install_dirs = []
  install_dir = os.path.join(parent_dir, "amber_library")
  if os.path.exists(install_dir): install_dirs.append(install_dir)
  install_dir = os.path.join(parent_dir, 'chem_data', 'geostd')
  if os.path.exists(install_dir): install_dirs.append(install_dir)
  if install_dirs:
    return install_dirs
  if verbose:
    print("""
    Couldn't find amber_library
      1. Set AMBER_LIBRARY_DIR in environment
      2. Add/link to $PHENIX/modules
    """)
  return None

def is_in_components_lib(residue_name):
  rc = path_in_components_lib(residue_name)
  if rc:
    return rc
  else:
    return False

def path_in_components_lib(residue_name):
  rds = repo_dir()
  if rds is None:
    return None
  # if is_energy_outlier(residue_name):
  #   return 0
  for rd in rds:
    preamble = os.path.join(rd,
                            residue_name[0].lower(),
                            residue_name.upper(),
                            )
    files = []
    for ext in ["frcmod", "mol2", "lib"]:
      af = "%s.%s" % (preamble, ext)
      if os.path.exists(af):
        files.append(af)
      af = os.path.join(os.path.dirname(af),
                        'data_%s' % os.path.basename(af))
      if os.path.exists(af):
        files.append(af)
    if len(files) == 2:
      return files
  return False

def run(only_code=None):
  print('Repo directory', repo_dir())
  answers = ['geostd','geostd','geostd', False, 'geostd','amber_library']
  for i, code in enumerate([
    "000",
    "HOH",
    "NWM",
    "NUC",
    'AUX',
    'CIT',
    ]):
    if only_code: code=only_code
    iic = is_in_components_lib(code)
    print(code, iic, path_in_components_lib(code))
    if only_code: break
    if answers[i]:
      assert iic[0].find(answers[i])>-1
    else:
      assert not iic

if __name__ == "__main__":
  run(*tuple(sys.argv[1:]))
