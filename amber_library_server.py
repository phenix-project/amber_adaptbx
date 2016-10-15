import os
import libtbx.load_env

parent_dir = os.path.dirname(libtbx.env.dist_path("amber_adaptbx"))


def is_energy_outlier(residue_name):
  rd = repo_dir()
  if rd is None:
    return None
  outliers = os.path.join(rd, "outliers_min_energy.dat")
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
  install_dir = os.path.join(parent_dir, "chem_data", "amber_library")
  if os.path.exists(install_dir):
    return install_dir
  if verbose:
    print """
    Couldn't find amberlibrary
      1. Set AMBER_LIBRARY_DIR in environment
      2. Add/link to $PHENIX/chem_data
    """
  return None


def is_in_components_lib(residue_name):
  rc = path_in_components_lib(residue_name)
  if rc:
    return rc
  else:
    return False


def path_in_components_lib(residue_name):
  rd = repo_dir()
  if rd is None:
    return None
  if is_energy_outlier(residue_name):
    return 0
  preamble = os.path.join(rd,
                          residue_name[0].lower(),
                          residue_name.upper(),
                          )
  files = []
  for ext in ["frcmod", "mol2"]:
    if not os.path.exists("%s.%s" % (preamble, ext)):
      break
    files.append("%s.%s" % (preamble, ext))
  else:
    assert len(files) == 2
    return files
  return False


def run():
  print 'Repo directory', repo_dir()
  for code in ["000",
               "HOH",
               "NWM",
               "NUC",
               ]:
    print code, is_in_components_lib(code), path_in_components_lib(code)

if __name__ == "__main__":
  run()
