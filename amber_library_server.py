import os, sys

import libtbx.load_env

parent_dir = os.path.dirname(libtbx.env.dist_path("elbow"))

def repo_dir():
  env_dir = os.environ.get("AMBER_LIBRARY_DIR", None)
  if env_dir is not None:
    return env_dir
  install_dir = os.path.join(parent_dir, "chem_data", "amberlibrary")
  if os.path.exists(install_dir):
    return install_dir
  return None

def is_in_components_lib(residue_name):
  if path_in_components_lib(residue_name): return True
  else: return False

def path_in_components_lib(residue_name):
  rd = repo_dir()
  preamble = residue_name = os.path.join(rd,
                                         residue_name[0].lower(),
                                         residue_name.upper(),
                                         )
  files = []
  for ext in ["frcmod", "mol2"]:
    if not os.path.exists("%s.%s" % (preamble, ext)):
      break
    files.append("%s.%s" % (preamble, ext))
  else:
    return files
  return False

def run():
  print 'Repo directory',repo_dir()
  for code in ["000",
               "HOH",
               "NWM",
               ]:
    print code, is_in_components_lib(code), path_in_components_lib(code)
  
if __name__=="__main__":
  run()#sys.argv[1])
  
