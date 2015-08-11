
import libtbx.load_env
from libtbx.utils import Sorry
import os.path
import sys


amber_dist = libtbx.env.dist_path("amber", default=None)
amberhome = os.environ.get("AMBERHOME", False)
if(amber_dist and 
   os.path.exists(amber_dist) and
   os.path.isdir(amber_dist) and
   amberhome
   ):
  print "AMBER SCONS: Amber linked. Will attempt MDGX compile."
  Import("env_base", "env_etc")
  ext_sources = [
    "ext.cpp",
    "getmdgxfrc.c",
  ]

  env_etc.amber_adaptbx_dist = libtbx.env.dist_path("amber_adaptbx")
  amber_dir = os.path.join(os.path.dirname(env_etc.amber_adaptbx_dist),
                           "amber")
  amber_src_dir = os.path.join(amber_dir, "src")

  env_etc.amber_common_includes = [
    amber_dir+"/include",
    env_etc.boost_include,
    env_etc.python_include,
    env_etc.libtbx_include,
    env_etc.scitbx_include,
    os.path.dirname(env_etc.amber_adaptbx_dist),
  ]

  Import("env_scitbx_boost_python_ext")
  env_amber_ext = env_scitbx_boost_python_ext.Clone()
  env_amber_ext.Append(CPPDEFINES='AMBERPHENIX')
  env_amber_ext.Append(LIBS=["fftw3", "netcdf", "mdgx",])
  env_amber_ext.Append(LIBPATH=[
    os.path.join("%s" % os.environ["AMBERHOME"], "lib")])
  env_etc.include_registry.append(
    env=env_amber_ext,
    paths=env_etc.amber_common_includes)
  Export("env_amber_ext")
  env_amber_ext.SharedLibrary(
    target="#lib/amber_adaptbx_ext",
    source=ext_sources)

else:
  print "AMBER SCONS: no Amber link. Only Sander will be available."
