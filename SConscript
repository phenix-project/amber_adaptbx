
import libtbx.load_env
from libtbx.utils import Sorry
import os.path
import sys

Import("env_base", "env_etc")
ext_sources = [
  "ext.cpp",
  "getmdgxfrc.c",
]

env_etc.amber_adaptbx_dist = libtbx.env.dist_path("amber_adaptbx")
amber_dir = os.path.join(os.path.dirname(env_etc.amber_adaptbx_dist),
  "amber")
amber_src_dir = amber_dir + "/src"

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
#env_rosetta_ext.Append(LIBPATH=["%s/lib" % rosetta_dir])
env_amber_ext.Append(LIBPATH=env_etc.libpath_python)
env_amber_ext.Append(CCFLAGS=["-I%s" % amber_dir])
env_amber_ext.Append(LIBS=["fftw3", "netcdf", "mdgx",])
env_etc.include_registry.append(
  env=env_amber_ext,
  paths=env_etc.amber_common_includes)
Export("env_amber_ext")
env_amber_ext.SharedLibrary(
  target="#lib/amber_adaptbx_ext",
  source=ext_sources)
