
import libtbx.load_env
from libtbx.utils import Sorry
import os.path
import sys


if os.path.isdir('/net/casegroup2/u2/pjanowsk/bin/phenix_bootstrap/modules/amber'):
  print "AMBER SCONS: Amber linked. Will attempt MDGX compile."
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
  env_amber_ext.Append(CPPDEFINES='AMBERPHENIX')
  env_amber_ext.Append(LIBS=["fftw3", "netcdf", "mdgx",])
  env_amber_ext.Append(LIBPATH=["%s/lib" % os.environ["AMBERHOME"]])
  env_etc.include_registry.append(
    env=env_amber_ext,
    paths=env_etc.amber_common_includes)
  Export("env_amber_ext")
  env_amber_ext.SharedLibrary(
    target="#lib/amber_adaptbx_ext",
    source=ext_sources)

else:
  print "AMBER SCONS: no Amber link. Only Sander will be available."