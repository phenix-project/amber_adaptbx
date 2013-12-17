
from __future__ import division
from libtbx import easy_run
import libtbx.load_env
import os.path
import iotbx.pdb

def exercise_vAla3 () :
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="amber_adaptbx/vAla3/vAla3_minimized.pdb",
    test=os.path.isfile)
  top_file = libtbx.env.find_in_repositories(
    relative_path="amber_adaptbx/vAla3/vAla3.prmtop",
    test=os.path.isfile)
  crd_file = os.path.splitext(top_file)[0] + ".rst7"
  cif_file = os.path.splitext(top_file)[0] + ".cif"

  args = [
    "phenix.pdbtools",
    pdb_file,
    cif_file,
    "output.file_name=vAla3_shaken.pdb",
    "sites.shake=0.1",
  ]
  rc = easy_run.fully_buffered(" ".join(args)).raise_if_errors().return_code
  assert (rc == 0)

  args = [
    "phenix.geometry_minimization",
    "vAla3_shaken.pdb",
    cif_file,
    "topology_file_name=\"%s\"" % top_file,
    "coordinate_file_name=\"%s\"" % crd_file,
    "use_amber=True"
  ]
  rc = easy_run.fully_buffered(" ".join(args)).raise_if_errors().return_code
  assert (rc == 0)

  pdb_inp = iotbx.pdb.input(file_name=pdb_file)
  symm = pdb_inp.crystal_symmetry()
  xray_structure = pdb_inp.xray_structure_simple(enable_scattering_type_unknown=True)
  sites_cart_inp=xray_structure.sites_cart()

  pdb_inp = iotbx.pdb.input(file_name='vAla3_shaken.pdb')
  symm = pdb_inp.crystal_symmetry()
  xray_structure = pdb_inp.xray_structure_simple(enable_scattering_type_unknown=True)
  sites_cart_shaken=xray_structure.sites_cart()

  pdb_inp = iotbx.pdb.input(file_name='vAla3_shaken_minimized.pdb')
  symm = pdb_inp.crystal_symmetry()
  xray_structure = pdb_inp.xray_structure_simple(enable_scattering_type_unknown=True)
  sites_cart_min=xray_structure.sites_cart()

  assert sites_cart_inp.rms_difference(sites_cart_shaken) >0.09
  assert sites_cart_inp.rms_difference(sites_cart_min) <0.06, \
    "RMSD of amber-minimized structure is %5.4f (<0.06 required to pass)." \
    % sites_cart_inp.rms_difference(sites_cart_min)


if (__name__ == "__main__") :
  exercise_vAla3()
  print "OK"
