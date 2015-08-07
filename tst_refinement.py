
from __future__ import division
from libtbx import easy_run
import libtbx.load_env
import os.path

def exercise_vAla3 () :
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="amber_adaptbx/test_vAla3/vAla3.pdb",
    test=os.path.isfile)
  top_file = os.path.splitext(pdb_file)[0] + ".prmtop"
  crd_file = os.path.splitext(pdb_file)[0] + ".rst7"
  cif_file = os.path.splitext(pdb_file)[0] + ".cif"
  args = [
    "phenix.fmodel",
    pdb_file,
    "high_resolution=1.5",
    "type=real",
    "label=F",
    "r_free_flags_fraction=0.1",
    "output.file_name=vAla3.mtz",
  ]
  rc = easy_run.fully_buffered(" ".join(args)).raise_if_errors().return_code
  assert (rc == 0)
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
    "phenix.refine",
    "vAla3_shaken.pdb",
    cif_file,
    "vAla3.mtz",
    "topology_file_name=\"%s\"" % top_file,
    "coordinate_file_name=\"%s\"" % crd_file,
    "use_amber=True",
    "wxc_scale=0.1",
    "--overwrite",
  ]
  rc = easy_run.fully_buffered(" ".join(args)).raise_if_errors().return_code
  assert (rc == 0)
  from phenix.refinement.runtime import extract_phenix_refine_r_factors
  (r_work, r_free) = extract_phenix_refine_r_factors(
    file_name="vAla3_shaken_refine_001.pdb")
  assert r_work.startswith("0.00") and r_free.startswith("0.00"), """
  r_work : %s != 0.00
  r_free : %s != 0.00
  """ % (r_work, r_free)

  remove_files = ['vAla3.mtz', 'vAla3_shaken.pdb', 
                  'vAla3_shaken_refine_001.eff',
                  'vAla3_shaken_refine_001.geo',
                  'vAla3_shaken_refine_001.log',
                  'vAla3_shaken_refine_001.mtz',
                  'vAla3_shaken_refine_001.pdb',
                  'vAla3_shaken_refine_002.def']
  for file in remove_files:
    if os.path.exists(file):
      os.remove(file)                  

if (__name__ == "__main__") :
  exercise_vAla3()
  print "OK"
