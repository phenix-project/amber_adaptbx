# LIBTBX_SET_DISPATCHER_NAME elbow.amber.tleap
import os, sys
import StringIO

from iotbx import pdb

from libtbx import easy_run
from libtbx.utils import Sorry

his_names = ["HIS", "HID", "HIE", "HIP"]

changing_residues = []
changing_residues += his_names

def write_tleap_cmd_file(pdb_filename, tleap_input="tleap.in"):
  cmd_template = """
source leaprc.ff12SB
loadamberparams frcmod.ionsjc_tip3p
loadAmberParams frcmod.tip3pf

loadAmberPrep nitrate.prepin
loadAmberParams nitrate.frcmod

loadAmberPrep acetate.prepin

x = loadPdb "%s"

bond x.6.SG x.127.SG
bond x.30.SG x.115.SG
bond x.64.SG x.80.SG
bond x.76.SG x.94.SG

set x box {27.240   31.870   34.230}
set default nocenter on
saveAmberParm x %s.prmtop %s.rst7
quit
"""
  cmd = cmd_template % (
    pdb_filename,
    pdb_filename,
    pdb_filename,
    )
  f=file(tleap_input, "wb")
  f.write(cmd)
  f.close()

def run_tleap_cmd_file(tleap_input):
  fatal_errors = [
    "Failed to generate parameters",
    "Unknown residue",
    ]
  non_fatal_errors = [
    "Could not open file",
    "duplicate [",
    ]
  cmd = "tleap -f %s" % tleap_input
  ero = easy_run.fully_buffered(command=cmd)
  err = StringIO.StringIO()
  ero.show_stdout(out=err)
  outl = ""
  for line in err.getvalue().split("\n"):
    for fe in non_fatal_errors:
      if line.find(fe)>-1:
        print line
    for fe in fatal_errors:
      if line.find(fe)>-1:
        raise Sorry(line)

def run_tleap(pdb_filename):
  tleap_input = "%s.in" % pdb_filename
  write_tleap_cmd_file(pdb_filename, tleap_input)
  run_tleap_cmd_file(tleap_input)

def adjust_his(atom_group):
  hydrogens = ["ND1", "NE2"]
  found = []
  for atom in atom_group.atoms():
    if atom.name.strip() in hydrogens:
      found.append(atom.name.strip())
  if len(found)==2:
    new_residue = "HIP"
  elif len(found)==1:
    if "ND1" in found:
      new_residue = "HID"
    elif "NE2" in found:
      new_residue = "HIE"
    else: assert 0
  else: assert 0
  atom_group.resname = new_residue

def convert_to_amber(hierarchy):
  # assert hydrogens
  for atom_group in hierarchy.atom_groups():
    if atom_group.resname in changing_residues:
      if atom_group.resname in his_names:
        adjust_his(atom_group)

def generate_altloc_filenames(hierarchy):
  conformer_hierarchy = hierarchy.deep_copy()
  for model in conformer_hierarchy.models():
    for chain in model.chains():
      conformers = chain.conformers()
      break
  for i in range(len(conformers)):
    altloc_hierarchy = hierarchy.deep_copy()
    for residue_group in altloc_hierarchy.residue_groups():
      blank = 0
      remove = []
      for j, atom_group in enumerate(residue_group.atom_groups()):
        if j==0 and not atom_group.altloc.strip():
          blank=1
          continue
        if i+blank!=j:
          remove.append(atom_group)
      if remove:
        for r in remove:
          residue_group.remove_atom_group(r)
    output_filename = "tleap_input_%d.pdb" % (i+1)
    f=file(output_filename, "wb")
    f.write(altloc_hierarchy.as_pdb_string(
      #crystal_symmetry=pdb_inp.crystal_symmetry()
      )
      )
    f.close()
    yield output_filename

def run(pdb_filename):
  # convert to tleap input PDB
  pdb_inp = pdb.input(pdb_filename)
  hierarchy = pdb_inp.construct_hierarchy()
  convert_to_amber(hierarchy)
  for i, tleap_input_filename in enumerate(generate_altloc_filenames(hierarchy)):
    print tleap_input_filename
    # run tleap
    run_tleap(tleap_input_filename)

if __name__=="__main__":
  run(sys.argv[1])
  
