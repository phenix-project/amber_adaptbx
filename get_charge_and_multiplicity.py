import os, sys

from iotbx import pdb

from elbow.command_line import builder

def generate_unknown_residue_names(pdb_filename):
  pdb_inp = pdb.input(pdb_filename)
  hierarchy = pdb_inp.construct_hierarchy()
  for atom_group in hierarchy.atom_groups():
    yield atom_group.resname

def get_charge_and_multiplicity(pdb_filename):
  results = {}
  for residue_name in generate_unknown_residue_names(pdb_filename):
    print residue_name
    if residue_name in results: continue
    mol = builder.run(chemical_component=residue_name,
                      no_output=True,
                      )
    mol.Multiplicitise()
    print mol.DisplayBrief()
    results[residue_name] = [mol.charge, mol.multiplicity]
  return results

def run(pdb_filename):
  print pdb_filename
  print get_charge_and_multiplicity(pdb_filename)

if __name__=="__main__":
  run(sys.argv[1])
