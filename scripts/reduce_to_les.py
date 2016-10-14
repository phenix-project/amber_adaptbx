import os
import sys
import argparse
from collections import defaultdict
import parmed as pmd

"""

Examples
--------
phenix.python reduce_to_les.py --asu-pdb=2igd.pdb --uc-pdb=2igd_uc_H.pdb --prmtop=2igdab.parm7 --rst7=2igdab.rst7

Output: 2igdab_fix.rst7, 4phenix_2igdab.pdb (ASU)

Notes
-----
2igd.pdb : original ASP pdb (downloaded from rcsb?)
2igd_uc_H.pdb : unitcell pdb, hydrogend were added by reduduce
2igdab.parm7 : LES parm7
2igdab.rst7 : LES rst7 having coordinates needed to be updated
"""

def main():
  # for command line
  parser = argparse.ArgumentParser(description='fix LES parm7 and rst7 files')
  parser.add_argument('--asu-pdb', help='original ASU pdb')
  parser.add_argument('--uc-pdb', help='unitcell pdb that have added hydrogens (by reduce)')
  parser.add_argument('--prmtop', help='prmtop for LES')
  parser.add_argument('--rst7', help='rst7 file for LES (need to be fixed)')
  args = parser.parse_args()

  uc_parm = pmd.load_file(args.uc_pdb)
  parm = pmd.load_file(args.prmtop, xyz=args.rst7)
  base_name = os.path.basename(args.prmtop).split('.')[0]
  rst7_fn = base_name + '_fix.rst7'

  parm = update_rst7_and_pdb_coordinates_LES(template_parm=uc_parm,
                          target_parm=parm,
                          rst7_fn=rst7_fn)

  # update symmetry and space group info
  asu_parm = pmd.load_file(args.asu_pdb)
  parm.box = asu_parm.box

  n_asu_res = len(asu_parm.residues)
  selection = ':1-{}'.format(str(n_asu_res))

  # pick the 1st ASU and correctly assing symmetry, space_group
  les_asu_parm = parm[selection]
  les_asu_parm.box = asu_parm.box
  les_asu_parm.symmetry = asu_parm.symmetry
  les_asu_parm.space_group = asu_parm.space_group
  # write ASU pdb with the same atom order as prmtop
  les_asu_parm.save('4phenix_{}.pdb'.format(base_name), overwrite=True)

  # write correct LES rst7 file
  parm.save(rst7_fn, overwrite=True)

def get_atom_dict_for_uc(uc_parm):
  """

  Parameters
  ----------
  uc_parm : parmed's Structure, read from unitcell pdb file
      the pdb file must have alternatve conformations
  """
  atom_dict = defaultdict(list)
  
  resids = range(len(uc_parm.residues))

  for resid in set(resids):
    # include atom in conformation A
    atom_dict[resid] = defaultdict(list)
    for atom in uc_parm.residues[resid].atoms:
      atom_dict[resid][atom.name].append(atom)

      if atom.other_locations:
        for alt, alt_atom in atom.other_locations.items():
          atom_dict[resid][atom.name].append(alt_atom)

  return atom_dict

def get_atom_dict_for_amber_parm(parm, resids):
  """

  Parameters
  ----------
  parm : parmed's Structure, read from AMBER's prmtop
  resids: list of residues having alternate atoms
      keys can be computed from get_atom_dict_for_uc method
  """
  atom_dict = defaultdict(list)

  for resid in resids:
    atom_dict[resid] = defaultdict(list)
    residue = parm.residues[resid]

    for atom in residue.atoms:
      atom_dict[resid][atom.name].append(atom)

  return atom_dict

def update_rst7_and_pdb_coordinates_LES(template_parm, target_parm, rst7_fn):
  '''

  Parameters
  ----------
  template_parm : parmed.Structure, created from unit cell pdb that has alternate atoms
      template_parm must have full atoms (including H)
  parm : parmed.Structure created from parm7 and rst7 files
  '''
  # dict of residue and alt atoms
  atom_dict_template = get_atom_dict_for_uc(template_parm)
  alt_resids = atom_dict_template.keys()
  atom_dict_amber_parm = get_atom_dict_for_amber_parm(target_parm, alt_resids)

  for resid in alt_resids:
    for atom_name in atom_dict_template[resid].keys():
      template_atoms = atom_dict_template[resid][atom_name]
      target_atoms = atom_dict_amber_parm[resid][atom_name]

      for atom, template_atom in zip(target_atoms, template_atoms):
        atom.xx = template_atom.xx
        atom.xy = template_atom.xy
        atom.xz = template_atom.xz
        atom.bfactor = template_atom.bfactor
        atom.occupancy = template_atom.occupancy

  pdb_fn = rst7_fn.split('.')[0] + '.pdb'

  # TODO: already updated for ParmEd, will remove
  label_alternates(target_parm)
  return target_parm

def label_alternates(parm):
  """

  Parameters
  ----------
  parm : parmed's Structure
      parm should have duplicated atom name. This function will
      assign label A, B, ...

  Returns
  -------
  out : updated `parm`

  Examples
  --------
  >>> parm = label_alternates(parm)
  >>> parm.write_pdb("corrected_label.pdb", altlocs='all', standard_resnames=True)
  """
  resids = range(len(parm.residues))
  atom_with_residue_dict = get_atom_dict_for_amber_parm(parm, resids=resids)

  possible_labels = list('ABCDEF')

  for _, adict in atom_with_residue_dict.items():
    for atom_name, atom_list in adict.items():
      if len(atom_list) > 1:
        # if there is alternate atom
        for atom, label in zip(atom_list, possible_labels):
          atom.altloc = label

  return parm

if __name__ == '__main__':
  main()
