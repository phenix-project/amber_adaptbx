import sys
import argparse
import parmed as pmd

"""
Examples
--------

  phenix.python make_4phenix_LES_pdb.py --asu-pdb=2igd.pdb  --prmtop=2igdab.parm7 --rst7=2igdab.min1.rst7 --output-pdb=4phenix_2igdab.pdb
"""

def main():
  parser = argparse.ArgumentParser(description='convert LES parm7 and rst7 to ASU pdb for phenix')
  parser.add_argument('--asu-pdb', help='original ASU pdb to give info about space group and number or residues in ASU')
  parser.add_argument('--prmtop', help='prmtop for LES')
  parser.add_argument('--rst7', help='rst7 file for LES')
  parser.add_argument('--output-pdb', help='pdb output name')
  args = parser.parse_args()

  asu_parm = pmd.load_file(args.asu_pdb)
  n_asu_residue = len(asu_parm.residues)
  selection = ':1-{}'.format(n_asu_residue)

  parm = pmd.load_file(args.prmtop, xyz=args.rst7)
  # strip atoms
  asu_new_parm = parm[selection]
  # need symmetry too? seems not
  # asu_new_parm.symmetry = asu_parm.symmetry
  asu_new_parm.space_group = asu_parm.space_group
  asu_new_parm.box = parm.box
  asu_new_parm.symmetry = asu_parm.symmetry
  # update occupancy
  # we can use original occupancy from ASU pdb file too.
  for atom in asu_new_parm.atoms:
      atom.occupancy = 1.0
  asu_new_parm.write_pdb(args.output_pdb, standard_resnames=True)

if __name__ == '__main__':
  main()
