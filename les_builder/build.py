import os
import iotbx.pdb
import parmed
from amber_adaptbx.les_builder import make_addles_input
from libtbx import easy_run
from amber_adaptbx.les_builder import reduce_to_les
from amber_adaptbx.utils import build_unitcell
from parmed.utils.six.moves import StringIO

import argparse

__all__ = ['LESBuilder']

"""
Command line
------------

  phenix.python les_build.py 2igd.pdb
  phenix.python les_build.py --help

Output example
--------------
4amber_2igd.LES.prmtop  4amber_2igd.LES.rst7  4phenix_2igd.LES.pdb

API
---

Check the example below.

Notes
-----

"""

EXPECTED_HEADER_TEMPLATE = """
file rprm name=(4amber_{base}.prmtop) read
file rcbd name=(4amber_{base}.rst7) read
file wprm name=(4amber_{base}.LES.prmtop) wovr
file wcrd name=(4amber_{base}.LES.rst7) wovr
""".strip()

class LESBuilder(object):
  """Class to build LES parm, rst7 and pdb file.

  Parameters
  ----------
  original_pdb_file : str, ASU pdb (e.g downloaded from rcsb)
  prmtop            : str, prmtop of the single conformation unitcell
                      (You can use AmberPrep to prepare this prmtop file)
  rst7_file         : str, corresponding rst7 filename
  unitcell_pdb_file : str or None, default None
    fully atom unitcell pdb filename. If given, LESBuilder will use it. 
    If not, it will call build_unitcell method from amber_adaptbx.utils.
  addles_input_file  : str, default ''
    addles input filename. If given, use this file instead of making a new one

  Examples
  --------
  >>> # build LES for pdb that have ligands
  >>> # use phenix.AmberPrep to build 2igda.parm7, 2igda.rst7
  >>> builder = LESBuilder('./2igd.pdb', prmtop='./2igda.parm7', rst7_file='./2igda.rst7')
  >>> builder.run()
  """

  def __init__(self, original_pdb_file, prmtop=None, rst7_file=None, unitcell_pdb_file=None, addles_input_file=''):
    self.original_pdb_file = original_pdb_file
    self.new_pdb_with_H = None

    # prmtop and rst7_file were built from phenix.AmberPrep.
    self.prmtop = prmtop
    self.rst7_file = rst7_file

    self.base = '.'.join(os.path.basename(
        self.original_pdb_file).split('.')[:-1])
    self._orig_structure = parmed.load_file(original_pdb_file)
    self.space_group = self._orig_structure.space_group
    self.symmetry = self._orig_structure.symmetry
    self.box = self._orig_structure.box
    self.n_asu_residues = len(self._orig_structure.residues)
    self.unitcell_pdb_file = unitcell_pdb_file
    self.addles_input_file = addles_input_file

  def build_LES_parm(self):
    # create addles.in, then run addles with this input.
    #   creates xxxxab.prmtop and xxxxab.rst7
    print "\n============================================================"
    print " Building the prmtop and rst7 file with alternate conformers"
    print "============================================================"
    if not self.addles_input_file:
      commands = make_addles_input.addles_input(self.original_pdb_file,
                                                self.prmtop,
                                                self.rst7_file)
      fn = 'addles.in'
      with open(fn, 'w') as fh:
        fh.write('\n'.join(commands))
    else:
      fn = self.addles_input_file
      self._check_valid_addles_input(fn)
      print('\n--> Using provided addles input: {}'.format(fn))
    command_addles = 'addles <{} > addles.log'.format(fn)
    print "\n~> %s\n" % command_addles
    easy_run.fully_buffered(command_addles)

  def update_LES_coordinates_from_uc(self):
    """label_alternates + update LES coordinates
    """
    uc_parm = parmed.load_file(self.new_pdb_with_H)
    parm = parmed.load_file('4amber_' + self.base + '.LES.prmtop',
                            '4amber_' + self.base + '.LES.rst7')
    # add space group
    parm.space_group = self.space_group
    parm.symmetry = self.symmetry
    parm.box = self.box
    reduce_to_les.update_rst7_and_pdb_coordinates_LES(template_parm=uc_parm,
                                                      target_parm=parm)
    parm.write_pdb('4amber_%s.LES.pdb' % self.base, standard_resnames=True)
    parm.save('4amber_%s.LES.rst7' % self.base, overwrite=True)

  def write_LES_asu_pdb(self):
    """
    write LES asu pdb for phenix: 4phenix_xxxx.LES.pdb
    """
    les_parm = parmed.load_file('4amber_%s.LES.pdb' % self.base)
    orig_pdb_parm = parmed.load_file(self.original_pdb_file)

    selection = ':1-{}'.format(self.n_asu_residues)

    final_pdb_asu_file = '4phenix_' + self.base + '.LES.pdb'
    final_parm = les_parm[selection]
    final_parm.symmetry = self.symmetry
    final_parm.box = self.box
    final_parm.space_group = self.space_group

    # update H occupancy and bfactors:
    for atom in final_parm.atoms:
      if atom.atomic_number == 1:
        atom.occupancy = atom.bond_partners[0].occupancy
        atom.bfactor = atom.bond_partners[0].bfactor

    # restore original residue numbers
    for residue, template_residue in zip(final_parm.residues, 
                                         orig_pdb_parm.residues):
      residue.number = template_residue.number
      residue.chain = template_residue.chain

    final_parm.write_pdb( final_pdb_asu_file,
                          standard_resnames=True, renumber=False )

    # run through phenix heierarchy to get atoms sorted phenix-style:
    pdb_input = iotbx.pdb.input(final_pdb_asu_file)
    pdb_hierarchy = pdb_input.construct_hierarchy(sort_atoms=True)
    pdb_hierarchy.write_pdb_file(file_name=final_pdb_asu_file,
            append_end=True,
            crystal_symmetry=pdb_input.crystal_symmetry() )

  def _check_valid_addles_input(self, fn):
    assert os.path.exists(fn), 'make sure {} exists'.format(fn)
    input = open(fn).read()
    expected_header = EXPECTED_HEADER_TEMPLATE.format(base=self.base)
    assert expected_header in input, 'addles input must have header\n\n{}\n'.format(expected_header)

  def _has_altlocs(self):
    for residue in self._orig_structure.residues:
      for atom in residue.atoms:
        if atom.other_locations:
          return True
    return False

  def run(self):

    if not self._has_altlocs():
      raise ValueError("pdb file should have altlocs for LES build")

    # build unitcell from asu pdb============================================
    #    (input is usually original pdb file; creates xxxx_uc.pdb)
    #    note: this file won't have missing atoms...
    if self.unitcell_pdb_file is None:
      self.unitcell_pdb_file = self.base + '_uc.pdb'
      build_unitcell(self.original_pdb_file, self.unitcell_pdb_file)

      # use reduce to add hydrogens to unitcell pdb==========================
      #   (input is xxxx_uc.pdb; creates xxxx_uc_H.pdb)
      self.new_pdb_with_H = self.base + '_uc_H.pdb'
      cmd = 'reduce -build -nuclear {} > {} 2>reduce_lesbuilder.log'.format(
          self.unitcell_pdb_file, self.new_pdb_with_H)
      print "\n~> %s\n" % cmd
      easy_run.fully_buffered(cmd)

    # use addles to construct LES parm7 and rst7 files========================
    #   (input is 4amber_xxxx.{prmtop,rst7};  output is
    #             4amber_xxxx.LES.{prmtop,rst7}
    self.build_LES_parm()

    # fix previous step: need to update rst7 coordinates since addles
    # uses the same coordinates for alternative atoms.
    # Also correctly label atom and residue names.============================
    #   (updates 4amber_xxxx.LES.rst7 in place; also creates 
    #    4amber_xxxx.LES.pdb for use in the next step)
    self.update_LES_coordinates_from_uc()

    # generate final ASU LES pdb for phenix (with unitcell, symmetry)========
    #   (input is 4amber_xxxx.LES.pdb; creates 4phenix_xxxx.LES.pdb )
    self.write_LES_asu_pdb()

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
      description='build LES parm7 and rst7 files from ASU pdb.')
  parser.add_argument('pdb', help='ASU pdb')
  parser.add_argument('--prmtop', help='prmtop for ASU (optional)')
  parser.add_argument('--rst7', help='rst7 file for ASU (optional)')
  args = parser.parse_args()
  les_builder = LESBuilder(args.pdb, prmtop=args.prmtop, rst7_file=args.rst7)
  les_builder.run()
