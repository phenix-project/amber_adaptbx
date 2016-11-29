import os
import iotbx.pdb
import parmed as pmd
from amber_adaptbx.les_builder import make_addles_input
from libtbx import easy_run
from amber_adaptbx.les_builder import reduce_to_les
from amber_adaptbx.utils import build_unitcell, write_standard_pdb
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
file rprm name=(4amber_{root_name}.prmtop) read
file rcbd name=(4amber_{root_name}.rst7) read
file wprm name=(4amber_{root_name}.LES.prmtop) wovr
file wcrd name=(4amber_{root_name}.LES.rst7) wovr
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

    self.root_name = '.'.join(os.path.basename(
        self.original_pdb_file).split('.')[:-1])
    asu_pdb_parm = pmd.load_file(original_pdb_file)
    self.space_group = asu_pdb_parm.space_group
    self.symmetry = asu_pdb_parm.symmetry
    self.box = asu_pdb_parm.box
    self.n_asu_residues = len(asu_pdb_parm.residues)
    self.unitcell_pdb_file = unitcell_pdb_file
    self.addles_input_file = addles_input_file

  def run(self):
    # main driver

    # build unitcell from asu pdb============================================
    #    (input is usually original pdb file; creates xxxx_uc.pdb)
    if self.unitcell_pdb_file is None:
      self.unitcell_pdb_file = self.root_name + '_uc.pdb'
      build_unitcell(self.original_pdb_file, self.unitcell_pdb_file)

      # use reduce to add hydrogens to unitcell pdb==========================
      #   (input is xxxx_uc.pdb; creates xxxx_uc_H.pdb)
      self.new_pdb_with_H = self.root_name + '_uc_H.pdb'
      command_add_hydrogens = 'reduce -build -nuclear {} > {} 2>reduce_lesbuilder.log'.format(
          self.unitcell_pdb_file, self.new_pdb_with_H)
      print "\n~> %s\n" % command_add_hydrogens
      easy_run.fully_buffered(command_add_hydrogens)

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
    uc_parm = pmd.load_file(self.new_pdb_with_H)
    parm = pmd.load_file('4amber_' + self.root_name + '.LES.prmtop',
                         '4amber_' + self.root_name + '.LES.rst7')
    # add space group
    parm.space_group = self.space_group
    parm.symmetry = self.symmetry
    parm.box = self.box
    reduce_to_les.update_rst7_and_pdb_coordinates_LES(template_parm=uc_parm,
                                                      target_parm=parm)
    parm.write_pdb('4amber_{}.LES.pdb'.format(self.root_name), standard_resnames=True)
    parm.save('4amber_' + self.root_name + '.LES.rst7', overwrite=True)

  def write_LES_asu_pdb(self):
    """write LES asu pdb for phenix: 4phenix_{code}.LES.pdb
    """
    les_parm = pmd.load_file('4amber_{}.LES.pdb'.format(self.root_name))
    orig_pdb_parm = pmd.load_file(self.original_pdb_file)

    selection = ':1-{}'.format(self.n_asu_residues)

    final_pdb_asu_file = '4phenix_' + self.root_name + '.LES.pdb'
    final_parm = les_parm[selection]
    final_parm.symmetry = self.symmetry
    final_parm.box = self.box
    final_parm.space_group = self.space_group
    # rename WAT and update H (HOH) occupancy
    for atom in final_parm.atoms:
      if atom.atomic_number == 1:
        atom.occupancy = atom.bond_partners[0].occupancy
        atom.bfactor = atom.bond_partners[0].bfactor
    # update HOH
    for residue in final_parm.residues:
      if residue.name.startswith('WAT'):
        residue.name = 'HOH'
    # update original resnum
    for residue, residue_template in zip(final_parm.residues, orig_pdb_parm.residues):
      residue.number = residue_template.number
    write_standard_pdb(final_parm, final_pdb_asu_file)
    self._construct_hierarchy(final_pdb_asu_file)

  def _construct_hierarchy(self, filename, sort_atoms=True):
    """ construct_hierarchy and overwrite the pdb filename

    Parameters
    ----------
    filename : str
      pdb filename
    sort_atoms : bool, default True
    """
    pdb_input = iotbx.pdb.input(filename)
    with open(filename, 'w') as fh:
      pdb_hierachy = pdb_input.construct_hierarchy(sort_atoms=True)
      fh.write(self._get_pdb_header())
      fh.write(pdb_hierachy.as_pdb_string())

  def _get_pdb_header(self):
    output = StringIO()
    parm = pmd.Structure()
    parm.space_group = self.space_group
    parm.box = self.box
    parm.symmetry = self.symmetry
    parm.write_pdb(output)
    output.seek(0)
    return output.read().strip("END")

  def _check_valid_addles_input(self, fn):
    assert os.path.exists(fn), 'make sure {} exists'.format(fn)
    input = open(fn).read()
    expected_header = EXPECTED_HEADER_TEMPLATE.format(root_name=self.root_name)
    assert expected_header in input, 'addles input must have header\n\n{}\n'.format(expected_header)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
      description='build LES parm7 and rst7 files from ASU pdb.')
  parser.add_argument('pdb', help='ASU pdb')
  parser.add_argument('--prmtop', help='prmtop for ASU (optional)')
  parser.add_argument('--rst7', help='rst7 file for ASU (optional)')
  args = parser.parse_args()
  les_builder = LESBuilder(args.pdb, prmtop=args.prmtop, rst7_file=args.rst7)
  les_builder.run()
