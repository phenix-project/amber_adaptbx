import os
import sys
import parmed as pmd
from amber_adaptbx.make_addles_input import addles_input
from libtbx import easy_run
from amber_adaptbx.scripts import reduce_to_les
import argparse

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
the *_fix.pdb output does not have REMARK 290 yet. I (Hai) introdued to ParmEd and waiting for Jason
to merge back to amber git repo.

"""


class LESBuilder(object):
  """Class to build LES parm, rst7 and pdb file.

  Parameters
  ----------
  original_pdb_file : str, ASU pdb (e.g downloaded from rcsb)
  prmtop : str, prmtop of the single conformation unitcell pdb, default None (optional)
      if not provided, LESBuilder will use `build_parm_single_conformation` to build prmtop.
      It's better to specify prmtop since `build_parm_single_conformation` does not handle ligand.
      DAC note: should we not make the single conformer versions mandatory?
      You can use AmberPrep to prepare this prmtop file
  rst7_file : str, corresponding rst7 filename, default None (optional)
     only needed if `prmtop` is given

  Examples
  --------
  >>> # build LES for pdb that does not have ligands
  >>> builder = LESBuilder('./2igd.pdb')
  >>> builder.build()

  >>> # build LES for pdb that have ligands
  >>> # use phenix.AmberPrep to build 2igda.parm7, 2igda.rst7
  >>> builder = LESBuilder('./2igd.pdb', prmtop='./2igda.parm7', rst7_file='./2igda.rst7')
  >>> builder.build()
  """

  def __init__(self, original_pdb_file, prmtop=None, rst7_file=None):
    self.original_pdb_file = original_pdb_file
    self.new_pdb_with_H = None

    # prmtop and rst7_file were built from phenix.AmberPrep.
    self.prmtop = prmtop
    self.rst7_file = rst7_file

    self.root_name = '.'.join(os.path.basename(
        self.original_pdb_file).split('.')[:-1])
    self.uc_parm_file = None
    self.rst7_fn_fix = None
    asu_pdb_parm = pmd.load_file(original_pdb_file)
    self.space_group = asu_pdb_parm.space_group
    self.symmetry = asu_pdb_parm.symmetry
    self.box = asu_pdb_parm.box
    self.n_asu_residues = len(asu_pdb_parm.residues)

  def run(self):
    # main driver
    # use UnitCell to build unitcell from asu pdb
    self.build_unitcell()
    # use reduce to add hydrogens to unitcell pdb
    self.add_hydrogens()
    # use tleap to build parm7 and rst7 files for single conformation
    self.build_parm_single_conformation()
    # use addles to construct LES parm7 and rst7 files
    self.build_LES_parm()
    # fix previous step: need to update rst7 coordinates since addles use the same coordinates
    # for alternative atoms. Also correctly label atom and residue names.
    self.update_LES_coordinates_from_uc()
    # generate final ASU pdb for phenix (with unitcell, symmetry)
    self.write_asu_LES()
    # rename files
    self.rename_4amber_4phenix()

  def build_unitcell(self):
    new_pdb = self.root_name + '_uc.pdb'
    command_build_unitcell = 'UnitCell -p {} -o {}'.format(
        self.original_pdb_file, new_pdb)
    # print(command_build_unitcell)
    easy_run.fully_buffered(command_build_unitcell)

  def add_hydrogens(self):
    new_pdb = self.root_name + '_uc.pdb'
    new_pdb_with_H = self.root_name + '_uc_H.pdb'
    command_add_hydrogens = 'reduce -build -nuclear {} > {} 2>reduce.log'.format(
        new_pdb, new_pdb_with_H)
    # print(command_add_hydrogens)
    easy_run.fully_buffered(command_add_hydrogens)
    self.new_pdb_with_H = new_pdb_with_H

  def build_parm_single_conformation(self):
    # Note: should require that next if statement fails(?)
    if self.prmtop is None and self.rst7_file is None:
      # only do this step if prmtop and rst7_file is not provided
      # assert 0
      leap_input = """
      logFile leap.log
      source leaprc.protein.ff14SB
      source leaprc.DNA.OL15
      source leaprc.RNA.OL3
      source leaprc.water.tip3p
      source leaprc.water.spce
      source leaprc.gaff2
      x = loadPdb {uc_pdb}
      set x box {{20.000   20.000   20.000}}
      set default nocenter on
      set default reorder_residues off
      saveAmberParm x {root_name}a.parm7 {root_name}a.rst7
      quit
      """.format(uc_pdb=self.new_pdb_with_H, root_name=self.root_name)

      leap_input = '\n'.join(line.strip()
                             for line in leap_input.split('\n'))

      with open('leap.in', 'w') as fh:
        fh.write(leap_input)
      print('building single conformation prmtop and rst7 files')
      easy_run.fully_buffered('tleap -f leap.in')

      # update unitcell info
      output_parm = '{}a.parm7'.format(self.root_name)
      uc_parm = pmd.load_file(self.original_pdb_file)
      parm = pmd.load_file(output_parm)
      parm.box = uc_parm.box
      parm.save(output_parm, overwrite=True)
      self.prmtop = output_parm
      self.rst7_file = '{}a.rst7'.format(self.root_name)
    else:
      print("Using given {} and {} files".format(self.prmtop, self.rst7_file))

  def build_LES_parm(self):
    # create addles.in
    print "\n============================================================"
    print " Building the prmtop and rst7 file with alternate conformers"
    print "============================================================"
    commands = addles_input(self.original_pdb_file,
                            self.prmtop,
                            self.rst7_file)
    fn = 'addles.in'
    with open(fn, 'w') as fh:
      fh.write('\n'.join(commands))
    command_addles = 'addles <{}'.format(fn)
    # print(command_addles)
    easy_run.fully_buffered(command_addles)

  def update_LES_coordinates_from_uc(self):
    """label_alternates + update LES coordinates
    """
    uc_parm = pmd.load_file(self.new_pdb_with_H)
    parm = pmd.load_file(self.root_name + 'ab.parm7', self.root_name + 'ab.rst7')
    # add space group
    parm.space_group = self.space_group
    parm.symmetry = self.symmetry
    parm.box = self.box
    self.rst7_fn_fix = '4amber_' + self.root_name + '.LES.rst7'
    reduce_to_les.update_rst7_and_pdb_coordinates_LES(template_parm=uc_parm,
                                                      target_parm=parm,
                                                      rst7_fn=self.rst7_fn_fix)
    parm.write_pdb('4amber_{}.pdb'.format(self.root_name), standard_resnames=True)
    parm.save(self.rst7_fn_fix, overwrite=True)

  def write_asu_LES(self):
    """write LES ASU pdb for phenix: 4phenix_{code}.LES.pdb
    """
    # print('write_asu_LES')
    les_parm = pmd.load_file('4amber_{}.pdb'.format(self.root_name))

    selection = ':1-{}'.format(self.n_asu_residues)

    final_pdb_asu_file = '4phenix_' + self.root_name + '.LES.pdb'
    final_parm = les_parm[selection]
    final_parm.symmetry = self.symmetry
    final_parm.box = self.box
    final_parm.space_group = self.space_group
    final_parm.write_pdb(final_pdb_asu_file, standard_resnames=True)

  def rename_4amber_4phenix(self):
    """rename
    """
    # print('rename files')
    original_name = self.root_name + 'ab.parm7'
    final_parm_name = '4amber_' + self.root_name + '.LES.prmtop'
    os.rename(original_name, final_parm_name)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
      description='build LES parm7 and rst7 files from ASU pdb. Note: Syntax might be changed')
  parser.add_argument('pdb', help='ASU pdb')
  parser.add_argument('--prmtop', help='prmtop for ASU (optional)')
  parser.add_argument('--rst7', help='rst7 file for ASU (optional)')
  args = parser.parse_args()
  les_builder = LESBuilder(args.pdb, prmtop=args.prmtop, rst7_file=args.rst7)
  les_builder.run()
