import os
import subprocess
import pytest
from numpy.testing import assert_almost_equal as aa_eq
import libtbx.load_env
from amber_adaptbx.tests.utils import get_fn, tempfolder
import parmed as pmd
import itertools

def test_ensure_getting_the_same_pdb_with_and_without_LES_build():
  pdb_fn = get_fn('2igd/2igd.pdb')
  command_non_les = [
          'phenix.AmberPrep',
          pdb_fn,
          'LES=False'
  ]

  command_les = [
          'phenix.AmberPrep',
          pdb_fn,
          'LES=True'
  ]

  output_pdb = '4phenix_2igd.pdb'
  with tempfolder():
    subprocess.check_call(command_non_les)
    parm_non_les_build = pmd.load_file(output_pdb)
  with tempfolder():
    subprocess.check_call(command_les)
    parm_les_build = pmd.load_file(output_pdb)

  aa_eq(parm_non_les_build.coordinates, parm_les_build.coordinates, decimal=3)
  for atom_from_non_les, atom_from_les in zip(parm_non_les_build, parm_les_build):
      assert atom_from_non_les.name == atom_from_les.name
      assert atom_from_non_les.bfactor == atom_from_les.bfactor
      assert atom_from_non_les.residue.name == atom_from_les.residue.name
      assert atom_from_non_les.occupancy == atom_from_les.occupancy

@pytest.mark.parametrize('code', ['1gdu', '2igd'])
def test_writing_HOH(code):
 # not WAT
 pdb_fn = get_fn('{code}/{code}.pdb'.format(code=code))
 command_les = [
         'phenix.AmberPrep',
         pdb_fn,
         'LES=True'
 ]
 output_pdb = '4phenix_{}.pdb'.format(code)
 output_les_pdb = '4phenix_{}.LES.pdb'.format(code)
 with tempfolder():
   subprocess.check_call(command_les)
   parm = pmd.load_file(output_pdb)
   assert 'HOH' in set(residue.name for residue in parm.residues)
   assert 'WAT' not in set(residue.name for residue in parm.residues)
   parm_les = pmd.load_file(output_les_pdb)
   assert 'HOH' in set(residue.name for residue in parm_les.residues)
   assert 'WAT' not in set(residue.name for residue in parm_les.residues)


@pytest.mark.parametrize('code', ['1gdu', '2igd'])
def test_HOH_occupancy(code):
 # H in HOH must have the same occupancy with O (HOH)
 pdb_fn = get_fn('{code}/{code}.pdb'.format(code=code))
 command_les = [
         'phenix.AmberPrep',
         pdb_fn,
         'LES=True'
 ]
 output_pdb = '4phenix_{}.pdb'.format(code)
 output_les_pdb = '4phenix_{}.LES.pdb'.format(code)
 with tempfolder():
   subprocess.check_call(command_les)
   parm = pmd.load_file(output_pdb)
   parm_les = pmd.load_file(output_les_pdb)
   wat_residues_les = [residue for residue in parm_les.residues if residue.name[:3] == 'HOH']
   wat_residues = [residue for residue in parm.residues if residue.name[:3] == 'HOH']

   for residue in wat_residues_les:
       # all atoms in residue should have the same occupancy
       assert len(set(atom.occupancy for atom in residue.atoms)) == 1
   for residue in wat_residues:
       # all atoms in residue should have the same occupancy
       assert len(set(atom.occupancy for atom in residue.atoms)) == 1

   if code in ['1gdu',]:
       # having some O (HOH) with occupancy < 1.0
       occupancy_set = set(atom.occupancy for atom in itertools.chain.from_iterable(wat_residues_les))
       assert len(occupancy_set - set([0.75, 1.0, 0.74, 0.94, 0.5, 0.88, 0.73,])) > 0
