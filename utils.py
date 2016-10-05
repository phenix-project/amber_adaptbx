import os
import numpy as np
import itertools
from math import acos, pi, sqrt
from contextlib import contextmanager
import tempfile
from shutil import rmtree
from scitbx.array_family import flex
from libtbx.utils import Sorry

__all__ = [
        'tempfolder',
        'expand_coord_to_unit_cell',
        'get_amber_structs',
        'bond_rmsd',
        'bond_rmsZ',
        'angle_rmsZ',
        'angle_rmsd',
        'is_prmtop_LES',
        'collapse_grad_to_asu',
        'check_file',
        'print_sites_cart',
        'build_unitcell',
]


@contextmanager
def tempfolder():
    my_temp = tempfile.mkdtemp()
    cwd = os.getcwd()
    os.chdir(my_temp)
    yield
    os.chdir(cwd)
    rmtree(my_temp)

def expand_coord_to_unit_cell(sites_cart, crystal_symmetry):
  sites_cart_uc = flex.vec3_double()
  cell = crystal_symmetry.unit_cell()
  sg = crystal_symmetry.space_group()
  for i, op in enumerate(sg.all_ops()):
    r = op.r().as_double()
    t = op.t().as_double()
    sites_frac = cell.fractionalize(sites_cart)
    sites_cart_uc.extend( cell.orthogonalize(r*sites_frac + t) )

  return sites_cart_uc

def bond_rmsd(parm, sites_cart, ignore_hd, get_deltas=False):
  if ignore_hd:
    bonds = parm.bonds_without_h
  else:
    bonds = itertools.chain(parm.bonds_inc_h, parm.bonds_without_h)
  bond_deltas = []
  for i, bond in enumerate(bonds):
    atom1= bond.atom1.idx
    atom2= bond.atom2.idx
    natoms=len(sites_cart)
    # in non-P1 space groups, amber topology knows entire unit cell bonds
    # only use bonds from 1st ASU
    if atom1 >= natoms or atom2 >=natoms:
      continue
    atom1 = sites_cart[atom1]
    atom2 = sites_cart[atom2]
    dx = atom1[0] - atom2[0]
    dy = atom1[1] - atom2[1]
    dz = atom1[2] - atom2[2]
    delta = bond.type.req - sqrt(dx*dx + dy*dy + dz*dz)
    bond_deltas.append(delta)
  bond_deltas = flex.double(bond_deltas)
  b_sq  = bond_deltas * bond_deltas
  b_ave = sqrt(flex.mean_default(b_sq, 0))
  b_max = sqrt(flex.max_default(b_sq, 0))
  b_min = sqrt(flex.min_default(b_sq, 0))
  if not get_deltas:
    return b_min, b_max, b_ave
  else:
    return (b_min, b_max, b_ave), bond_deltas

def bond_rmsZ(parm, sites_cart, ignore_hd, get_deltas=False):
  if ignore_hd:
    bonds = parm.bonds_without_h
  else:
    bonds = itertools.chain(parm.bonds_inc_h, parm.bonds_without_h)
  bond_Zs = []
  for i, bond in enumerate(bonds):
    atom1= bond.atom1.idx
    atom2= bond.atom2.idx
    natoms=len(sites_cart)
    if atom1 >= natoms or atom2 >=natoms:
      continue
    atom1 = sites_cart[atom1]
    atom2 = sites_cart[atom2]
    dx = atom1[0] - atom2[0]
    dy = atom1[1] - atom2[1]
    dz = atom1[2] - atom2[2]
    Z = sqrt(bond.type.k)*(bond.type.req - sqrt(dx*dx + dy*dy + dz*dz))
    bond_Zs.append(Z)
  bond_Zs = flex.double(bond_Zs)
  b_sq  = bond_Zs * bond_Zs
  b_ave = sqrt(flex.mean_default(b_sq, 0))
  b_max = sqrt(flex.max_default(b_sq, 0))
  b_min = sqrt(flex.min_default(b_sq, 0))
  if not get_deltas:
    return b_min, b_max, b_ave
  else:
    return (b_min, b_max, b_ave), bond_Zs

def angle_rmsd(parm, sites_cart, ignore_hd, get_deltas=False):
  if ignore_hd:
    angles = parm.angles_without_h
  else:
    angles = itertools.chain(parm.angles_inc_h, parm.angles_without_h)
  angle_deltas = []
  for i, angle in enumerate(angles):
    # in non-P1 space groups, amber topology knows entire unit cell angles
    # only use angles from 1st ASU
    atom1= angle.atom1.idx
    atom2= angle.atom2.idx
    atom3= angle.atom3.idx
    natoms=len(sites_cart)
    if atom1 >= natoms or atom2 >=natoms or atom3 >=natoms:
      continue
    atom1 = sites_cart[atom1]
    atom2 = sites_cart[atom2]
    atom3 = sites_cart[atom3]
    a = [ atom1[0]-atom2[0], atom1[1]-atom2[1], atom1[2]-atom2[2] ]
    b = [ atom3[0]-atom2[0], atom3[1]-atom2[1], atom3[2]-atom2[2] ]
    a = flex.double(a)
    b = flex.double(b)
    delta = angle.type.theteq - acos(a.dot(b)/(a.norm()*b.norm()))*180/pi
    assert abs(delta)<360
    angle_deltas.append(delta)
  angle_deltas= flex.double(angle_deltas)
  a_sq  = angle_deltas * angle_deltas
  a_ave = sqrt(flex.mean_default(a_sq, 0))
  a_max = sqrt(flex.max_default(a_sq, 0))
  a_min = sqrt(flex.min_default(a_sq, 0))
  if not get_deltas:
    return (a_min, a_max, a_ave)
  else:
    return (a_min, a_max, a_ave), angle_deltas

def angle_rmsZ(parm, sites_cart, ignore_hd, get_deltas=False):
  if ignore_hd:
    angles = parm.angles_without_h
  else:
    angles = itertools(parm.angles_inc_h, parm.angles_without_h)
  angle_Zs = []
  for i, angle in enumerate(angles):
    atom1= angle.atom1.idx
    atom2= angle.atom2.idx
    atom3= angle.atom3.idx
    natoms=len(sites_cart)
    if atom1 >= natoms or atom2 >=natoms or atom3 >=natoms:
      continue
    atom1 = sites_cart[atom1]
    atom2 = sites_cart[atom2]
    atom3 = sites_cart[atom3]
    a = [ atom1[0]-atom2[0], atom1[1]-atom2[1], atom1[2]-atom2[2] ]
    b = [ atom3[0]-atom2[0], atom3[1]-atom2[1], atom3[2]-atom2[2] ]
    a = flex.double(a)
    b = flex.double(b)
    Z = sqrt(angle.type.k)*(angle.type.theteq - acos(a.dot(b)/(a.norm()*b.norm()))*180/pi)
    angle_Zs.append(Z)
  angle_Zs= flex.double(angle_Zs)
  a_sq  = angle_Zs * angle_Zs
  a_ave = sqrt(flex.mean_default(a_sq, 0))
  a_max = sqrt(flex.max_default(a_sq, 0))
  a_min = sqrt(flex.min_default(a_sq, 0))
  if not get_deltas:
    return (a_min, a_max, a_ave)
  else:
    return (a_min, a_max, a_ave), angle_Zs

def collapse_grad_to_asu(gradients_uc, crystal_symmetry):
  # does this code really work?
  cell = crystal_symmetry.unit_cell()
  sg = crystal_symmetry.space_group()
  n_symop = sg.n_smx()
  n_asu_atoms = int(gradients_uc.size() / n_symop)
  gradients = flex.vec3_double(n_asu_atoms)
  for i, op in enumerate(sg.all_ops()):
    inv_rotn = op.r().inverse().as_double()
    start = i*n_asu_atoms
    end = (i+1)*n_asu_atoms
    g_frac = cell.fractionalize(gradients_uc[start:end])
    gradients += inv_rotn * (g_frac-t)
  gradients = gradients * (1.0/n_symop)
  return gradients

def is_prmtop_LES(parm_file_name):
  with open(parm_file_name) as f:
    for line in f:
      if "FLAG LES_TYPE" in line:
        return True
    return False

def check_file(s, file_name):
  if not os.path.exists(file_name):
    raise Sorry("Filename %s does not exist" % file_name)

def get_amber_structs(parm_file_name, rst_file_name):
  # I do not know where this will be used. Is it old code for mdgx?
  import boost.python
  ext = boost.python.import_ext("amber_adaptbx_ext")
  return ext.uform(parm_file_name, rst_file_name)

def print_sites_cart(sites_cart):
  for atom in sites_cart:
     print("%8.3f%8.3f%8.3f"%(atom[0], atom[1], atom[2]))

def write_rst7_from_geometry_manager(geom, crystal_symmetry, filename):
  # @Nigel: Please use this to write phenix pdb file to rst7 file after doing
  # phenix.refine or phenix.geometry_minimization minimise=phenix_all
  """
  Parameters
  ----------
  geom : amber_adaptbx.geometry_manager object 
  crystal_symmetry : ?
  filename : str
    output rst7 filename

  Examples
  --------
  >>> write_rst7_from_geometry_manager(geom, '4amber_minimized.rst7')
  """
  parm = geom.amber_structs.parm
  old_coords = parm.coordinates
  sites_cart_uc = expand_coord_to_unit_cell(geom.sites_cart, crystal_symmetry)
  sites_cart_uc_arr = np.asarray(sites_cart_uc, dtype='f8')
  indices_p2a = geom.order_converter['p2a']
  # reorder
  parm.coordinates = sites_cart_uc_arr[indices_p2a]
  # filename should have .rst7 ext if not specifying format
  parm.save(filename, format='restrt')
  # restore
  parm.coordinates = old_coords

def build_unitcell(asu_pdb_file, output_file):
  '''build unitcell from asu pdb file

  Parameters
  ----------
  asu_pdb_file : str, ASU pdb file name
  output_file : str, unitcell pdb file name
  '''
  import parmed as pmd
  from iotbx import pdb

  asu_pdb = pdb.input(asu_pdb_file)
  symmetry = asu_pdb.crystal_symmetry()
  asu_site_carts = asu_pdb.atoms().extract_xyz()

  # compute unitcell's coordinates from asu's coordinates
  uc_site_carts = expand_coord_to_unit_cell(asu_site_carts, symmetry)
  space_group = asu_pdb.crystal_symmetry().space_group()

  # number of ASU in unitcell
  n_asu_pdb = len(space_group.all_ops())

  asu_structure = pmd.load_file(asu_pdb_file)
  asu_n_residues = len(asu_structure.residues)
  # build residue template from parmed's Structure
  # Note: all ASUs have the same coordinates, we will update later
  uc_structure = asu_structure * n_asu_pdb

  # add TER to seperate each asu
  for index in range(n_asu_pdb):
     res_ter = uc_structure.residues[asu_n_residues*(index+1) - 1]
     res_ter.ter = True

  # TODO: write to StringIO
  # Does iotbx.pdb can read file handler?
  with tempfolder():
    fn = 'tmp.pdb'
    uc_structure.save(fn, overwrite=True)
    uc_pdb = pdb.input(fn)
    # update expanded coordinates
    uc_pdb.atoms().set_xyz(uc_site_carts)
  uc_pdb.write_pdb_file(output_file)
