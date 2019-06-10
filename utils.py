import os
import sys
import numpy as np
import itertools
import time
from math import acos, pi, sqrt
from contextlib import contextmanager
import tempfile
from shutil import rmtree
from scitbx.array_family import flex
from libtbx.utils import Sorry
from libtbx import easy_run

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

class extreme(dict):
  def __init__(self, size=10):
    self.size=size
    self[-1e9]=None

  def __repr__(self):
    keys = self.keys()
    keys.sort()
    keys.reverse()
    outl = getattr(self, 'header', '')
    for key in keys:
      item = self[key]
      for i in range(9):
        atom = getattr(item, 'atom%d' % i, None)
        if atom:
          outl += ' %4s %3s %7s %1s %1s -' % (atom.name,
                                              atom.residue.name,
                                              atom.idx,
                                              atom.residue.chain,
                                              atom.altloc,
            )
      outl = outl[:-1]
      eq = getattr(item.type, 'theteq', None)
      if eq is None:
        eq = getattr(item.type, 'req', None)
      if eq: outl += ' %7.3f' % eq
      outl += ' %7.3f' % item.model 
      outl += ' %7.3f' % key
      outl += '\n'
    return outl

  def process(self, value, item):
    if value>min(self.keys()):
      self[value]=item
    if min(self.keys())==-1e9: del self[-1e9]
    if len(self)>self.size:
      del self[min(self.keys())]

def print_cmd(cmd, verbose=False):
  print "\n~> %s\n" % cmd

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

def bond_rmsd(parm,
              sites_cart,
              ignore_hd,
              get_deltas=False,
              get_extremes=False,
              verbose=False,
  ):
  if verbose: print "starting bond_rmsd: %s" % time.strftime("%H:%M:%S")
  ignore_hd=True   # dac timing test
  if ignore_hd:
    bonds = parm.bonds_without_h
  else:
    bonds = itertools.chain(parm.bonds_inc_h, parm.bonds_without_h)
  bond_deltas = []
  bond_extremes = extreme()
  bond_extremes.header  = '  Bond deltas from Amber ideals\n'
  bond_extremes.header += '    Atoms %s ideal   model   delta\n' % (' '*36)
  # save coordinates here since calling parm.coordinates is time consumming
  parm_coordinates = parm.coordinates
  for i, bond in enumerate(bonds):
    atom1_idx= bond.atom1.idx
    atom2_idx= bond.atom2.idx
    natoms=len(sites_cart)
    # in non-P1 space groups, amber topology knows entire unit cell bonds
    # only use bonds from 1st ASU
    if atom1_idx >= natoms or atom2_idx >=natoms:
      continue
    atom1 = parm_coordinates[atom1_idx]
    atom2 = parm_coordinates[atom2_idx]
    dx = atom1[0] - atom2[0]
    dy = atom1[1] - atom2[1]
    dz = atom1[2] - atom2[2]
    delta = bond.type.req - sqrt(dx*dx + dy*dy + dz*dz)
    #print "bond deltas:  %6d %6d %6d  %7.2f" % ( i, atom1_idx, atom2_idx, delta )
    if get_extremes:
      bond.model = sqrt(dx*dx + dy*dy + dz*dz)
      bond_extremes.process(delta, bond)
    bond_deltas.append(delta)
  bond_deltas = flex.double(bond_deltas)
  b_sq  = bond_deltas * bond_deltas
  b_ave = sqrt(flex.mean_default(b_sq, 0))
  b_max = sqrt(flex.max_default(b_sq, 0))
  b_min = sqrt(flex.min_default(b_sq, 0))
  if verbose: print "done with bond_rmsd: %s" % time.strftime("%H:%M:%S")
  if not get_deltas:
    return b_min, b_max, b_ave
  else:
    return (b_min, b_max, b_ave), bond_deltas, bond_extremes

def bond_rmsZ(parm, sites_cart, ignore_hd, get_deltas=False):
  if ignore_hd:
    bonds = parm.bonds_without_h
  else:
    bonds = itertools.chain(parm.bonds_inc_h, parm.bonds_without_h)
  bond_Zs = []
  # save coordinates here since calling parm.coordinates is time consumming
  parm_coordinates = parm.coordinates
  for i, bond in enumerate(bonds):
    atom1_idx= bond.atom1.idx
    atom2_idx= bond.atom2.idx
    natoms=len(sites_cart)
    # in non-P1 space groups, amber topology knows entire unit cell bonds
    # only use bonds from 1st ASU
    if atom1_idx >= natoms or atom2_idx >=natoms:
      continue
    atom1 = parm_coordinates[atom1_idx]
    atom2 = parm_coordinates[atom2_idx]
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

def angle_rmsd(parm,
               sites_cart,
               ignore_hd,
               get_deltas=False,
               get_extremes=False,
               verbose=False,
  ):
  if verbose: print "starting angle_rmsd: %s" % time.strftime("%H:%M:%S")
  ignore_hd=True   # dac timing test
  if ignore_hd:
    angles = parm.angles_without_h
  else:
    angles = itertools.chain(parm.angles_inc_h, parm.angles_without_h)
  angle_deltas = []
  angle_extremes = extreme()
  angle_extremes.header  = '  Angle deltas from Amber ideals\n'
  angle_extremes.header += '    Atoms %s ideal   model   delta\n' % (' '*59)
  # save coordinates here since calling parm.coordinates is time consumming
  parm_coordinates = parm.coordinates
  for i, angle in enumerate(angles):
    # in non-P1 space groups, amber topology knows entire unit cell angles
    # only use angles from 1st ASU
    atom1_idx= angle.atom1.idx
    atom2_idx= angle.atom2.idx
    atom3_idx= angle.atom3.idx
    natoms=len(sites_cart)
    if atom1_idx >= natoms or atom2_idx >=natoms or atom3_idx >=natoms:
      continue
    atom1 = parm_coordinates[atom1_idx]
    atom2 = parm_coordinates[atom2_idx]
    atom3 = parm_coordinates[atom3_idx]
    a = [ atom1[0]-atom2[0], atom1[1]-atom2[1], atom1[2]-atom2[2] ]
    b = [ atom3[0]-atom2[0], atom3[1]-atom2[1], atom3[2]-atom2[2] ]
    a = flex.double(a)
    b = flex.double(b)
    acosarg = a.dot(b)/(a.norm()*b.norm())
    if acosarg >= 1.0: acosarg = 0.9999999
    if acosarg <= -1.0: acosarg = -0.9999999
    delta = angle.type.theteq - acos(acosarg)*180/pi
    assert abs(delta)<360
    if get_extremes:
      angle.model = acos(acosarg)*180/pi
      angle_extremes.process(delta, angle)
    angle_deltas.append(delta)
  angle_deltas= flex.double(angle_deltas)
  a_sq  = angle_deltas * angle_deltas
  a_ave = sqrt(flex.mean_default(a_sq, 0))
  a_max = sqrt(flex.max_default(a_sq, 0))
  a_min = sqrt(flex.min_default(a_sq, 0))
  if verbose: print "done with angle_rmsd: %s" % time.strftime("%H:%M:%S")
  if not get_deltas:
    return (a_min, a_max, a_ave)
  else:
    return (a_min, a_max, a_ave), angle_deltas, angle_extremes

def angle_rmsZ(parm, sites_cart, ignore_hd, get_deltas=False):
  if ignore_hd:
    angles = parm.angles_without_h
  else:
    angles = itertools(parm.angles_inc_h, parm.angles_without_h)
  angle_Zs = []
  # save coordinates here since calling parm.coordinates is time consumming
  parm_coordinates = parm.coordinates
  for i, angle in enumerate(angles):
    atom1_idx= angle.atom1.idx
    atom2_idx= angle.atom2.idx
    atom3_idx= angle.atom3.idx
    natoms=len(sites_cart)
    if atom1_idx >= natoms or atom2_idx >=natoms or atom3_idx >=natoms:
      continue
    atom1 = parm_coordinates[atom1_idx]
    atom2 = parm_coordinates[atom2_idx]
    atom3 = parm_coordinates[atom3_idx]
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
    # Hai: where is "t"?
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
  return
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

def write_rst7_from_GeometryManager(geom, crystal_symmetry, filename):
  # @Nigel: Please use this to write phenix pdb file to rst7 file after doing
  # phenix.refine or phenix.geometry_minimization minimise=phenix_all
  """
  Parameters
  ----------
  geom : amber_adaptbx.GeometryManager object 
  crystal_symmetry : ?
  filename : str
    output rst7 filename

  Examples
  --------
  >>> write_rst7_from_GeometryManager(geom, '4amber_minimized.rst7')
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

def build_unitcell(asu_pdb_file, output_file, use_amber_unitcell=False):
  '''build unitcell from asu pdb file

  Parameters
  ----------
  asu_pdb_file : str, ASU pdb file name
  output_file : str, unitcell pdb file name
  use_amber_unitcell : bool, default False
      if True, use amber UnitCell program to construct unitcell.
      (require asu_pdb_file to have remark 290)
      If False, use phenix stuff
  '''
  if use_amber_unitcell:
    with open(asu_pdb_file ) as fh:
      if 'REMARK 290   SMTRY' not in fh.read():
        print('UnitCell program requires pdb file to have "REMARK 290   SMTRY" line')
        raise ValueError('Should use option: use_amber_unitcell=False')
    cmd = os.path.join( os.environ["LIBTBX_BUILD"], '..', 'conda_base', 
                  'bin', 'UnitCell' )
    cmd += " -p %s -o %s" % (asu_pdb_file, output_file)
    print_cmd(cmd)
    ero = easy_run.fully_buffered(cmd)
    ero.show_stdout()
    ero.show_stderr()
  else:
    import iotbx.pdb
    import cctbx
    pdb_inp = iotbx.pdb.input(asu_pdb_file)
    cs = pdb_inp.crystal_symmetry()
    ph = pdb_inp.construct_hierarchy()
    ph_p1 = ph.expand_to_p1(crystal_symmetry=cs)
    abc = cs.unit_cell().parameters()[:6]
    cs_p1 = cctbx.crystal.symmetry(abc, "P 1")
    tmp_file = output_file + 'xxx'
    ph_p1.write_pdb_file(tmp_file, crystal_symmetry=cs_p1)
    # make column 21 blank, since parmed is pretty strict about this:
    fout = file(output_file,"wb")
    with open(tmp_file) as fin:
       for line in fin:
          if line[0:4] == "ATOM" or line[0:6] == "HETATM" \
             or line[0:6] == "ANISOU":
             fout.write( line[0:20] + " " + line[21:] )
          else:
             fout.write(line)
    os.remove(tmp_file)
       

# atom ordering
from collections import OrderedDict

def get_indices(ids_dict, big_arr):
  """

  Examples
  --------
  >>> ids_dict = {'1.2 4.6 7.9': 0, '3.5 7.9 9.1': 1}
  >>> arr = [[1.234, 4.56, 7.89],
  ...        [3.466, 7.893, 9.134]]
  >>> get_indices(ids_dict, arr)
  array([0, 1])
  """
  msg = "sites_cart (read from pdb) and amber_coords (read from rst7) do not have matched value"
  string_list = []
  for index, arr in enumerate(big_arr):
    try:
      string_list.append(ids_dict[' '.join(str(round3(i)) for i in arr)])
    except KeyError:
      print("**********ATTENTION***************")
      print(msg)
      print('atom index = {}, coordinates in rst7 file = {}'.format(index, arr))
      attempted_key = ' '.join(str(round3(i)) for i in arr)
      print('argument to sites_cart_ids dict: {}'.format(attempted_key))
      print ids_dict
      raise KeyError(msg)
  return np.array(string_list)

def round3(a):
  """

  Examples
  --------
  >>> round3(1.2350000..)
  1.235
  """
  #  format the input number in the way that would have been done
  #    when 4phenix_xxxx.pdb was made:
  y = '%8.3f' % a 

  if y == -0.000:
    # avoid confusion between '-0.000' and '0.000'
    y = 0.000
  return y

def make_dict(big_arr):
  """
  
  Examples
  --------
  >>> make_dict([[1.234, 4.56, 7.89],
  ...            [3.466, 7.893, 9.134]])
  {'1.2 4.6 7.9': 0, '3.5 7.9 9.1': 1}
  """
  return OrderedDict((' '.join(str(round3(i)) for i in arr), idx) for (idx, arr) in enumerate(big_arr))

def get_indices_convert_dict_from_array(sites_cart, amber_coords):
  """this method will be use in amber_adaptbx/__init__.py

  Parameters
  ----------
  sites_cart : flex.vec3_double
  amber_coords : 2D numpy array

  Returns
  -------
  dict(a2p=, p2a=)

  Examples
  --------
  >>> # make a fake phenix's sites_cart (flex.vec3_double)
  >>> # make a fake amber's coordinates (2D array)
  """
  # 2D array
  old_arr = np.asarray(amber_coords)
  new_arr = np.asarray(sites_cart)

  old_ids = make_dict(old_arr)
  new_ids = make_dict(new_arr)

  return {'p2a': get_indices(new_ids, old_arr),
          'a2p': get_indices(old_ids, new_arr)}
