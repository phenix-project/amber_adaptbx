from __future__ import division
from libtbx import group_args
import sys, os
import itertools
import iotbx.pdb
import argparse
from scitbx.array_family import flex
import scitbx.restraints
from libtbx.utils import Sorry
#~ import boost.python
#~ ext = boost.python.import_ext("amber_adaptbx_ext")
try:
  import sander, sanderles
except:
  raise Sorry('Unable to import "sander". Check that $AMBERHOME is set correctly to the Amber directory.')
try:
  from parmed.amber.readparm import AmberParm, Rst7  #post AmberTools15
  import parmed
except ImportError:
  try:
    from chemistry.amber.readparm import AmberParm, Rst7 #up to AmberTools15
  except ImportError:
    raise ImportError("could not import parmed modules. Check path.")

from amber_adaptbx.amber_phenix_reorder import (
    initialize_order_converter, reorder_coords_phenix_to_amber,
    reorder_force_amber_to_phenix, get_indices_convert_dict_from_array
    )

master_phil_str = """
  use_amber = False
    .type = bool
    .help = Use Amber for all the gradients in refinement
    .short_caption = Enable Amber
  topology_file_name = None
    .type = path
    .help = A topology file needed by Amber. Can be generated using phenix.AmberPrep.
    .style = bold input_file
  coordinate_file_name = None
    .type = path
    .help = A coordinate file needed by Amber. Can be generated using phenix.AmberPrep.
    .style = bold input_file
  order_mapping_file_name = None
    .type = path
    .help = A pickled file that maps amber's coordinates to phenix's coordinates. If given, reuse this map.
    .style = bold input_file
  wxc_factor = .1
    .type = float
    .style = hidden
  automatic_wxc_scale = False
    .type = bool
    .style = hidden
    .help = Use the ratio of the restraints gradident norm and the Amber \
            gradient norm to set wxc_scale
  print_amber_energies = False
    .type = bool
    .help = Print details of Amber energies during refinement
  md_engine = *sander
    .type = choice
    .help = TODO: Place holder. Will remove this option later. @Dave: do not remove for now.
"""

class geometry_manager(object):
  COUNT = 0
  # all objects share the same order_converter
  order_converter = None

  def __init__(self,
        sites_cart=None,
        energy_components=None,
        gradients=None,
        number_of_restraints=0,
        gradients_factory=flex.vec3_double,
        amber_structs=None):
    self.sites_cart = sites_cart
    self.energy_components = energy_components
    self.gradients_factory = gradients_factory
    self.number_of_restraints=number_of_restraints
    self.amber_structs=amber_structs

    # order_converter is a Python dict that map amber atom order to phenix order
    # assign later

    if geometry_manager.COUNT == 0:
      # compute order_converter from original sites_cart or load from file
      initialize_order_converter(self)

    # increase COUNT to avoid recompute order_converter
    geometry_manager.COUNT += 1

    if self.energy_components is None:
      self.energy_components = flex.double([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])

  def energies_sites(self,
        crystal_symmetry,
        log=None,
        print_amber_energies=False,
        compute_gradients=False):
    # if log==None: assert 0
    #Expand sites_cart to unit cell
    sites_cart_uc=expand_coord_to_unit_cell(self.sites_cart, crystal_symmetry)

    sander_coords = reorder_coords_phenix_to_amber(sites_cart_uc, self.order_converter['p2a'])
    if self.amber_structs.is_LES:
      sanderles.set_positions(sander_coords)
      ene, frc = sanderles.energy_forces()
    else:
      sander.set_positions(sander_coords)
      ene, frc = sander.energy_forces()
    frc = reorder_force_amber_to_phenix(frc, self.order_converter['a2p'])
    if (compute_gradients) :
      gradients_uc=flex.vec3_double(flex.double(frc)) * -1
      gradients = gradients_uc[0:self.sites_cart.size()]
      # gradients = collapse_grad_to_asu(gradients_uc, crystal_symmetry)
    else :
      gradients = self.gradients_factory(
        flex.double(self.sites_cart.size() * 3,0))
    result = energies(
      compute_gradients=compute_gradients,
      gradients=gradients,
      gradients_size=None,
      gradients_factory=None,
      normalization=False)
    result.number_of_restraints = self.number_of_restraints
    result.residual_sum = ene.tot
    ptrfunc = self.amber_structs.parm.ptr
    nbond = ptrfunc('nbonh') + ptrfunc('nbona')
    nangl = ptrfunc('ntheth') + ptrfunc('ntheta')
    nmphi = ptrfunc('nphih') + ptrfunc('nphia')
    result.energy_components = [ene.tot, ene.bond, ene.angle, ene.dihedral,
                                ene.elec + ene.elec_14, ene.vdw + ene.vdw_14,
                                nbond, nangl, nmphi]
    result.finalize_target_and_gradients()
    if log==None: log=sys.stdout
    # following forces printing of Amber energies; placeholder until
    #    this can become an input keyword
    print_amber_energies=True
    if print_amber_energies:
      print >>log, "    Amber total energy: %0.2f" %(result.residual_sum)
      print >>log, "      bonds (n=%d): %0.2f" %(result.energy_components[6],
                                           result.energy_components[1])
      print >>log, "      angles (n=%d): %0.2f" %(result.energy_components[7],
                                           result.energy_components[2])
      print>>log,  "      dihedrals (n=%d): %0.2f" %(result.energy_components[8],
                                           result.energy_components[3])
      print>>log,  "      electrostatics: %0.2f" %(result.energy_components[4])
      print>>log,  "      van der Waals: %0.2f" %(result.energy_components[5])

    return result

class energies (scitbx.restraints.energies) :
  def __init__ (self, *args, **kwds) :
    scitbx.restraints.energies.__init__(self, *args, **kwds)
    self.energy_components = None
    self.amber=True
    # import code; code.interact(local=dict(globals(), **locals()))
    # sys.exit()

  def show(self):
    print "    Amber total energy: %0.2f" %(self.residual_sum)
    print "      bonds (n=%d): %0.2f" %(self.energy_components[6],
                                             self.energy_components[1])
    print "      angles (n=%d): %0.2f" %(self.energy_components[7],
                                             self.energy_components[2])
    print "      dihedrals (n=%d): %0.2f" %(self.energy_components[8],
                                             self.energy_components[3])
    print "      electrostatics: %0.2f" %(self.energy_components[4])
    print "      van der Waals: %0.2f" %(self.energy_components[5])
    return 0

  def get_grms(self):
    from math import sqrt
    gradients_1d = self.gradients.as_double()
    grms = sum(gradients_1d**2)
    grms /= gradients_1d.size()
    grms = sqrt(grms)
    return grms

  def get_gnorm(self):
    from math import sqrt
    gradients_1d = self.gradients.as_double()
    grms = sum(gradients_1d**2)
    grms = sqrt(grms)
    return grms

  # this calculates both bond and angle rms but is called
  def angle_deviations(self, sites_cart, parm, ignore_hd=False, get_deltas=False):
    return angle_rmsd(parm, sites_cart, ignore_hd, get_deltas)

  def bond_deviations(self, sites_cart, parm, ignore_hd=False, get_deltas=False):
    return bond_rmsd(parm, sites_cart, ignore_hd, get_deltas)

  def n_angle_proxies(self, parm, ignore_hd=False):
    if not ignore_hd:
      return self.energy_components[7]
    else:
      return parm.ptr('ntheta')

  def n_bond_proxies(self, parm, ignore_hd=False):
    if not ignore_hd:
      return self.energy_components[6]
    else:
      return parm.ptr('nbona')

def print_sites_cart(sites_cart):
        for atom in sites_cart:
                print("%8.3f%8.3f%8.3f"%(atom[0], atom[1], atom[2]))

def get_amber_structs (parm_file_name, rst_file_name):
        return ext.uform(parm_file_name, rst_file_name)

def check_file(s,filename):
  if filename is None:
    raise Sorry("Filename %s is None. Please set this parameter" % s)
  if not os.path.exists(filename):
    raise Sorry("Filename %s does not exist" % filename)

def check_file(s,file_name):
  if file_name is None:
    raise Sorry('Parameter %s not set. Please supply filename.' % (s))
  if not os.path.exists(file_name):
    raise Sorry("Filename %s does not exist" % file_name)

class sander_structs ():
  def __init__ (self, parm_file_name, rst_file_name, ridingH=True):
    check_file("amber.topology_file_name", parm_file_name)
    check_file("amber.coordinate_file_name", rst_file_name)
    self.md_engine = 'sander'
    self.parm = AmberParm(parm_file_name)
    self.rst = Rst7.open(rst_file_name)
    self.ridingH = ridingH
    self.is_LES = is_prmtop_LES(parm_file_name)

    self.order_converter = None
    self.order_map_file_name = None

    if self.is_LES:
      self.inp = sanderles.pme_input()
    else:
      self.inp = sander.pme_input()
    parm = parmed.load_file(parm_file_name, rst_file_name)
    # use initial_coordinates for mapping with phenix's sites_cart
    self.initial_coordinates = parm.coordinates

def is_prmtop_LES(parm_file_name):
  with open(parm_file_name) as f:
    for line in f:
      if "FLAG LES_TYPE" in line:
        return True
    return False

def expand_coord_to_unit_cell(sites_cart, crystal_symmetry):
  sites_cart_uc = flex.vec3_double()
  cell = crystal_symmetry.unit_cell()
  sg = crystal_symmetry.space_group()
  for i, op in enumerate(sg.all_ops()):
    #~ rotn = op.r().as_double()
    #~ tln = cell.orthogonalize(op.t().as_double())
    #~ # import code; code.interact(local=dict(globals(), **locals()))
    #~ # sys.exit()
    #~ sites_cart_uc.extend( (rotn * sites_cart) + tln)

    r = op.r().as_double()
    t = op.t().as_double()
    sites_frac = cell.fractionalize(sites_cart)
    sites_cart_uc.extend( cell.orthogonalize(r*sites_frac + t) )

  return sites_cart_uc

def collapse_grad_to_asu(gradients_uc, crystal_symmetry):
  cell = crystal_symmetry.unit_cell()
  sg = crystal_symmetry.space_group()
  n_symop = sg.n_smx()
  n_asu_atoms = int(gradients_uc.size() / n_symop)
  gradients = flex.vec3_double(n_asu_atoms)
  for i, op in enumerate(sg.all_ops()):
    inv_rotn = op.r().inverse().as_double()
    tln = op.t().as_double()
    start = i*n_asu_atoms
    end = (i+1)*n_asu_atoms
    g_frac = cell.fractionalize(gradients_uc[start:end])
    gradients += inv_rotn * (g_frac-t)
  gradients = gradients * (1.0/n_symop)
  return gradients

def bond_rmsd(parm, sites_cart, ignore_hd, get_deltas=False):
  from math import acos, pi, sqrt
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
  from math import acos, pi, sqrt
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
  from math import acos, pi, sqrt
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
  from math import acos, pi, sqrt
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

