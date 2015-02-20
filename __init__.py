from __future__ import division
from libtbx import group_args
import sys, os
import iotbx.pdb
import argparse
from scitbx.array_family import flex
import scitbx.restraints
import boost.python
ext = boost.python.import_ext("amber_adaptbx_ext")
import sander
from chemistry.amber.readparm import AmberParm, Rst7


master_phil_str = """
  use_amber = False
    .type = bool
  topology_file_name = None
    .type = path
  coordinate_file_name = None
    .type = path
  wxc_factor = None
    .type = float
  use_sander = True 
    .type = bool
"""

class geometry_manager(object):

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

    if self.energy_components is None:
      self.energy_components = flex.double([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
  def energies_sites(self,
        crystal_symmetry,
        compute_gradients=False):
    # import code; code.interact(local=dict(globals(), **locals()))
    #Expand sites_cart to unit cell
    sites_cart_uc=expand_coord_to_unit_cell(self.sites_cart, crystal_symmetry)

    if hasattr(self.amber_structs,'parm'):
      # print "\n\nUSING SANDER\n\n"
      sander_coords = list(sites_cart_uc.as_double())
      #~ import code; code.interact(local=dict(globals(), **locals()))
      #~ sys.exit()
      sander.set_positions(sander_coords)
      ene, frc = sander.energy_forces()
      # sander.cleanup()
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

    else:
      # print "\n\nUSING MDGX\n\n"
      #Convert flex arrays to C arrays
      sites_cart_c=ext.ExtractVec(sites_cart_uc.as_double())
      gradients_c=ext.ExtractVec(flex.double(sites_cart_uc.size() * 3, 0))
      energy_components_c=ext.ExtractVec(self.energy_components)

      # Call c++ interface to call mdgx to calculate new gradients and target
      ext.callMdgx(sites_cart_c, gradients_c, energy_components_c, self.amber_structs)
      if (compute_gradients) :
        # import code; code.interact(local=dict(globals(), **locals()))
        # sys.exit()
        gradients_uc = self.gradients_factory(gradients_c) * -1
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
      result.residual_sum = float(energy_components_c[0])
      result.energy_components = list(energy_components_c)
      result.finalize_target_and_gradients()
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

class sander_structs ():
  def __init__ (self, parm_file_name, rst_file_name, ridingH=True):
    self.parm = AmberParm(parm_file_name)
    self.rst = Rst7.open(rst_file_name)
    self.inp = sander.pme_input()
    self.ridingH = ridingH

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
    bonds = parm.bonds_inc_h + parm.bonds_without_h
  bond_deltas = []
  for i, bond in enumerate(bonds):
    atom1= bond.atom1.starting_index
    atom2= bond.atom2.starting_index
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
    delta = bond.bond_type.req - sqrt(dx*dx + dy*dy + dz*dz)
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
    bonds = parm.bonds_inc_h + parm.bonds_without_h
  bond_Zs = []
  for i, bond in enumerate(bonds):
    atom1= bond.atom1.starting_index
    atom2= bond.atom2.starting_index
    natoms=len(sites_cart)
    if atom1 >= natoms or atom2 >=natoms:
      continue
    atom1 = sites_cart[atom1]
    atom2 = sites_cart[atom2]
    dx = atom1[0] - atom2[0]
    dy = atom1[1] - atom2[1]
    dz = atom1[2] - atom2[2]
    Z = sqrt(bond.bond_type.k)*(bond.bond_type.req - sqrt(dx*dx + dy*dy + dz*dz))
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
    angles = parm.angles_inc_h + parm.angles_without_h
  angle_deltas = []
  for i, angle in enumerate(angles):
    # in non-P1 space groups, amber topology knows entire unit cell angles
    # only use angles from 1st ASU
    atom1= angle.atom1.starting_index
    atom2= angle.atom2.starting_index
    atom3= angle.atom3.starting_index
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
    delta = angle.angle_type.theteq - acos(a.dot(b)/(a.norm()*b.norm()))
    delta *= 180/pi
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
    angles = parm.angles_inc_h + parm.angles_without_h
  angle_Zs = []
  for i, angle in enumerate(angles):
    atom1= angle.atom1.starting_index
    atom2= angle.atom2.starting_index
    atom3= angle.atom3.starting_index
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
    Z = sqrt(angle.angle_type.k)*(angle.angle_type.theteq - acos(a.dot(b)/(a.norm()*b.norm())))
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


def run(pdb,prmtop, crd):

  #===================================================================#
  #                                                                   #
  #  BEFORE C++                                                       #
  #                                                                   #
  #===================================================================#

  #file i/o
  pdb_file = os.path.abspath(pdb)
  pdb_inp = iotbx.pdb.input(file_name=pdb_file)
  pdb_atoms = pdb_inp.atoms_with_labels()
  symm = pdb_inp.crystal_symmetry()
  xray_structure = pdb_inp.xray_structure_simple(enable_scattering_type_unknown=True)


  #     initiate flex arrays for coordinates, gradients, energy
  sites_cart=xray_structure.sites_cart()
  gradients=flex.double(len(sites_cart)*3)
  target=flex.double([6.7,1.0,2.0,3.0,4.0,5.0,0.0,0.0,0.0,0.0])
  print "Number of atom sites: %d " %sites_cart.size()
  print "\nGradients and target BEFORE C call:"
  print list(gradients[1:10])
  print target[0]

  #===================================================================#
  #                                                                   #
  #  CALL C++                                                         #
  #                                                                   #
  #===================================================================#

  U=ext.uform(prmtop, crd)

  #Convert flex arrays to C arrays
  sites_cart_c=ext.ExtractVec(sites_cart.as_double())
  gradients_c=ext.ExtractVec(gradients)
  target_c=ext.ExtractVec(target)

  # Call c++ interface to call mdgx to calculate new gradients and target
  ext.callMdgx(sites_cart_c, gradients_c, target_c, U)

  # Convert back into python types (eg. into flex arrays for phenix to use)
  gradients=flex.vec3_double(gradients_c)*-1

  target= flex.double(target_c)

  #===================================================================#
  #                                                                   #
  #  AFTER C++                                                        #
  #                                                                   #
  #===================================================================#


  print "\nGradients and target AFTER C call:"
  print list(gradients[0:10])
  print target[0]
  print target[9]

  print "Amber_total_energy: %7.6f"             %(target[0])
  print "  bonds (n= ): %7.6f"                  %(target[1])
  print "  angles (n= ): %7.6f"                         %(target[2])
  print "  dihedrals (n= ): %7.6f"              %(target[3])
  print "  electrostatics: %7.6f"               %(target[4])
  print "  vanderWaals: %7.6f"                  %(target[5])

  return 0

if __name__ == "__main__" :
        parser = argparse.ArgumentParser()
        parser.add_argument("pdb", help="name of pdb file")
        parser.add_argument("prmtop", help="name of topology file")
        parser.add_argument("crd", help="name of coordinate file")
        args = parser.parse_args()
        run(args.pdb,args.prmtop, args.crd)
