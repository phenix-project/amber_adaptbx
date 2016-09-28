from __future__ import division
from math import sqrt
from libtbx import group_args
import sys, os
import itertools
import iotbx.pdb
import argparse
from scitbx.array_family import flex
import scitbx.restraints
from libtbx.utils import Sorry
try:
  import sander, sanderles
except:
  raise Sorry('Unable to import "sander". Check that $AMBERHOME is set correctly to the Amber directory.')

# require the most updated ParmEd version
# AmberTools >= 16
from parmed.amber.readparm import AmberParm, Rst7
import parmed

from amber_adaptbx.amber_phenix_reorder import (
    initialize_order_converter, reorder_coords_phenix_to_amber,
    reorder_force_amber_to_phenix, get_indices_convert_dict_from_array
)

from amber_adaptbx.utils import (
        tempfolder,
        expand_coord_to_unit_cell,
        get_amber_structs,
        bond_rmsd,
        bond_rmsZ,
        angle_rmsZ,
        angle_rmsd,
        is_prmtop_LES,
        collapse_grad_to_asu,
        check_file,
        print_sites_cart,
)

__all__ = ['geometry_manager', 'energies', 'sander_structs']

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
    self.amber_structs.sander_engine.set_positions(sander_coords)
    ene, frc = self.amber_structs.sander_engine.energy_forces()
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
    if print_amber_energies or 1:
      if 1:
        print >>log, """  Amber total: %0.2f bonds (n=%d): %0.2f angles (n=%d): %0.2f diheds (n=%d): %0.2f elec.: %0.2f vdW: %0.2f""" %(
          result.residual_sum,
          result.energy_components[6],
          result.energy_components[1],
          result.energy_components[7],
          result.energy_components[2],
          result.energy_components[8],
          result.energy_components[3],
          result.energy_components[4],
          result.energy_components[5])
      else:
        print >>log, """  Amber total energy: %0.2f
    bonds (n=%d): %0.2f
    angles (n=%d): %0.2f
    dihedrals (n=%d): %0.2f
    electrostatics: %0.2f
    van der Waals: %0.2f""" %(result.residual_sum,
                              result.energy_components[6],
                              result.energy_components[1],
                              result.energy_components[7],
                              result.energy_components[2],
                              result.energy_components[8],
                              result.energy_components[3],
                              result.energy_components[4],
                              result.energy_components[5])
    #  assert 0
    #else:
    #  from libtbx.introspection import show_stack
    #  show_stack()
    #  assert 0

    return result

class energies(scitbx.restraints.energies) :
  def __init__ (self, *args, **kwds) :
    scitbx.restraints.energies.__init__(self, *args, **kwds)
    self.energy_components = None
    self.amber=True
    # import code; code.interact(local=dict(globals(), **locals()))
    # sys.exit()

  def show(self):
    assert 0 # not writing to log...
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
    gradients_1d = self.gradients.as_double()
    grms = sum(gradients_1d**2)
    grms /= gradients_1d.size()
    grms = sqrt(grms)
    return grms

  def get_gnorm(self):
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

class sander_structs(object):
  def __init__ (self, parm_file_name, rst_file_name, ridingH=True):
    check_file("amber.topology_file_name", parm_file_name)
    check_file("amber.coordinate_file_name", rst_file_name)
    self.md_engine = 'sander'
    self.parm = parmed.load_file(parm_file_name, xyz=rst_file_name)
    # where do we need this self.rst?
    self.rst = Rst7.open(rst_file_name)
    self.ridingH = ridingH
    self.is_LES = is_prmtop_LES(parm_file_name)
    self.sander_engine = sanderles if self.is_LES else sander

    self.order_converter = None
    self.order_map_file_name = None
    self.inp = self.sander_engine.pme_input()
