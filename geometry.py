from __future__ import division
import sys
from scitbx.array_family import flex

from amber_adaptbx.les_builder.amber_phenix_reorder import (
    initialize_order_converter, reorder_coords_phenix_to_amber,
    reorder_force_amber_to_phenix
)

from amber_adaptbx.utils import (
    expand_coord_to_unit_cell,
)

from amber_adaptbx.energy import energies

__all__ = ['geometry_manager']


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
    self.number_of_restraints = number_of_restraints
    self.amber_structs = amber_structs

    # order_converter is a Python dict that map amber atom order to phenix order
    # assign later

    if geometry_manager.COUNT == 0:
      # compute order_converter from original sites_cart or load from file
      initialize_order_converter(self)

    # increase COUNT to avoid recompute order_converter
    geometry_manager.COUNT += 1

    if self.energy_components is None:
      self.energy_components = flex.double([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

  def energies_sites(self,
                     crystal_symmetry,
                     log=None,
                     print_amber_energies=False,
                     compute_gradients=False):
    # if log is None: assert 0
    # Expand sites_cart to unit cell
    sites_cart_uc = expand_coord_to_unit_cell(self.sites_cart, crystal_symmetry)

    sander_coords = reorder_coords_phenix_to_amber(sites_cart_uc, self.order_converter['p2a'])
    self.amber_structs.sander_engine.set_positions(sander_coords)
    ene, frc = self.amber_structs.sander_engine.energy_forces()
    frc = reorder_force_amber_to_phenix(frc, self.order_converter['a2p'])
    if (compute_gradients):
      gradients_uc = flex.vec3_double(flex.double(frc)) * -1
      gradients = gradients_uc[0:self.sites_cart.size()]
    else:
      gradients = self.gradients_factory(
          flex.double(self.sites_cart.size() * 3, 0))
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
    if log is None:
      log = sys.stdout
    # following forces printing of Amber energies; placeholder until
    #    this can become an input keyword
    if print_amber_energies or 1:
      print >>log, """  Amber total: %0.2f bonds (n=%d): %0.2f angles (n=%d): %0.2f diheds (n=%d): %0.2f elec.: %0.2f vdW: %0.2f""" % (
          result.residual_sum,
          result.energy_components[6],
          result.energy_components[1],
          result.energy_components[7],
          result.energy_components[2],
          result.energy_components[8],
          result.energy_components[3],
          result.energy_components[4],
          result.energy_components[5])

    return result
