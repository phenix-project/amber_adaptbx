from __future__ import print_function
from io import StringIO
import sys

from libtbx.utils import Sorry
from scitbx.array_family import flex
from libtbx import adopt_init_args
from libtbx.str_utils import make_header

import amber_adaptbx
from amber_adaptbx.les_builder.amber_phenix_reorder import (
  initialize_order_converter, reorder_coords_phenix_to_amber,
  reorder_force_amber_to_phenix
)

from amber_adaptbx.utils import expand_coord_to_unit_cell

from amber_adaptbx.energy import energies

from cctbx.geometry_restraints.manager import manager as standard_manager

class manager(standard_manager):
  COUNT = 0
  # all objects share the same order_converter
  order_converter = None
  n_sites_cart = None

  def __init__(self,
               params,
               energy_components=None,
               gradients_factory=flex.vec3_double,
               log=StringIO()):
    self.gradients_factory = gradients_factory
    adopt_init_args(self, locals(), exclude=["log"])

    self.initialisation(params, log)

    from amber_adaptbx import interface
    self.amber_structs, self.sander = interface.get_amber_struct_object(params)
    compute_gradients=False
    if manager.COUNT == 0:
      # compute order_converter from original sites_cart or load from file
      initialize_order_converter(self)
    self.print_amber_energies = params.amber.print_amber_energies
    # self.qmmask = params.amber.qmmask
    # sites_cart = self.pdb_hierarchy.atoms().extract_xyz()
    # print('sites_cart',len(sites_cart))
    # geometry = self.energies_sites(
    #   sites_cart = self.pdb_hierarchy.atoms().extract_xyz(),
    #   compute_gradients = compute_gradients,
    #   log=log,
    #   print_amber_energies=self.print_amber_energies,
    #   #qmmask=self.qmmask,
    #   )

    # increase COUNT to avoid recompute order_converter
    manager.COUNT += 1

    if self.energy_components is None:
      self.energy_components = flex.double([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

  def __repr__(self):
    return 'Amber manager'

  def initialisation(self, params, log=None):
    make_header("Initializing Amber", out=log)
    error = '''

    no filename for %s provided

    use

      %s=<filename>.%s

    '''
    print("  topology    : %s" % params.amber.topology_file_name, file=log)
    if not params.amber.topology_file_name:
      raise Sorry(error % ('topology', 'amber.topology_file_name', 'prmtop'))
    if params.amber.topology_file_name.endswith('rst7'):
      raise Sorry('possible wrong format - need .prmtop file')
    print("  atom order  : %s" % params.amber.order_file_name, file=log)
    if not params.amber.order_file_name:
      raise Sorry(error % ('order', 'amber.order_file_name', 'order'))
    if params.amber.coordinate_file_name or 1:
      print("  coordinates : %s" % params.amber.coordinate_file_name, file=log)
      if not params.amber.coordinate_file_name:
        raise Sorry(error % ('coordinate', 'amber.coordinate_file_name', 'rst7'))
    make_header('...', out=log)

  def energies_sites(self,
                     sites_cart,
                     flags=None,
                     custom_nonbonded_function=None,
                     compute_gradients=False,
                     gradients=None,
                     disable_asu_cache=False,
                     normalization=False,
                     external_energy_function=None,
                     extension_objects=[],
                     site_labels=None,
                     log=None):
    result = standard_manager.energies_sites(
      self,
      sites_cart,
      flags=flags,
      custom_nonbonded_function=custom_nonbonded_function,
      compute_gradients=False, #compute_gradients,
      gradients=gradients,
      disable_asu_cache=disable_asu_cache,
      normalization=normalization,
      external_energy_function=external_energy_function,
      extension_objects=extension_objects,
      site_labels=site_labels,
      )
    # Expand sites_cart to unit cell
    sites_cart_uc = expand_coord_to_unit_cell(
      sites_cart,
      self.crystal_symmetry)
    sander_coords = reorder_coords_phenix_to_amber(sites_cart_uc,
                                                   self.order_converter['p2a'])
    self.amber_structs.sander_engine.set_positions(sander_coords)

    ene, frc = self.amber_structs.sander_engine.energy_forces()
    if self.amber_structs.writer is not None:
      self.amber_structs.writer.add_coordinates(sander_coords)
      self.amber_structs.writer.add_box(self.amber_structs.parm.box)
    frc = reorder_force_amber_to_phenix(frc, self.order_converter['a2p'])
    if (compute_gradients):
      gradients_uc = flex.vec3_double(flex.double(frc)) * -1
      gradients = gradients_uc[0:sites_cart.size()]
    else:
      # WHY????
      gradients = self.gradients_factory(
          flex.double(sites_cart.size() * 3, 0))
    result.gradients=gradients
    result.residual_sum = ene.tot
    ptrfunc = self.amber_structs.parm.ptr
    nbond = ptrfunc('nbonh') + ptrfunc('nbona')
    nangl = ptrfunc('ntheth') + ptrfunc('ntheta')
    nmphi = ptrfunc('nphih') + ptrfunc('nphia')
    result.energy_components = [ene.tot, ene.bond, ene.angle, ene.dihedral,
                                ene.elec + ene.elec_14 + ene.gb, ene.vdw + ene.vdw_14,
                                nbond, nangl, nmphi]
    result.finalize_target_and_gradients()
    #
    # to usurp a test in statistics.py
    #
    result.bond_residual_sum=result.target
    if log is None:
      log = sys.stdout
    if self.print_amber_energies:
      import decimal
      outl = []
      def _outl(f):
        outl.append('%0.1f' % f)
      for i in range(9):
        if i: _outl(result.energy_components[i])
        else: _outl(result.residual_sum)
      energies = 'Amber: %10s %10s %10s %10s %10s %10s' % (outl[0],
                        outl[1], outl[2], outl[3], outl[4], outl[5],)
      print(energies, file=log)

      #print result.bond_deviations(sites_cart,
      #                             self.amber_structs.parm,
      #                             ignore_hd=False,
      #                             verbose=1,
      #)
      #print result.angle_deviations(sites_cart,
      #                             self.amber_structs.parm,
      #                             ignore_hd=False,
      #                             verbose=1,
      #)
      # print(len(result.gradients),list(result.gradients)[:10])

    return result

  def select(self, selection=None, iselection=None):
    if selection is not None: self.n_sites_cart = len(selection)
    if iselection is not None: assert self.n_sites_cart
    return_standard_grm = False
    if selection and False in selection: return_standard_grm = True
    if iselection and self.n_sites_cart!=len(iselection):
      return_standard_grm = True
    if return_standard_grm:
      result = self.standard_geometry_restraints_manager.select(selection=selection,
                                                                iselection=iselection)
      return result
    return self

  def cleanup(self):
    make_header('Cleaning up - Amber')
    if self.sander and self.amber_structs:
      if self.amber_structs.is_LES:
        import sanderles; sanderles.cleanup()
      else:
        import sander; sander.cleanup()

def digester(standard_geometry_restraints_manager,
             params,
             log=StringIO(),
             ):
  sgrm = standard_geometry_restraints_manager
  agrm = manager(params, log=log)
  for attr, value in list(vars(sgrm).items()):
    if attr.startswith('__'): continue
    setattr(agrm, attr, value)
  agrm.standard_geometry_restraints_manager = sgrm
  return agrm
