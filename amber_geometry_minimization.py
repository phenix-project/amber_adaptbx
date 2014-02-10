from __future__ import division
import iotbx.phil
from cctbx import geometry_restraints
import cctbx.geometry_restraints.lbfgs
import scitbx.lbfgs
import sys
import amber_adaptbx.lbfgs
from amber_adaptbx import get_amber_structs

class lbfgs(amber_adaptbx.lbfgs.lbfgs):
  def __init__(self,
        sites_cart,
        geometry_restraints_manager,
        geometry_restraints_flags,
        lbfgs_termination_params,
        mdgx_structs,
        sites_cart_selection=None,
        lbfgs_exception_handling_params=None,
        grmsd_termination_cutoff=0,
        site_labels=None):
    self.grmsd_termination_cutoff = grmsd_termination_cutoff
    amber_adaptbx.lbfgs.lbfgs.__init__(self,
      sites_cart=sites_cart,
      geometry_restraints_manager=geometry_restraints_manager,
      geometry_restraints_flags=geometry_restraints_flags,
      lbfgs_termination_params=lbfgs_termination_params,
      sites_cart_selection=sites_cart_selection,
      lbfgs_exception_handling_params=lbfgs_exception_handling_params,
      site_labels=site_labels,
      mdgx_structs=mdgx_structs)

  def callback_after_step(self, minimizer):
    self.apply_shifts()
    if self.rmsd_gradient <self.grmsd_termination_cutoff:
      return True

class run(object):
  def __init__(self,
               restraints_manager,
               pdb_hierarchy,
               max_number_of_iterations       = 500,
               number_of_macro_cycles         = 5,
               selection                      = None,
               bond                           = False,
               nonbonded                      = False,
               angle                          = False,
               dihedral                       = False,
               chirality                      = False,
               planarity                      = False,
               generic_restraints             = False,
               grmsd_termination_cutoff       = 0,
               alternate_nonbonded_off_on     = False,
               log                            = None,
               prmtop                         = None,
               ambcrd                           = None):
    self.log = log
    if self.log is None:
      self.log = sys.stdout
    self.pdb_hierarchy = pdb_hierarchy
    self.minimized = None
    self.restraints_manager = restraints_manager
    assert max_number_of_iterations+number_of_macro_cycles > 0
    assert [bond,nonbonded,angle,dihedral,chirality,planarity].count(False) < 6
    if(alternate_nonbonded_off_on and number_of_macro_cycles % 2 != 0):
      number_of_macro_cycles += 1
    import scitbx.lbfgs
    lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
      max_iterations = max_number_of_iterations)
    exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
      ignore_line_search_failed_step_at_lower_bound = True)
    geometry_restraints_flags = geometry_restraints.flags.flags(
      bond               = bond,
      nonbonded          = nonbonded,
      angle              = angle,
      dihedral           = dihedral,
      chirality          = chirality,
      planarity          = planarity,
      reference_dihedral = True,
      bond_similarity    = True,
      generic_restraints = True)

    mdgx_structs=get_amber_structs (prmtop, ambcrd)

    self.show(mdgx_structs)
    for i_macro_cycle in xrange(number_of_macro_cycles):
      print >> self.log, "  macro-cycle:", i_macro_cycle
      if(alternate_nonbonded_off_on and i_macro_cycle<=number_of_macro_cycles/2):
        geometry_restraints_flags.nonbonded = bool(i_macro_cycle % 2)
      sites_cart = self.pdb_hierarchy.atoms().extract_xyz()
      self.minimized = lbfgs(
        sites_cart                      = sites_cart,
        geometry_restraints_manager     = restraints_manager.geometry,
        geometry_restraints_flags       = geometry_restraints_flags,
        lbfgs_termination_params        = lbfgs_termination_params,
        lbfgs_exception_handling_params = exception_handling_params,
        sites_cart_selection            = selection,
        grmsd_termination_cutoff        = grmsd_termination_cutoff,
        site_labels                     = None,
        mdgx_structs                    = mdgx_structs)
      self.pdb_hierarchy.atoms().set_xyz(sites_cart)
      self.show(mdgx_structs)
      geometry_restraints_flags.nonbonded = nonbonded
      lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
          max_iterations = max_number_of_iterations)

  def show(self, mdgx_structs):
    import amber_adaptbx as amber
    amber_geometry_manager=amber.geometry_manager(
       sites_cart=self.pdb_hierarchy.atoms().extract_xyz(),
       mdgx_structs=mdgx_structs)
    amber_geometry=amber_geometry_manager.energies_sites()
    amber_geometry.show()
