from __future__ import division
from __future__ import print_function
from cctbx import geometry_restraints
import sys
import amber_adaptbx.lbfgs

try:
    xrange
except NameError:
    xrange = range

class lbfgs(amber_adaptbx.lbfgs.lbfgs):
  def __init__(self,
        sites_cart,
        geometry_restraints_manager,
        geometry_restraints_flags,
        lbfgs_termination_params,
        amber_structs,
        sites_cart_selection=None,
        lbfgs_exception_handling_params=None,
        grms_termination_cutoff=0,
        site_labels=None):
    self.grms_termination_cutoff = grms_termination_cutoff
    amber_adaptbx.lbfgs.lbfgs.__init__(self,
      sites_cart=sites_cart,
      geometry_restraints_manager=geometry_restraints_manager,
      geometry_restraints_flags=geometry_restraints_flags,
      lbfgs_termination_params=lbfgs_termination_params,
      sites_cart_selection=sites_cart_selection,
      lbfgs_exception_handling_params=lbfgs_exception_handling_params,
      site_labels=site_labels,
      amber_structs=amber_structs)

  def callback_after_step(self, minimizer):
    self.apply_shifts()
    if self.grms <self.grms_termination_cutoff:
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
               parallelity                    = False,
               grms_termination_cutoff        = 0,
               alternate_nonbonded_off_on     = False,
               log                            = None,
               prmtop                         = None,
               ambcrd                           = None):
    # parallelity needs to be propagated - NWM
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
      bond_similarity    = True)

    import sander, sanderles
    amber_structs = amber_adaptbx.SanderStruct(
      parm_file_name=prmtop,
      rst_file_name=ambcrd)
    amber_structs.sander_engine.setup(amber_structs.parm,
             amber_structs.parm.coordinates,
             amber_structs.parm.box,
             amber_structs.inp)
    self.show(amber_structs, log=log)

    for i_macro_cycle in xrange(number_of_macro_cycles):
      print("  macro-cycle:", i_macro_cycle, file=self.log)
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
        grms_termination_cutoff         = grms_termination_cutoff,
        site_labels                     = None,
        amber_structs                   = amber_structs)
      self.pdb_hierarchy.atoms().set_xyz(sites_cart)
      self.show(amber_structs)
      geometry_restraints_flags.nonbonded = nonbonded
      lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
          max_iterations = max_number_of_iterations)

    amber_structs.sander_engine.cleanup()

  def show(self, amber_structs, log=None):
    import amber_adaptbx as amber
    amber_geometry_manager=amber.geometry_manager(
       sites_cart=self.pdb_hierarchy.atoms().extract_xyz(),
       amber_structs=amber_structs)
    amber_geometry=amber_geometry_manager.energies_sites(
       self.restraints_manager.geometry.crystal_symmetry,
       log=log,
       )
    amber_geometry.show()
