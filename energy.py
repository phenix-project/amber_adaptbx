from __future__ import division
from __future__ import print_function
from math import sqrt
import scitbx.restraints
from libtbx.utils import Sorry
try:
  import sander
  import sanderles
  import parmed # require version in AmberTools >= 16
except ImportError as e:
  sander = None
  sanderles = None
  parmed = None
  
from amber_adaptbx.utils import (
    bond_rmsd,
    angle_rmsd,
    is_prmtop_LES,
    check_file,
)

__all__ = ['energies', 'SanderStruct']


class energies(scitbx.restraints.energies):

  def __init__(self, *args, **kwds):
    scitbx.restraints.energies.__init__(self, *args, **kwds)
    self.energy_components = None
    self.amber = True
    # import code; code.interact(local=dict(globals(), **locals()))
    # sys.exit()

  def show(self):
    # assert 0 # not writing to log...
    print("    Amber total energy: %0.2f" % (self.residual_sum))
    print("      bonds (n=%d): %0.2f" % (self.energy_components[6],
                                         self.energy_components[1]))
    print("      angles (n=%d): %0.2f" % (self.energy_components[7],
                                          self.energy_components[2]))
    print("      dihedrals (n=%d): %0.2f" % (self.energy_components[8],
                                             self.energy_components[3]))
    print("      electrostatics: %0.2f" % (self.energy_components[4]))
    print("      van der Waals: %0.2f" % (self.energy_components[5]))
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
  def angle_deviations(self,
                       sites_cart,
                       parm,
                       ignore_hd=False,
                       get_deltas=False,
                       get_extremes=False,
                       verbose=False,
                       ):
    return angle_rmsd(parm, sites_cart, ignore_hd, get_deltas, get_extremes, verbose)

  def bond_deviations(self,
                      sites_cart,
                      parm,
                      ignore_hd=False,
                      get_deltas=False,
                      get_extremes=False,
                      verbose=False
                      ):
    return bond_rmsd(parm, sites_cart, ignore_hd, get_deltas, get_extremes, verbose)

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


class SanderStruct(object):

  def __init__(self, parm_file_name, rst_file_name, ridingH=True,
        igb=0, cut=8.0 ):
    check_file("amber.topology_file_name", parm_file_name)
    check_file("amber.coordinate_file_name", rst_file_name)
    self.md_engine = 'sander'
    self.parm = parmed.load_file(parm_file_name, xyz=rst_file_name)
    # where do we need this self.rst?
    self.ridingH = ridingH
    self.is_LES = is_prmtop_LES(parm_file_name)
    self.sander_engine = sanderles if self.is_LES else sander

    self.order_converter = None
    self.order_map_file_name = None
    if igb == 0:
       self.inp = self.sander_engine.pme_input()
    else:
       self.inp = self.sander_engine.gas_input(igb)
       self.inp.ntb = 1
       self.inp.cut = cut
       self.inp.rgbmax = 15.0
       print( "Setting Amber igb, cut, rgbmax: %d, %6.2f  15.00\n" % 
            ( igb, cut ) )
    self.writer = None
