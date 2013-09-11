#! /usr/bin/env phenix.python
from libtbx import group_args
import sys, os
import iotbx.pdb
import argparse
from scitbx.array_family import flex
#~ sys.path.append('/home/pjanowsk/amberSD/AmberTools/src/xtalutil/Phenix')
#~ import amber_adaptbx_ext as PAI
import boost.python
ext = boost.python.import_ext("amber_adaptbx_ext")


class geometry_manager(object):

  def __init__(self,
        prmtop=None,
        ambcrd=None,
        sites_cart=None,
        energy_components=flex.double([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]),
        gradients=None,
        number_of_restraints=100):
    self.prmtop = prmtop
    self.ambcrd = ambcrd
    self.sites_cart = sites_cart
    self.energy_components = energy_components
    self.gradients=flex.double(len(sites_cart)*3)
    self.number_of_restraints=number_of_restraints


  def energies_sites(self):
    #Convert flex arrays to C arrays
    sites_cart_c=ext.ExtractVec(self.sites_cart.as_double())
    gradients_c=ext.ExtractVec(self.gradients)
    energy_components_c=ext.ExtractVec(self.energy_components)

    # Call c++ interface to call mdgx to calculate new gradients and target
    ext.callMdgx(sites_cart_c, gradients_c, energy_components_c, self.prmtop, self.ambcrd)
    # Convert back into python types (eg. into flex arrays for phenix to use)
    self.gradients=flex.vec3_double(gradients_c)
    #~ print "\nGRADIENTS"  
    #~ print list(self.gradients[0:3])
    self.energy_components=flex.double(energy_components_c)
    self.residual_sum= float(energy_components_c[0])
    return self

  def show(self):
    print "\n\n"
    print "Amber_total_energy: %7.6f" 		%(self.residual_sum)
    print "  bonds (n=%d): %7.6f" 			%(self.energy_components[6], self.energy_components[1])
    print "  angles (n=%d): %7.6f" 			%(self.energy_components[7], self.energy_components[2])
    print "  dihedrals (n=%d): %7.6f" 		%(self.energy_components[8], self.energy_components[3])
    print "  electrostatics: %7.6f" 		%(self.energy_components[4])
    print "  vanderWaals: %7.6f" 			%(self.energy_components[5])
    print "\n\n"
    return 0

def print_sites_cart(sites_cart):
	for atom in sites_cart:
		print("%8.3f%8.3f%8.3f"%(atom[0], atom[1], atom[2]))




