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
from amber_adaptbx_ext import *

class geometry_manager(object):

  def __init__(self,
        prmtop=None,
        ambcrd=None,
        sites_cart=None,
        residual_sum=flex.double([0.0]),
        gradients=None,
        number_of_restraints=100):
    self.prmtop = prmtop
    self.ambcrd = ambcrd
    self.sites_cart = sites_cart
    self.residual_sum = residual_sum
    self.gradients=flex.double(len(sites_cart)*3)
    self.number_of_restraints=number_of_restraints


  def energies_sites(self):
    #Convert flex arrays to C arrays
    sites_cart_c=ExtractVec(self.sites_cart.as_double())
    gradients_c=ExtractVec(self.gradients)
    target_c=ExtractVec(self.residual_sum)
    
    # Call c++ interface to call mdgx to calculate new gradients and target
    callMdgx(sites_cart_c, gradients_c, target_c, self.prmtop, self.ambcrd)
    # Convert back into python types (eg. into flex arrays for phenix to use)
    self.gradients=flex.vec3_double(gradients_c)
    self.residual_sum= float(target_c[0])
    return self

  def show(self):
    print "\n\namber show\n\n"
    return 0

def print_sites_cart(sites_cart):
	for atom in sites_cart:
		print("%8.3f%8.3f%8.3f"%(atom[0], atom[1], atom[2]))




