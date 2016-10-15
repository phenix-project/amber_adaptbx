from __future__ import print_function
import sys
import numpy as np

"""
phenix.python print_atom_map_LES.py amber_phenix_atom_order_map.txt
"""

try:
    saved_map_file = sys.argv[1]
except IndexError:
    raise IndexError("must provide pickled file: amber_phenix_atom_order_map.txt")

arr = np.loadtxt(saved_map_file, dtype='i4').transpose()
atom_order_dict = dict(a2p=arr[0], p2a=arr[1])

print('phenix to amber order', atom_order_dict['p2a'].tolist())
print('amber to phenix order', atom_order_dict['a2p'].tolist())
