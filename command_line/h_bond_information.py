from __future__ import division, print_function

from iotbx.cli_parser import run_program
from amber_adaptbx.programs import h_bond_information

if __name__ == '__main__':
  run_program(program_class=h_bond_information.Program)
