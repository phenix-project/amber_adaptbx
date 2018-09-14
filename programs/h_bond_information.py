import os, sys

from libtbx import easy_run
import libtbx.load_env
from libtbx.utils import Sorry

from amber_adaptbx.program_template import ProgramTemplate

amber_dir = libtbx.env.dist_path("amber")

input_str = '''parm %(preamble)s.pdb
trajin %(preamble)s.pdb
hbond avgout %(preamble)s.hbond.dat
go
'''

def get_number_of_residues(hierarchy):
  atoms = hierarchy.select(hierarchy.get_peptide_c_alpha_selection()).atoms()
  return len(atoms)

class Program(ProgramTemplate):

  description = '''
  Enumerate the H-bonds in a model
    - requires hydrogen atoms be present
    - requires AmberTools to be installed
'''
  datatypes = ['phil', 'model']

  master_phil_str = '''

  output
  {
    show_number_of_h_bonds_per_1000 = True
      .type = bool
    list_h_bonds = False
      .type = bool
  }
  '''

  def validate(self):
    pass

  def run(self):
    model = self.data_manager.get_model()
    if not model.has_hd():
      raise Sorry('Model has no hydrogen atoms! Please add H/D and try again.')

    cpptraj_exe = os.path.join(amber_dir,
                               'bin',
                               'cpptraj'
                               )
    preamble = self.data_manager.get_model_names()[0].split('.')[0]
    print '''
    Running cpptraj executable
      %s
    on input
"""
%s
"""
    ''' % ( cpptraj_exe,  input_str % locals())
    rc = easy_run.go(command=cpptraj_exe,
                     stdin_lines=input_str % locals(),
                     )
    rc.show_stdout()

    print '  Writing H-bonds to %(preamble)s.hbond.dat' % locals()

    f=file('%(preamble)s.hbond.dat' % locals(), 'rb')
    lines = f.read()
    f.close()

    if self.params.output.list_h_bonds:
      print
      print lines

    nr = get_number_of_residues(model.get_hierarchy())

    if self.params.output.show_number_of_h_bonds_per_1000:
      h_bonds = len(lines.splitlines())-1
      atoms = len(model.get_hierarchy().atoms())
      print '_'*80
      print '''
      Number of H-bonds found : %7d
      Number of atoms         : %7d
      H-bonds per 1000 atoms  : %10.2f
      H-bonds per residues    : %10.2f
      ''' % (h_bonds,
             atoms,
             h_bonds/atoms*1000,
             h_bonds/nr,
             )


