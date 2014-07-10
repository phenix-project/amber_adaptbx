import os, sys
from cctbx.array_family import flex

master_phil_str = """
  use_afitt = False
    .type = bool
  ligand_name = None
    .type = str  
  ff = 'mmff'
    .type = str  
"""
  
#~ class afitt_object:
    #~ def __init__(
            #~ self,
            #~ cif_object=None,
            #~ header_object=None):
      #~ self.cif_object = cif_object
      #~ self.header_object = header_object
  
#~ def read_cif_afitt(ligand_name):
  #~ from mmtbx import monomer_library
  #~ import mmtbx.monomer_library.server
  #~ mon_lib_srv = monomer_library.server.server()
  #~ get_func = getattr(mon_lib_srv, "get_comp_comp_id", None)
  #~ if (get_func is not None): 
    #~ ml=get_func(comp_id=ligand_name)
  #~ else:  
    #~ ml=mon_lib_srv.get_comp_comp_id_direct(comp_id=ligand_name)
  #~ bonds = []
  #~ for bond in ml.bond_list:
    #~ bonds.append([bond.atom_id_1, bond.atom_id_2])
  #~ for atom1, atom2 in bonds:
    #~ print atom1, atom2
  #~ return ml  

#~ [(2, 0, 'coval'), (1, 0, 'coval'), (3, 2, 'coval'), (4, 2, 'coval'), (5, 2, 'coval')]


class afitt_object:
  def __init__(self, ligand_path, pdb_hierarchy, ff='mmff'):
    self.n_atoms = 0
    self.resname = None
    self.chain = None
    self.number = 1 
    self.charge = 0
    self.atom_charges = None
    self.atom_elements = None
    self.bonds = None
    self.nbonds = 0
    self.sites_cart_ptrs = None
    self.total_model_atoms = 0
    self.ff = ff
    
    cif_object = self.read_cif_file(ligand_path)
    self.process_cif_object(cif_object, pdb_hierarchy)

  def read_cif_file(self, ligand_path):
    from iotbx import cif
    cif_object = cif.reader(file_path=ligand_path, strict=False).model()
    return cif_object
    
  def get_sites_cart_pointers(self, atom_ids, pdb_hierarchy, resname):
    sites_cart_ptrs=[0]*len(atom_ids)
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        for conformer in chain.conformers():
          for residue in conformer.residues():
            if residue.resname == resname:
              self.chain = chain.id
              for atom in residue.atoms():
                for atom_id in atom_ids:
                  if atom.name.strip() == atom_id.strip():
                    loc=atom_ids.index(atom_id)
                    sites_cart_ptrs[loc] = atom.i_seq
    return sites_cart_ptrs                
  
  def process_cif_object(self, cif_object, pdb_hierarchy):
    self.resname = cif_object['comp_list']['_chem_comp.id'][0]
    comp_rname='comp_%s' %self.resname
    assert comp_rname == cif_object.keys()[1]
    self.n_atoms = \
      int(cif_object['comp_list']['_chem_comp.number_atoms_all'][0])
    self.atom_charges =  \
      [float(i) for i in cif_object[comp_rname]['_chem_comp_atom.partial_charge']] 
    self.atom_elements = \
      [i for i in cif_object[comp_rname]['_chem_comp_atom.type_symbol']] 
    atom_ids = \
      [i for i in cif_object[comp_rname]['_chem_comp_atom.atom_id']] 
    bond_atom_1 = \
      [atom_ids.index(i) for i in cif_object[comp_rname]['_chem_comp_bond.atom_id_1']] 
    bond_atom_2 = \
      [atom_ids.index(i) for i in cif_object[comp_rname]['_chem_comp_bond.atom_id_2']]
    bond_dict={'single':1, 'double':2, 'triple':3, 'aromatic':4, 'coval':1}         
    bond_type = \
      [bond_dict[i] for i in cif_object[comp_rname]['_chem_comp_bond.type']] 
    self.bonds = zip(bond_atom_1, bond_atom_2, bond_type)
    self.charge = sum( self.atom_charges )
    self.nbonds = len(self.bonds)
    self.sites_cart_ptrs = self.get_sites_cart_pointers(
          atom_ids, 
          pdb_hierarchy,
          self.resname)
    self.total_model_atoms=pdb_hierarchy.atoms_size()
    #~ import code; code.interact(local=dict(globals(), **locals()))      

  def make_afitt_input(self, sites_cart, afitt_input):
    f=open(afitt_input,'w')
    f.write('%d\n' %self.n_atoms)
    f.write('residue_type %s chain %s number %d total_charge %d\n' 
            %(self.resname, self.chain, self. number, self.charge))
    assert len(self.atom_elements) == len(self.sites_cart_ptrs)
    for atom,ptr in zip(self.atom_elements, self.sites_cart_ptrs):
      f.write('%s   %20.16f   %20.16f   %20.16f\n' %(atom, 
            sites_cart[ptr][0], sites_cart[ptr][1], sites_cart[ptr][2]) )
    f.write('bond_table_nbonds %d\n' %self.nbonds)
    for bond in self.bonds:
      f.write('%d %d %d\n' %(bond[0], bond[1], bond[2]))       
    f.close()
  
def call_afitt(afitt_input, afitt_output, ff):
  from subprocess import Popen
  #~ os.system('buster_helper_%s <%s >afitt_out' %(ff, file_name))
  sts = Popen("buster_helper_%s" %ff + " <%s" %afitt_input + 
              " >%s" %afitt_output + " 2>trash", shell=True).wait()
  sts = Popen("rm trash", shell=True)

def process_afitt_output(afitt_output, geometry, afitt_object):
  ptrs = afitt_object.sites_cart_ptrs
  afitt_gradients = flex.vec3_double()
  with open(afitt_output, 'r') as afitt_o:
    for line in afitt_o:
      if line.startswith('ENERGYTAG'):
         afitt_energy=float(line.split()[1])
      elif line.startswith('GRADIENTTAG'):
         afitt_gradients.append (
            (float(line.split()[1]), 
             float(line.split()[2]), 
             float(line.split()[3]) ) ) 
  geometry.residual_sum += afitt_energy
  if (geometry.gradients is not None):  
    assert afitt_gradients.size() == len(ptrs)
    for afitt_gradient, ptr in zip(afitt_gradients, ptrs):
      geometry.gradients[ptr] = afitt_gradient
  return geometry


  

