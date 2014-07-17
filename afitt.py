import os, sys
from cctbx.array_family import flex

master_phil_str = """
  use_afitt = False
    .type = bool
  ligand_file = None
    .type = str
  ligand_names = None
    .type = str
  ff = 'mmff'
    .type = str  
  scale = 'gnorm'
    .type = str  
"""
  

class afitt_object:
  def __init__(self, ligand_path,ligand_names, pdb_hierarchy, 
                ff='mmff', scale='gnorm'):
    self.n_atoms = []
    self.resname = ligand_names
    self.res_ids = [] #[chain, altloc,resseq]
    self.charge = []
    self.partial_charges = []
    self.atom_elements = []
    self.bonds = []
    self.nbonds = []
    self.sites_cart_ptrs = []
    self.formal_charges = []
    self.total_model_atoms = 0
    self.ff = ff
    self.scale = scale
    
    cif_object = self.read_cif_file(ligand_path)
    self.process_cif_object(cif_object, pdb_hierarchy)
    

  def read_cif_file(self, ligand_path):
    from iotbx import cif
    cif_object = cif.reader(file_path=ligand_path, strict=False).model()
    return cif_object
    
  def get_sites_cart_pointers(self, atom_ids, pdb_hierarchy, chain_id, altloc, resseq):
    sites_cart_ptrs=[0]*len(atom_ids)
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        if chain.id == chain_id:
          for conformer in chain.conformers():
            if conformer.altloc == altloc:
              for residue in conformer.residues():
                if residue.resseq == resseq:
                  for atom in residue.atoms():
                    for atom_id in atom_ids:
                      if atom.name.strip() == atom_id.strip():
                        loc=atom_ids.index(atom_id)
                        sites_cart_ptrs[loc] = atom.i_seq
    return sites_cart_ptrs

  def get_res_ids(self, pdb_hierarchy, resname):
    ids=[]
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        for conformer in chain.conformers():
          for residue in conformer.residues():
            if residue.resname == resname:
              ids.append([chain.id,conformer.altloc,residue.resseq])
    return ids
  
  def process_cif_object(self, cif_object, pdb_hierarchy):
    for res in self.resname:
      for i, id in enumerate(cif_object['comp_list']['_chem_comp.id']):
        if res == id:
          self.n_atoms.append(
            int(cif_object['comp_list']['_chem_comp.number_atoms_all'][i]) )
      comp_rname='comp_%s' %res
      assert cif_object.has_key(comp_rname), "Residue %s not in cif file!" %res
      self.partial_charges.append(
        [float(i) for i in cif_object[comp_rname]['_chem_comp_atom.partial_charge']]
        )
      self.atom_elements.append(
        [i for i in cif_object[comp_rname]['_chem_comp_atom.type_symbol']] 
        )
      atom_ids = \
        [i for i in cif_object[comp_rname]['_chem_comp_atom.atom_id']] 
      bond_atom_1 = \
        [atom_ids.index(i) for i in cif_object[comp_rname]['_chem_comp_bond.atom_id_1']] 
      bond_atom_2 = \
        [atom_ids.index(i) for i in cif_object[comp_rname]['_chem_comp_bond.atom_id_2']]
      bond_dict={'single':1, 'double':2, 'triple':3, 'aromatic':4, 'coval':1}         
      bond_type = \
        [bond_dict[i] for i in cif_object[comp_rname]['_chem_comp_bond.type']] 
      self.bonds.append( zip(bond_atom_1, bond_atom_2, bond_type) )
      self.charge.append( sum(self.partial_charges[-1]) )
      self.nbonds.append ( len(self.bonds[-1]) )
      res_ids = self.get_res_ids(pdb_hierarchy, res)
      self.res_ids.append(res_ids)
      this_res_sites_cart_ptrs=[]
      for residue_instance in self.res_ids[-1]:
        this_res_sites_cart_ptrs.append( self.get_sites_cart_pointers(
                                          atom_ids, 
                                          pdb_hierarchy,
                                          chain_id=residue_instance[0],
                                          altloc=residue_instance[1],
                                          resseq=residue_instance[2]) 
                                        )
      self.sites_cart_ptrs.append( this_res_sites_cart_ptrs )                           
      if cif_object[comp_rname].has_key('_chem_comp_atom.formal_charge'):
        self.formal_charges.append(
          [float(i) for i in cif_object[comp_rname]['_chem_comp_atom.formal_charge']] 
          )
      else:
        self.formal_charges.append([])
    self.total_model_atoms=pdb_hierarchy.atoms_size()
    #~ import code; code.interact(local=dict(globals(), **locals()))      
    
  def make_afitt_input(self, sites_cart, afitt_input, resname_i, instance_i):
    r_i=resname_i
    i_i=instance_i
    sites_cart_ptrs=self.sites_cart_ptrs[r_i][i_i]
    elements=self.atom_elements[r_i]
    f=open(afitt_input,'w')
    f.write('%d\n' %self.n_atoms[r_i])
    f.write('residue_type %s chain %s number %d total_charge %d\n' 
            %(self.resname[r_i], self.res_ids[r_i][i_i][0],1,self.charge[r_i] ))
    assert len(elements) ==  len(sites_cart_ptrs), \
            "No. of atoms in residue %s, instance %d does not equal to \
            number of atom seq pointers." %(self.resname[resname_i], instance_i)
    for atom,ptr in zip(elements, sites_cart_ptrs):
      f.write('%s   %20.16f   %20.16f   %20.16f\n' %(atom, 
            sites_cart[ptr][0], sites_cart[ptr][1], sites_cart[ptr][2]) )
    f.write('bond_table_nbonds %d\n' %self.nbonds[r_i])
    for bond in self.bonds[r_i]:
      f.write('%d %d %d\n' %(bond[0], bond[1], bond[2]))
    if self.formal_charges[r_i]:
      f.write("formal charges\n")
      for fcharge in self.formal_charges[r_i]:
        f.write ('%d\n' %fcharge)
    f.close()
  
def call_afitt(afitt_input, afitt_output, ff):
  from subprocess import Popen
  #~ os.system('buster_helper_%s <%s >afitt_out' %(ff, file_name))
  sts = Popen("buster_helper_%s" %ff + " <%s" %afitt_input + 
              " >%s" %afitt_output + " 2>trash", shell=True).wait()
  sts = Popen("rm trash", shell=True)

def process_afitt_output(afitt_output, geometry, afitt_object, 
                          resname_i, instance_i):
  r_i=resname_i
  i_i=instance_i                            
  ptrs = afitt_object.sites_cart_ptrs[r_i][i_i]
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
  ### debug_stuff
  print ("AFITT_ENERGY %s_%d: %10.4f\n" 
                  %(afitt_object.resname[r_i], 
                    int(afitt_object.res_ids[r_i][i_i][2]), 
                    afitt_energy ))
  ### end_debug
  geometry.residual_sum += afitt_energy
  
  #~ import inspect
  #~ for i in inspect.stack():
    #~ print i[1], i[2], i[4]
  #~ print "\n\n\n\n"  
  
  if (geometry.gradients is not None):
    assert afitt_gradients.size() == len(ptrs)
    if afitt_object.scale == 'gnorm':
      from math import sqrt
      phenix_norm=0
      afitt_norm=0
      for afitt_gradient, ptr in zip(afitt_gradients, ptrs):
        phenix_norm += geometry.gradients[ptr][0]**2+geometry.gradients[ptr][1]**2+geometry.gradients[ptr][2]**2
        afitt_norm += afitt_gradient[0]**2+afitt_gradient[1]**2+afitt_gradient[2]**2
      phenix_norm = sqrt(phenix_norm)
      afitt_norm = sqrt(afitt_norm)
      gr_scale = phenix_norm/afitt_norm
      ### debug_stuff
      print ("GRNORM_RATIO %s_%d: %10.4f\n" 
                    %(afitt_object.resname[r_i], 
                      int(afitt_object.res_ids[r_i][i_i][2]),
                      gr_scale ))
  
      ### end_debug
    elif afitt_object.scale == 'noafitt':
      gr_scale = None
    else:
      gr_scale = float(afitt_object.scale)    
    
    ### debug_stuff
    print_gradients = False
    if print_gradients:
      print("\n\nGRADIENTS BEFORE AFTER AFITT\n")
      print "NORMS: %10.4f         %10.4f\n" %(phenix_norm, afitt_norm)
      for afitt_gradient, ptr in zip(afitt_gradients, ptrs):
        print "(%10.4f %10.4f %10.4f) (%4.4f %4.4f %4.4f)" \
            %(geometry.gradients[ptr][0], geometry.gradients[ptr][1], geometry.gradients[ptr][2],
            afitt_gradient[0], afitt_gradient[1], afitt_gradient[2])
    ### end_debug
    if gr_scale:
      for afitt_gradient, ptr in zip(afitt_gradients, ptrs):
        scaled_gradient = (afitt_gradient[0]*gr_scale, 
                         afitt_gradient[1]*gr_scale,
                         afitt_gradient[2]*gr_scale)
        geometry.gradients[ptr] = scaled_gradient
  return geometry

def get_afitt_energy(cif_file, ligand_names, pdb_hierarchy, ff, sites_cart):
  afitt_o = afitt_object(
                cif_file,
                ligand_names,
                pdb_hierarchy,
                ff)
  afitt_input='afitt_in'
  afitt_output='afitt_out'
  energies=[]
  for resname_i,resname in enumerate(afitt_o.resname):
    for instance_i, instance in enumerate(afitt_o.res_ids[resname_i]):
      #~ import code; code.interact(local=dict(globals(), **locals())) 
      afitt_o.make_afitt_input(sites_cart, afitt_input, resname_i, instance_i)
      call_afitt(afitt_input, afitt_output, ff)  
      with open(afitt_output, 'r') as afitt_out:
        for line in afitt_out:
          if line.startswith('ENERGYTAG'):
            energy=float(line.split()[1])
      energies.append([resname, int(instance[2]), energy] )
  return energies

def run(pdb_file, cif_file, ligand_names, ff='mmff'):
  import iotbx.pdb
  pdb_inp = iotbx.pdb.input(file_name=pdb_file)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  pdb_hierarchy.atoms().reset_i_seq()
  xrs = pdb_hierarchy.extract_xray_structure()
  sites_cart=xrs.sites_cart()
  energies = get_afitt_energy(cif_file, ligand_names, pdb_hierarchy, ff, sites_cart)
  for energy in energies:
    print "%s_%d AFITT_ENERGY: %10.4f" %(energy[0], energy[1], energy[2])
 
if (__name__ == "__main__"):
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument("pdb_file", help="pdb file")
  parser.add_argument("cif_file", help="cif file", default=0)
  parser.add_argument("ligand_names", help="3-letter ligand names separated by commas")
  parser.add_argument("-ff", help="afitt theory: mmff, pm3 or am1", default='mmff')
  args = parser.parse_args()
  ligand_names=args.ligand_names.split(',')
  run(args.pdb_file, args.cif_file, ligand_names, args.ff)

