
import iotbx.pdb
import sys

#======================================================================#
# When tleap creates an Amber topology and ambpdb recreates a pdb with #
# atom order that matches the topology, information about chain name,  #
# b-factors and occupancy is lost. This is a quick jiffy script (will  #
# most likely fail on more complex cases like multiple residues with   #
# same seq id) to replace the missing information using the original   #
# file.                                                                #
########################################################################



def run (pre_amber_pdb_file, post_amber_pdb_file, output_pdb_file="out.pdb"):
  pdb_pre=iotbx.pdb.input(file_name=pre_amber_pdb_file)
  pdb_h_pre = pdb_pre.construct_hierarchy()
  pdb_post=iotbx.pdb.input(file_name=post_amber_pdb_file)
  pdb_h_post = pdb_post.construct_hierarchy()

  for atom_group in pdb_h_post.atom_groups():
    if atom_group.resname in ['HID','HIP','HIE']:
      atom_group.resname = "HIS"
    elif atom_group.resname == "CYX":
      atom_group.resname = "CYS"
    elif atom_group.resname in ["WAT"]:
      atom_group.resname = "HOH"
      
  # even though we're not modifying pdb_pre, we need to change residue 
  # names to ensure matching below
  for atom_group in pdb_h_pre.atom_groups():
    if atom_group.resname in ['HID','HIP','HIE']:
      atom_group.resname = "HIS"
    elif atom_group.resname == "CYX":
      atom_group.resname = "CYS"
    elif atom_group.resname in ["WAT"]:
      atom_group.resname = "HOH"      

  #match residues based on resseq and resname
  #match atoms based on name an i_seq
  for chain_post in pdb_h_post.chains():
    for resi_post in chain_post.conformers()[0].residues():
      for atom_post in resi_post.atoms():
        for chain_pre in pdb_h_pre.chains():
          for resi_pre in chain_pre.conformers()[0].residues():
              if resi_pre.resseq==resi_post.resseq and resi_pre.resname.strip()==resi_post.resname.strip():
                for atom_pre in resi_pre.atoms():
                  if atom_pre.name == atom_post.name and atom_pre.i_seq==atom_post.i_seq:
                    atom_post.b=atom_pre.b
                    atom_post.occ=atom_pre.occ
                    chain_post.id=chain_pre.id

  pdb_h_post.write_pdb_file(file_name=output_pdb_file, append_end=True,crystal_symmetry=pdb_pre.crystal_symmetry())
    
if __name__=="__main__":
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument("pre_amber_pdb_file", help="name of pre tleap file")
  parser.add_argument("post_amber_pdb_file", help="name of post tleap file")
  parser.add_argument("output_filename", help="output filename", default="out.pdb")
  args = parser.parse_args()
  run(args.pre_amber_pdb_file, args.post_amber_pdb_file, args.output_filename)
