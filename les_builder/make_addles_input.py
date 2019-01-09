import os
import sys
import math
import parmed
import iotbx.pdb

"""
Example: phenix.python make_addles_input.py 2igd.pdb

2igd.pdb is original ASU pdb file
"""

__all__ = ['main']

def get_LES_residue_dict(parm):
  """
  
  Parameters
  ----------
  parm : parmed.Structure

  Returns
  -------
  out : {resid : number_of_conformers}
      resid is 1-based index
  """
  holder = {}
  holderbb = {}  # flag for N,H,C,O
  holderca = {}  # flag for CA alternates
  holdercb = {}  # flag for side-chain alternates that include CB
  holdercg = {}  # flag for side-chain alternates that include CG or OG
  holder_atoms = {}  # place to store atom names in this residue that
                     # have alternate locations
  for atom in parm.atoms:
    if atom.other_locations:
      holder[atom.residue.idx+1] = len(atom.other_locations)
      if atom.residue.idx+1 in holder_atoms:
         holder_atoms[atom.residue.idx+1] += " " + atom.name
      else:
         holder_atoms[atom.residue.idx+1] = atom.residue.name + ": " + atom.name
      if atom.name in ["N", "C", "O", "H" ]:
         holderbb[atom.residue.idx+1] = 1
      elif atom.name in ["CA" ]:
        # print "found CA altlocs in residue %d" % int(atom.residue.idx+1)
        for key in atom.other_locations:
           cadist = math.sqrt( (atom.xx - atom.other_locations[key].xx)**2 +
                               (atom.xy - atom.other_locations[key].xy)**2 +
                               (atom.xz - atom.other_locations[key].xz)**2 )
           # print "parent: %8.3f %8.3f %8.3f" % ( atom.xx, atom.xy, atom.xz )
           # print "child : %8.3f %8.3f %8.3f" % ( 
           #         atom.other_locations[key].xx,
           #         atom.other_locations[key].xy,
           #         atom.other_locations[key].xz )
           # print "distance: %8.3f" % cadist
           break
        holder_atoms[atom.residue.idx+1] += "%8.3f" % cadist
        holderca[atom.residue.idx+1] = cadist

      elif atom.name in ["CB" ]:
        # print "found CB altlocs in residue %d" % int(atom.residue.idx+1)
        for key in atom.other_locations:
           cbdist = math.sqrt( (atom.xx - atom.other_locations[key].xx)**2 +
                               (atom.xy - atom.other_locations[key].xy)**2 +
                               (atom.xz - atom.other_locations[key].xz)**2 )
           # print "parent: %8.3f %8.3f %8.3f" % ( atom.xx, atom.xy, atom.xz )
           # print "child : %8.3f %8.3f %8.3f" % ( 
           #         atom.other_locations[key].xx,
           #         atom.other_locations[key].xy,
           #         atom.other_locations[key].xz )
           # print "distance: %8.3f" % cbdist
           break
        holder_atoms[atom.residue.idx+1] += "%8.3f" % cbdist
        holdercb[atom.residue.idx+1] = cbdist

      elif atom.name in ["CG", "OG", "CG1", "CG2" ]:
        # print "found CG altlocs in residue %d" % int(atom.residue.idx+1)
        for key in atom.other_locations:
           cgdist = math.sqrt( (atom.xx - atom.other_locations[key].xx)**2 +
                               (atom.xy - atom.other_locations[key].xy)**2 +
                               (atom.xz - atom.other_locations[key].xz)**2 )
           # print "parent: %8.3f %8.3f %8.3f" % ( atom.xx, atom.xy, atom.xz )
           # print "child : %8.3f %8.3f %8.3f" % ( 
           #         atom.other_locations[key].xx,
           #         atom.other_locations[key].xy,
           #         atom.other_locations[key].xz )
           # print "distance: %8.3f" % cgdist
           break
        holder_atoms[atom.residue.idx+1] += "%8.3f" % cgdist
        holdercg[atom.residue.idx+1] = cgdist

  return holder, holderbb, holderca, holdercb, holdercg, holder_atoms
  
def addles_input(pdb_fn='2igd.pdb', prmtop=None, rst7_file=None):
  root_name = os.path.basename(pdb_fn).split('.')[0]
  uc_parm = root_name + 'a.parm7' if prmtop is None else prmtop
  uc_rst7 = root_name + 'a.rst7' if rst7_file is None else rst7_file
  
  uc_les_parm = '4amber_' + root_name + '.LES.prmtop'
  uc_les_rst7 = '4amber_' + root_name + '.LES.rst7'
  
  parm = parmed.load_file(pdb_fn)
  
  crystal_symmetry = iotbx.pdb.input(pdb_fn).crystal_symmetry()
  n_asu = len(crystal_symmetry.space_group().all_ops()) 
  
  n_asu_residues = len(parm.residues)
  # n_uc_residues = n_asu * n_asu_residues
  # print('n_asu', n_asu, 'n_asu_residues', n_asu_residues)
  # print('n_uc_residues', n_uc_residues)
  
  commands = []
  
  header = """
  file rprm name=({uc_parm}) read
  file rcbd name=({uc_rst7}) read
  file wprm name=({uc_les_parm}) wovr
  file wcrd name=({uc_les_rst7}) wovr
  action
  omas
  """.format(uc_parm=uc_parm,
             uc_rst7=uc_rst7,
             uc_les_parm=uc_les_parm,
             uc_les_rst7=uc_les_rst7)

  header = '\n'.join((line.strip() for line in header.split('\n')))
  commands.append(header.strip())
  
  (holder, holderbb, holderca, holdercb, holdercg, holder_atoms) = get_LES_residue_dict(parm)
  
  for chain_id in range(n_asu):
    commands.append('~ protein chain: {}'.format(chain_id+1))
  
    for resid in sorted(holder.keys()):
      n_conformers = holder[resid] + 1
      line_template = ''
      if resid in holderca.keys():
        # choose #sia if all CA atoms are within 0.15 Ang, and there are
        # no backbone alternates:
        if holderca[resid] > 0.15 or resid in holderbb:
           line_template = 'spac numc={numc} pick #cca {resid} {resid} done'
        else:
           line_template = 'spac numc={numc} pick #sia {resid} {resid} done'
      elif resid in holdercb.keys():
        # choose #sid if all CB atoms are within 0.15 Ang
        if holdercb[resid] > 0.15:
           line_template = 'spac numc={numc} pick #sia {resid} {resid} done'
        else:
           line_template = 'spac numc={numc} pick #sid {resid} {resid} done'
      elif resid in holdercg.keys():
        # choose #sig if all CG atoms are within 0.15 Ang
        if holdercg[resid] > 0.15:
           line_template = 'spac numc={numc} pick #sid {resid} {resid} done'
        else:
           line_template = 'spac numc={numc} pick #sig {resid} {resid} done'
      elif not resid in holderbb:
        line_template = 'spac numc={numc} pick #sig {resid} {resid} done'

      if line_template:
         commands.append(line_template.format(numc=n_conformers,
                                           resid=resid+chain_id*n_asu_residues))
         commands.append("~ " + holder_atoms[resid])
  
  commands.append("*EOD\n")
  return commands

def main(pdb_fn='2igd.pdb'):
  for cm in addles_input(pdb_fn=pdb_fn):
    print(cm)

if __name__ == '__main__':
  pdb_fn = sys.argv[1]

  main(pdb_fn=pdb_fn)
