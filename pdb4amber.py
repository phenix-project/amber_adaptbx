#!/usr/bin/env python
# Romain M. Wolf, NIBR Basel, December 2013

#  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#


# PDB analyzer to prepare protein(-ligand) PDB files for Amber simulations.

import os, sys
from optparse import OptionParser
from math import sqrt

# Global constants
RESPROT = ('ALA', 'ARG', 'ASN', 'ASP',
           'CYS', 'GLN', 'GLU', 'GLY',
           'HIS', 'ILE', 'LEU', 'LYS',
           'MET', 'PHE', 'PRO', 'SER',
           'THR', 'TRP', 'TYR', 'VAL',
           'HID', 'HIE', 'HIN', 'HIP',
           'CYX', 'ASH', 'GLH', 'LYH',
           'ACE', 'NME')

#=============================================
def pdb_read(pdbin):
#=============================================
# only records starting with the following strings are kept...
  ACCEPTED =  ('ATOM  ','HETATM', 'TER   ', 'END   ', 'MODEL ', 'ENDMDL')
# records starting with the following strings are considered 'dividers'
# and they are cleaned for the rest of the line...
  DIVIDERS =  ('TER   ', 'END   ', 'MODEL ', 'ENDMDL')
  records = open(pdbin, 'r')
  print "\n=================================================="
  print "Summary of pdb4amber for file %s"%pdbin
  print "=================================================="
  '''
  This PDB reader first splits PDB lines (records) into individual (named) fields,
  then later re-assembles the fields in records. This might be done more elegantly,
  but we keep this scheme for now for clarity (given the messy PDB file format)!

    record_type     0    atom_number     1    blank1           2
    atom_name       3    alt_loc_ind     4    residue_name     5
    blank2          6    chain_id        7    residue_number   8
    insertion_code  9    blank3         10    x               11
    y               12   z              13    occupancy       14
    bfactor         15   blank4         16    element         17
    charge          18

  '''
  record_type = []; atom_number = []; blank1 = []
  atom_name = []; alt_loc_ind = []; residue_name = []; blank2 = []
  chain_id = []; residue_number = []; insertion_code = []; blank3 = []
  x = []; y = []; z = []; occupancy = []
  bfactor = []; blank4 = []; element = []; charge = []
# keep everything in 'ACCEPTED'
  lines = records.readlines()
  for line in lines:
    # make all lines 80 characters long (fill up with blanks if necessary)
    # so that the line parser will not fail on shorter lines...
    line = (line.rstrip() + (80-len(line)) * ' ')
    if '%-6s' % line[0:6] not in ACCEPTED:
      continue
# make clean divider lines without additional stuff that might hurt later
    elif line[0:6] in DIVIDERS:
      line = line[0:6]
      line = (line.rstrip() + (80-len(line)) * ' ')
      print line
    else:
      pass
# split the line into records
    record_type.append('%-6s' % line[0:6])
    atom_number.append(line[6:11])
    blank1.append(line[11:12])
    atom_name.append(line[12:16])
    alt_loc_ind.append(line[16:17]); residue_name.append(line[17:20])
    blank2.append(line[20:21]);
    chain_id.append(line[21:22])
    residue_number.append(line[22:26])
    insertion_code.append(line[26:27])
    blank3.append(line[27:30])
    x.append(line[30:38]); y.append(line[38:46]); z.append(line[46:54])
    occupancy.append(line[54:60]); bfactor.append(line[60:66])
    blank4.append(line[66:76]); element.append(line[76:78])
    charge.append(line[78:80])
  insert_resnums = []; insert_resnames = []; chains = [];
  recordlist = []
  for i, record in enumerate(record_type):
# determine insertion code
    if insertion_code[i] != ' ' and residue_number[i]+insertion_code[i] not in insert_resnums:
      insert_resnums.append(residue_number[i]+insertion_code[i])
      insert_resnames.append(residue_name[i])
# forget insertion code
      insertion_code[i] = ' '
# determine chain id and record it, if not yet found before
    if chain_id[i] != ' ' and chain_id[i] not in chains:
      chains.append(chain_id[i])
    record = [record_type[i], atom_number[i], blank1[i], atom_name[i], alt_loc_ind[i],
                  residue_name[i], blank2[i], chain_id[i], residue_number[i],
                  insertion_code[i], blank3[i], x[i], y[i], z[i],
                  occupancy[i], bfactor[i], blank4[i], element[i], charge[i]]
# append the accepted record to the overlall record list
    recordlist.append(record)
# report findings so far to the screen
  if chains:
    print "\n----------Chains"
    print "The following (original) chains have been found:"
    for chain in chains:
      print chain
  if insert_resnums:
    print "\n----------Insertions"
    print "The following (original-number) residues were considered as insertions:"
    print "They will be renumbered 'normally' in the final 1-N sequence." 
    for i in range(0,len(insert_resnums)):
      print "%s%s"%(insert_resnames[i],(insert_resnums[i]))

  return(recordlist)

#==================================================
def pdb_write(recordlist, filename):
#==================================================
# uses a record list as created in pdb_read and writes it out to the filename
# using the format below
  pdbout = open(filename, 'w')
  format = "%6s%5s%1s%4s%1s%3s%1s%1s%4s%1s%3s%8s%8s%8s%6s%6s%10s%2s%2s\n"
  for i, record in enumerate(recordlist):
    pdbout.write(format % tuple(record))
  pdbout.close()

#==================================================
def prot_only(recordlist, filename):
#==================================================
# this strips any residues not recogized by Amber libraries...
# in a personalized Amber installation with additional libraries, 
# you might consider extending this list
  global RESPROT
  protlist=[]
  for record in recordlist:
    if record[5] not in RESPROT:
      continue
    else:
      protlist.append(record)
  return(protlist)

#==================================================
def remove_hydrogens(recordlist):
#==================================================
  nohlist = []

  for record in recordlist:
    if record[3][0] == 'H' or record[3][1] == 'H':
      continue
    else:
      nohlist.append(record)

# return the record list with all hydrogens removed
  return(nohlist)

#========================================
def remove_water(recordlist, filename):
#========================================
# removes all water molecules of option -d was specified
  drylist = []; waterlist = []; nwaters = 0
# format = "%6s%5s%1s%4s%1s%3s%1s%1s%4s%1s%3s%8s%8s%8s%6s%6s%10s%2s%2s\n"
# watpdb = open(filename+'_water.pdb', 'w')
  for record in recordlist:
# if oxygen, then count water
    if (record[5] == 'HOH' or record[5] == 'WAT') and 'O' in record[3]:
      nwaters +=1
      waterlist.append(record)
#     watpdb.write(format % tuple(record))
      continue
# if not oxygen, just remove, but do not count, since this is probably hydrogen
    elif  record[5] == 'HOH' or record[5] == 'WAT':
      continue
    else:
      drylist.append(record)

# report the water removal to the screen
  print "\n---------- Water"
  print "%d water molecules have been removed"%nwaters
  print "and stored in the file %s_water.pdb"%filename
# return the dry record list with all water removed
  pdb_write(waterlist, filename+'_water.pdb')
  return(drylist)

#========================================
def remove_altloc(recordlist):
#========================================
  noaltlist = []
  altloc_resnum = []; altloc_resname = []

  for record in recordlist:
# we accept only altlocs 'A' and '1'
    if record[4] != ' ' and record[4] != 'A' and record[4] != '1':
      if record[8] not in altloc_resnum:
        altloc_resnum.append(record[8])
        altloc_resname.append(record[5])
      continue

    else:
      record[4] = ' '
      noaltlist.append(record)

  if altloc_resnum:
    print "\n---------- Alternate Locations (Original Residues!)"
    print "The following residues had alternate locations:"

    for i, rname in enumerate(altloc_resname):
      print "%s_%d"%(rname, int(altloc_resnum[i]))

    print "The alternate coordinates have been discarded."
    print "Only the first occurrence for each atom was kept."

  return(noaltlist)

#==================================================
def atom_wrap(recordlist):
#==================================================
# !!! this function should always be called !!!

# wraps 4-letter hydrogens
  wraplist = []

  for record in recordlist:
    if record[0] != 'ATOM  ' and record[0] != 'HETATM':
      wraplist.append(record)
      continue

# shifts 3-letter atoms if needed
    elif record[3][0] != ' ' and record[3][3] == ' ':
      atomname = record[3][3] + record[3][0:3]
      record[3] = atomname
      wraplist.append(record)
      continue

    else:
      wraplist.append(record)
      continue

  return(wraplist)

#========================================
def renumber(recordlist, filename):
#========================================
  table = open('%s_renum.txt'%filename, 'w')
  renumbered = []; current = -100; iatom = 1
  original = []; oriresname = []; final = []; finresname = []

  for record in recordlist:
    if not 'ATOM' in record[0] and not 'HETATM' in record[0]:
      renumbered.append(record)

    elif current == -100:
      actualnum = record[8]+record[9]
      actualname = record[5]
      if record[3] == ' CA ' or record[3] == 'CH3':
        original.append(record[8]+record[9])
        oriresname.append(record[5])

      record[8] = 1
      record[9] = ' '
      current = 1
      record[1] = iatom
      iatom = iatom + 1

      if record[3] == ' CA ' or record[3] == ' CH3':
        final.append(record[8])
        finresname.append(record[5])

      renumbered.append(record)

    elif record[8]+record[9] == actualnum and record[5] == actualname:

      if record[3] == ' CA ' or record[3] == ' CH3':
        original.append(record[8]+record[9])
        oriresname.append(record[5])

      record[8] = current
      record[9] = ' '
      record[1] = iatom
      iatom = iatom + 1

      if record[3] == ' CA ' or record[3] == ' CH3':
        final.append(record[8])
        finresname.append(record[5])
      renumbered.append(record)

    elif record[8]+record[9] != actualnum or record[5] != actualname :
      actualnum = record[8]+record[9]
      actualname = record[5]
      if record[3] == ' CA ' or record[3] == ' CH3':
        original.append(record[8]+record[9])
        oriresname.append(record[5])

      current = current + 1
      record[8] = current
      record[9] = ' '
      record[1] = iatom
      iatom = iatom + 1

      if record[3] == ' CA ' or record[3] == ' CH3':
        final.append(record[8])
        finresname.append(record[5])
      renumbered.append(record)


  for i in range(0, len(original)):
    table.write("%3s %5s    %3s %5s\n"%(oriresname[i], (original[i]), \
                                    finresname[i], (final[i])))

  return(renumbered)

#========================================
def non_standard(recordlist, filename):
#========================================
# define the common AA and less common AA names that make up proteins
# and that are recognized by Amber routines in ATOM (or HETATM) records
  RES = ('A', 'A3', 'A5', 'ACE', \
         'ALA', 'AN', 'ARG', 'ASH', \
         'ASN', 'ASP', 'Br-', 'C', \
         'C3', 'C5', 'CN', 'CYM', \
         'CYS', 'CYX', 'Cl-', 'Cs+', \
         'DA', 'DA3', 'DA5', 'DAN', \
         'DC', 'DC3', 'DC4', 'DC5', \
         'DCN', 'DG', 'DG3', 'DG5', \
         'DGN', 'DT', 'DT3', 'DT5', \
         'DTN', 'F-', 'G', 'G3', \
         'G5', 'GLH', 'GLN', 'GLU', \
         'GLY', 'GN', 'HID', 'HIE', \
         'HIP', 'HIS', 'HOH', 'HYP', \
         'I-', 'ILE', 'K+', 'LEU', \
         'LYN', 'LYS', 'Li+', 'MET', \
         'Mg+', 'NHE', 'NME', 'Na+', \
         'OHE', 'PHE', 'PL3', 'PRO', \
         'Rb+', 'SER', 'SPC', 'SPF', \
         'SPG', 'T4E', 'THR', 'TP3', \
         'TP4', 'TP5', 'TPF', 'TRP', \
         'TYR', 'U', 'U3', 'U5', \
         'UN', 'VAL', 'WAT', 'U5', \
         'UN', 'VAL', 'WAT')

  hetero = open(filename+'_nonprot.pdb', 'w')
  format = "%6s%5s%1s%4s%1s%3s%1s%1s%4s%1s%3s%8s%8s%8s%6s%6s%10s%2s%2s\n"
  ns_resname = []

  for record in recordlist:
    if record[5].strip() not in RES and record[5] != '   ':
      hetero.write(format % tuple(record))
      if record[5] not in ns_resname:
          ns_resname.append(record[5])

  if ns_resname:
    print "\n---------- Non-Standard Residues"
    print "The following non-standard residue names in the original PDB file"
    print "are not recognized by Amber and have been written to the separate"
    print "file %s_nonprot.pdb"%filename
    print "\n".join(ns_resname)

  return(ns_resname)

#========================================
def non_standard_elbow(recordlist):
#========================================
# define the common AA and less common AA names that make up proteins
# and that are recognized by Amber routines in ATOM (or HETATM) records
  RES = ('A', 'A3', 'A5', 'ACE', \
         'ALA', 'AN', 'ARG', 'ASH', \
         'ASN', 'ASP', 'Br-', 'C', \
         'C3', 'C5', 'CN', 'CYM', \
         'CYS', 'CYX', 'Cl-', 'Cs+', \
         'DA', 'DA3', 'DA5', 'DAN', \
         'DC', 'DC3', 'DC4', 'DC5', \
         'DCN', 'DG', 'DG3', 'DG5', \
         'DGN', 'DT', 'DT3', 'DT5', \
         'DTN', 'F-', 'G', 'G3', \
         'G5', 'GLH', 'GLN', 'GLU', \
         'GLY', 'GN', 'HID', 'HIE', \
         'HIP', 'HIS', 'HOH', 'HYP', \
         'I-', 'ILE', 'K+', 'LEU', \
         'LYN', 'LYS', 'Li+', 'MET', \
         'Mg+', 'NHE', 'NME', 'Na+', \
         'OHE', 'PHE', 'PL3', 'PRO', \
         'Rb+', 'SER', 'SPC', 'SPF', \
         'SPG', 'T4E', 'THR', 'TP3', \
         'TP4', 'TP5', 'TPF', 'TRP', \
         'TYR', 'U', 'U3', 'U5', \
         'UN', 'VAL', 'WAT', 'U5', \
         'UN', 'VAL', 'WAT')

  
  format = "%6s%5s%1s%4s%1s%3s%1s%1s%4s%1s%3s%8s%8s%8s%6s%6s%10s%2s%2s\n"
  ns_resname = []

  for record in recordlist:
    if record[5].strip() not in RES:
      if record[5].strip() not in ns_resname:
        ns_resname.append(record[5].strip())
        try: f.close()
        except: pass  
        f=open('4antechamber_%s.pdb' %(record[5].strip()), 'w')
        resid=record[8]
        f.write(format % tuple(record))
      else: 
        if record[8]==resid:
          f.write(format % tuple(record))
        else:
          if not f.closed: f.close()  

  return ns_resname


#========================================
def find_his(recordlist):
#========================================
  nhis = 0; hisresname = []; hisresnum = []
  histidines = ('HIS', 'HID', 'HIP', 'HIN', 'HIE')
  for record in recordlist:
    if record[5] in histidines and record[3] == ' CA ':
      hisresname.append(record[5])
      hisresnum.append(record[8])
      nhis += 1

  if nhis > 0:
    print "\n---------- Histidines (Renumbered Residues!)"
    print "The following %d histidines are found in the PDB file: "%nhis

    for i, rname in enumerate(hisresname):
      print '%s_%d' % (rname, int(hisresnum[i]))

    print "If HIS, Amber will consider them as HIE (epsilon-HIS) by default."
    print "You might need to check their tautomerism or protonation state"
    print "and change them to HID (delta-HIS) or HIP (protonated HIS)"

  return(recordlist)

#========================================
def find_disulfide(recordlist, filename):
#========================================
  cys_residues = [];  cys_sgx = []; cys_sgy = []; cys_sgz = []
  cyx_residues = []; ncys = 0; ncyx = 0

  print "\n---------- Cysteines in Disulfide Bonds (Renumbered Residues!)"
  for record in recordlist:

    if 'SG' in record[3] and ('CYS' in record[5] or 'CYX' in record[5]):
      cys_residues.append(record[8])
      cys_sgx.append(record[11])
      cys_sgy.append(record[12])
      cys_sgz.append(record[13])
      ncys += 1

  if ncys > 0:
    
    sslink = open('%s_sslink'%filename, 'w')
    dist = [[0 for i in range(ncys)] for j in range(ncys)]

    for i in range(0, ncys-1):
      for j in range(i+1, ncys):
        dx = float(cys_sgx[i]) - float(cys_sgx[j])
        dx2 = dx*dx
        dy = float(cys_sgy[i]) - float(cys_sgy[j])
        dy2 = dy*dy
        dz = float(cys_sgz[i]) - float(cys_sgz[j])
        dz2 = dz*dz
        dist[i][j] = sqrt(dx2 +dy2 +dz2)
        if dist[i][j] < 2.5 and dist[i][j] > 0.1:
          cyx_residues.append(cys_residues[i])
          cyx_residues.append(cys_residues[j])
          print("CYS_%s - CYS_%s: S-S distance = %f Ang."%(cys_residues[i], cys_residues[j], dist[i][j]))
          sslink.write('%s %s\n'%(cys_residues[i], cys_residues[j]))
          ncyx += 1

# rename the CYS to CYX for disulfide-involved cysteines
    for record in recordlist:
      if record[8] in cyx_residues:
        record[5] = 'CYX'
      else:
        continue
  if ncyx:
    print "The above CYS have been renamed to CYX in the new PDB file."
    print "The created file '%s_sslink' can be used with --disul in 'pytleap'"%filename
    print "to generate the correct parameter-topology file."

  else:
    print "No disulfide bonds have been detected."
  return(recordlist)

#========================================
def find_gaps(recordlist):
#========================================
  global RESPROT
  ca_atoms = []; ca_resnum = []; ca_resname = [];
  ca_x = []; ca_y = []; ca_z = []
  gaplist = [];

  for record in recordlist:
    if ('CA' in record[3] or 'CH3' in record[3]) and record[5] in RESPROT:
      ca_atoms.append(record[1])
      ca_resnum.append(record[8])
      ca_resname.append(record[5])
      ca_x.append(float(record[11]))
      ca_y.append(float(record[12]))
      ca_z.append(float(record[13]))

  nca = len(ca_atoms)
  ngaps = 0

  for i in range(nca-1):
    dx = ca_x[i] - ca_x[i+1]
    dx2 = dx*dx
    dy = ca_y[i] - ca_y[i+1]
    dy2 = dy*dy
    dz = ca_z[i] - ca_z[i+1]
    dz2 = dz*dz
    gap = sqrt(dx2 +dy2 +dz2)

    if gap > 5.0:
      gaprecord = (gap, ca_resname[i], int(ca_resnum[i]), ca_resname[i+1], int(ca_resnum[i+1]))
      gaplist.append(gaprecord)
      ngaps += 1

  if ngaps > 0:
    print "\n---------- Gaps (Renumbered Residues!)"
    format = "gap of %lf A between %s_%d and %s_%d"

    for i, gaprecord in enumerate(gaplist):
      print (format % tuple(gaprecord))

    print "You MUST (!!!) insert a TER record between the residues listed above and"
    print "consider to introduce caps (ACE and NME) at the dangling N- and C-terminals."

  return()

#========================================
def find_incomplete(recordlist):
#========================================
# finds residues with missing heavy atoms in the following list of residues;
# dictionary with number of heavy atoms:
  #paj complete this list
  HEAVY = {'ALA':5,  'ARG':11, 'ASN':8,  'ASP':8, \
           'CYS':6,  'GLN':9,  'GLU':9,  'GLY':4, \
           'HIS':10, 'ILE':8,  'LEU':8,  'LYS':9, \
           'MET':8,  'PHE':11, 'PRO':7,  'SER':6, \
           'THR':7,  'TRP':14, 'TYR':12, 'VAL':7, \
           'HID':10, 'HIE':10, 'HIN':10, 'HIP':10, \
           'CYX':6,  'ASH':8,  'GLH':9,  'LYH':9}
  print '\n---------- Missing Heavy Atoms (Renumbered Residues!)' 
  resnum = []
  resname = []
  resheavy = []
  flag = 0
  nheavy = 1
  length = len(recordlist)
  for i, record in enumerate(recordlist):
    if i == length-1:
      for j, res in enumerate(resnum):
        missing = HEAVY[resname[j]] - resheavy[j]
        if missing > 0:
          flag = 1
          print "%s_%s misses %d heavy atom(s)"%(resname[j], res, missing) 
      if flag == 0:
        print "None"
      return()
    if not HEAVY.has_key(record[5]):
      continue
    if recordlist[i+1][8] == record[8]:
      nheavy += 1
    else:
      resnum.append(record[8])
      resname.append(record[5])
      resheavy.append(nheavy)
      nheavy = 1
      continue
  return()

def run(arg_pdbout, arg_pdbin, arg_nohyd, arg_dry, arg_prot, arg_elbow=False):
  filename, extension = os.path.splitext(arg_pdbout)
  pdbin = arg_pdbin
  recordlist = pdb_read(pdbin)
  
  # wrap all atom names to pure standard (always):======================
  recordlist = atom_wrap(recordlist)
  # remove alternate locations and keep only the first one:=============
  recordlist = remove_altloc(recordlist)
  # remove hydrogens if option -y is used:==============================
  if arg_nohyd:
    recordlist = remove_hydrogens(recordlist)
  # find non-standard Amber residues:===================================
  non_standard(recordlist, filename)
  if arg_elbow:
    ns_names=non_standard_elbow(recordlist)
  # keep only protein:==================================================
  if arg_prot:
    recordlist = prot_only(recordlist, filename)
  # remove water if -d option used:=====================================
  if arg_dry:
    recordlist = remove_water(recordlist, filename)
  # renumber atoms and residues:========================================
  recordlist = renumber(recordlist, filename)
  # after this call, residue numbers refer to the ***new*** PDB file
  #=====================================================================
  # find histidines that might have to be changed
  recordlist = find_his(recordlist)
  # find possible S-S in the final protein:=============================
  recordlist = find_disulfide(recordlist, filename)
  # find possible gaps:==================================================
  find_gaps(recordlist)
  # count heavy atoms
  find_incomplete(recordlist)
  # =====================================================================
  # make final output to new PDB file
  pdb_write(recordlist, arg_pdbout)
  print ""
  try: ns_names
  except: return
  else: return ns_names
    
#========================================main===========================
if __name__ ==  "__main__":
  parser = OptionParser()
  parser.add_option("-i","--in", metavar = "FILE", dest = "pdbin",
                    help = "PDB input file                      (no default)")
  parser.add_option("-o","--out", metavar = "FILE", dest = "pdbout",
                    help = "PDB output file                     (no default)")
  parser.add_option("-y","--nohyd", action = "store_true", dest = "nohyd",
                    help = "remove all hydrogen atoms           (default: no)")
  parser.add_option("-d","--dry", action = "store_true", dest = "dry",
                    help = "remove all water molecules          (default: no)")
  parser.add_option("-p", "--prot", action = "store_true", dest = "prot",
                    help = "keep only Amber-compatible residues (default: no)")

  if len(sys.argv) == 1:
    print("---------------------------------------------------------")
    print(" pdb4amber version 0.5 (December 2013)")
    print("---------------------------------------------------------")
    parser.print_help()
    sys.exit(-1)
  else:pass

  (opt, args) = parser.parse_args()

  if not opt.pdbin or not opt.pdbout:
    print "You must specify an input file and output file"
    print "using the option -i (or --in) and -o (or --out)"
    sys.exit(-1)

  if opt.pdbin == opt.pdbout:
    print "The input and outout file names cannot be the same!"

  run(opt.pdbout, opt.pdbin, opt.nohyd, opt.dry, opt.prot)


