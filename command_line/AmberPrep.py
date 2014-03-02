# LIBTBX_SET_DISPATCHER_NAME phenix.AmberPrep

import os, sys
from iotbx import pdb
from amber_adaptbx import pdb4amber
from elbow.command_line import builder
from libtbx import easy_run
import StringIO
from amber_adaptbx import fix_ambpdb
from amber_adaptbx import amber_library_server

# get box info & space group info
def initializePdb(pdb_filename):
  pdb_inp = pdb.input(pdb_filename)
  cryst1=pdb_inp.crystal_symmetry_from_cryst1()
  pdb_hierarchy=pdb_inp.construct_hierarchy()
  pdbstring=pdb_hierarchy.as_pdb_string()
  f=open('init.pdb','w')
  f.write(pdbstring)
  f.close()
  return cryst1
  
# run pdb4amber
def run_pdb4amber(pdb_filename):
  ns_names=pdb4amber.run('4tleap.pdb', pdb_filename, arg_elbow=True)
  return ns_names
  
# run elbow and antechamber
def run_elbow_antechamber(ns_names):
  for residue_name in ns_names:
    #~ if amber_library_server.is_in_components_lib(residue_name):
    if False:
      print "%s already in amber monomer library. Skipping elbow/antechamber run for this residue.\n" %residue_name
      continue
    else:
      mol = builder.run(chemical_component=residue_name,
                      no_output=True, silent=True)
      mol.WritePDB('4antechamber_%s.pdb' %residue_name)                  
      mol.Multiplicitise()
      print mol.DisplayBrief()
      cmd='antechamber -i 4antechamber_%s.pdb -fi pdb -o %s.mol2 -fo mol2 \
          -nc %d -m %d -s 2 -pf y -c bcc -at gaff' \
          %(residue_name, residue_name, mol.charge, mol.multiplicity)
      print cmd
      ero=easy_run.fully_buffered(cmd)
      stdo = StringIO.StringIO()
      ero.show_stdout(out=stdo)
      for line in stdo.getvalue().splitlines():
        if line.find('APS')>-1:
          print line
        if line.find('Error')>-1:
          print line  
      cmd='parmchk2 -i %s.mol2 -f mol2 -o %s.frcmod' %(residue_name, residue_name)
      print cmd
      easy_run.fully_buffered(cmd)
  return 0

# prepare tleap input (set box) or run pytleap?
def run_tleap(pdb_filename,ns_names):
  f=open('tleap.in','w')
  f.write('source leaprc.ff12SB\n')
  f.write('source leaprc.gaff\n')
  #f.write('loadoff atomic_ions.lib\n')
  f.write('loadamberparams frcmod.ionslrcm_iod\n')
  f.write('loadamberparams frcmod.ionsjc_spce\n')
  f.write('WAT = SPC\n')
  f.write('HOH = SPC\n')
  f.write('loadAmberParams frcmod.spce\n')
  for res in ns_names:
    if amber_library_server.is_in_components_lib(res):
        res_path=amber_library_server.path_in_components_lib(res)
        f.write('%s = loadmol2 %s\n' %(res,res_path[1]))
        f.write('loadAmberParams %s\n' %(res_path[0]))
    else:
      f.write('%s = loadmol2 %s.mol2\n' %(res,res))
      f.write('loadAmberParams %s.frcmod\n' %res)
  f.write('x = loadpdb 4tleap.pdb\n')
  #~ if os.path.isfile('4tleap_sslink'):
    #~ input = open('4tleap_sslink', 'r')
    #~ for line in input.readlines():
      #~ cys1 = line.split()[0]
      #~ cys2 = line.split()[1]
      #~ f.write("bond x.%s.SG x.%s.SG\n" %(cys1, cys2))
  f.write('set x box {20.000   20.000   20.000}\n')
  f.write('set default nocenter on\n')
  f.write('saveAmberParm x %s.prmtop %s.rst7\n' %(pdb_filename, pdb_filename))
  f.write('quit\n')   
  f.close()
  ero=easy_run.fully_buffered('tleap -f tleap.in')
  ero.show_stdout()
  ero.show_stderr()
  return 0

#change box size in rst7 file
def run_ChBox(base,cryst1):
  uc = cryst1.unit_cell().parameters()
  cmd="ChBox -c %s.rst7 -o chbox.rst7" % base
  cmd += " -X %s -Y %s -Z %s -al %s -bt %s -gm %s " % tuple(uc)
  print cmd
  ero=easy_run.fully_buffered(cmd)
  ero.show_stdout()
  ero.show_stderr()
  cmd='mv chbox.rst7 %s.rst7' %base
  ero=easy_run.fully_buffered(cmd)
  return 0

# make pdb 
def run_ambpdb(base):
  cmd='ambpdb -p %s.prmtop <%s.rst7 >new.pdb' %(base, base)
  print cmd
  ero=easy_run.fully_buffered(cmd)
  ero.show_stdout()
  ero.show_stderr()
  return 0

#add cryst1 and smtry
def finalizePdb(pdb_filename,cryst,base):
  pdb_inp = pdb.input(pdb_filename)
  pdb_hierarchy=pdb_inp.construct_hierarchy()
  pdbstring=pdb_hierarchy.as_pdb_string(crystal_symmetry=cryst)
  f=open('4phenix_'+base+'.pdb','w')
  f.write(pdbstring)
  f.close()
  #~ import code; code.interact(local=locals())
  return 0

#fix residue names for phenix, add original Bfactors
def run(pdb_filename):
  base = os.path.basename(pdb_filename).split('.')[0]
  cryst1=initializePdb(pdb_filename)
  ns_names=run_pdb4amber('init.pdb')
  run_elbow_antechamber(ns_names)
  run_tleap(base,ns_names)
  run_ChBox(base,cryst1)
  run_ambpdb(base)
  fix_ambpdb.run('4tleap.pdb', 'new.pdb', 'new2.pdb' )
  finalizePdb('new2.pdb',cryst1, base)
  
  
if __name__=="__main__":
  run(sys.argv[1])
