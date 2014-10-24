# LIBTBX_SET_DISPATCHER_NAME phenix.AmberPrep

import os, sys
import iotbx.pdb
from libtbx import phil
import libtbx.phil.command_line
from iotbx import pdb
from amber_adaptbx import pdb4amber
from elbow.command_line import builder
from libtbx import easy_run
import StringIO
from amber_adaptbx import fix_ambpdb
from amber_adaptbx import amber_library_server
from libtbx.utils import Sorry

master_phil_string = """
amber_prep
  .caption = Prepare Amber files
{
  input
  {
    pdb_file_name = None
      .type = path
    nproc = 1
      .type = int
      .short_caption = Number processes to use
  }
  actions
  {
    minimise = amber_all amber_h phenix_all *off
      .type = choice
    clean = on *off
      .type = choice
    redq = False
      .type = bool
  }
  output
  {
    file_name = None
      .type = path
  }
}
"""

master_params = master_phil_string # need for auto documentation
master_phil = phil.parse(master_phil_string,
                         process_includes=True,
                         )
def setup_parser():
  from libtbx.option_parser import OptionParser
  usage="""
  phenix.AmberPrep 3a37.pdb minimise=phenix_all
"""
  parser = OptionParser(
    prog="phenix.amber_prep",
    version="""
  up-to-date version
""",
    usage=usage,
    )
  # Input options
  parser.add_option("",
                    "--show_defaults",
                    dest="show_defaults",
                    default=False,
                    action="store_true",
                    help="Display defaults",
                    )
  return parser

def setup_options_args(rargs):
  rargs = list(rargs)
  parser = setup_parser()
  (options, args) = parser.parse_args(args=rargs)
  if options.show_defaults:
    tmp = phil.parse(master_phil_string,
                     process_includes=True,
                     )
    tmp.show()
    sys.exit()
  if len(args)==0:
    parser.print_help()
    sys.exit()

  #
  argument_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_phil,
    home_scope="amber_prep")
  #
  pdbs = []
  phils = []
  phil_args = []
  for arg in args:
    if os.path.isfile(arg) :
      if iotbx.pdb.is_pdb_file(arg):
        pdbs.append(arg)
      else:
        try : 
          file_phil = phil.parse(file_name=arg)
        except RuntimeError :
          pass
        else :
          phils.append(file_phil)
    else :
      phil_args.append(arg)
      phils.append(argument_interpreter.process(arg))
  working_phil = master_phil.fetch(sources=phils)
  assert len(pdbs)==1
  working_phil.show()
  working_params = working_phil.extract()
  working_params.amber_prep.input.pdb_file_name = pdbs[0]
  #check_working_params(working_params)
  preamble = 'test' #get_output_preamble(working_params)
  preamble = preamble.replace(".pdb","")
  print "  Writing effective parameters to %s.eff\n" % preamble
  #working_phil.format(python_object=working_params).show()
  print "#phil __ON__"
  master_phil.fetch_diff(source=master_phil.format(
      python_object=working_params)).show()
  print "#phil __OFF__\n\n"
  f=file("%s.eff" % preamble, "wb")
  f.write(working_phil.format(python_object=working_params).as_str())
  f.close()
  return working_params

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

def validatePdb(pdb_hierarchy):
  #pdb_hierarchy.show()
  from mmtbx import conformation_dependent_library
  gaps=[]
  for three in conformation_dependent_library.generate_protein_threes(
      hierarchy=pdb_hierarchy,
      geometry=None,
      include_non_linked=True,
      ):
    if not three.are_linked():
      print three.show()
      gaps.append(three)
  if gaps:
    return gaps
  return False

def get_molecule_from_hierarchy(hierarchy, resname):
  # only works for non altloc files
  from elbow.chemistry.MoleculeClass import MoleculeClass
  mol = MoleculeClass()
  for residue_group in hierarchy.residue_groups():
    for atom_group in residue_group.atom_groups():
      if atom_group.resname==resname:
        for atom in atom_group.atoms():
          mol.AddAtom(atom.element, xyz=atom.xyz)
          mol[-1].name = atom.name
        break
  return mol
  
# run pdb4amber
def run_pdb4amber(pdb_filename):
  ns_names=pdb4amber.run('4tleap.pdb', pdb_filename, arg_elbow=True)
  return ns_names
  
# run elbow and antechamber
def _run_elbow_antechamber(pdb_hierarchy, residue_name, debug=False):
  debug=1
  pdb_mol = get_molecule_from_hierarchy(pdb_hierarchy, residue_name)
  names = []
  for atom in pdb_mol:
    names.append(atom.name.strip())
  pdb_set = set(names)
  print >> sys.stderr, "\n=================================================="
  print >> sys.stderr, "Running elbow/antechamber for %s " %residue_name
  print >> sys.stderr, "=================================================="
  if debug:
    import pickle
    pf = "%s.pickle" % residue_name
    if os.path.exists(pf):
      f=file(pf, "rb")
      mol=pickle.load(f)
      f.close()
    else:
      mol = builder.run(chemical_component=residue_name,
                        no_output=True,
                        silent=True)
      mol.WritePickle()
  else:
    # this needs to be changed for novel ligands
    mol = builder.run(chemical_component=residue_name,
                      no_output=True,
                      silent=True)
  names = []
  for atom in mol:
    names.append(atom.name.strip())
  cc_set_h = set(names)
  names = []
  for atom in mol:
    if atom.isH(): continue
    names.append(atom.name.strip())
  cc_set = set(names)
  if ( pdb_set.symmetric_difference(cc_set) and
       pdb_set.symmetric_difference(cc_set_h)
       ):
    from elbow.utilities import molecule_superposition
    pdb_mol.Bondise(add_bond_order=False)
    rc = molecule_superposition.run(mol,
                                    pdb_mol,
                                    use_hydrogens=False,
                                    seed_with_unique_atoms=True,
      )
    #molecule_superposition.print_return_list(rc)
    for atom1, atom2 in rc:
      atom1.name = atom2.name

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
  cmd='antechamber -i sqm.pdb -fi pdb -o %s.mol2 -fo mol2 \
      -nc %s -m %d -s 2 -pf y -c bcc -at gaff' \
      %(residue_name, mol.charge, mol.multiplicity)    
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

def run_elbow_antechamber(pdb_hierarchy, ns_names, nproc=1, debug=False):
  if nproc>1:
    print "\n\tParallel processes not implemented\n"
  for residue_name in ns_names:
    if amber_library_server.is_in_components_lib(residue_name):
      print "%s already in amber monomer library. Skipping elbow/antechamber run for this residue.\n" %residue_name
      continue
    else:
      _run_elbow_antechamber(pdb_hierarchy, residue_name, debug=debug)
  return 0

# prepare tleap input (set box) or run pytleap?
def run_tleap(input_pdb, output_base,ns_names, reorder_residues, logfile, redq=False):
  f=open('tleap.in','w')
  if redq:
    f.write('source leaprc.ff14SB.redq\n')
  else:
    f.write('source leaprc.ff14SB\n')
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
  f.write('x = loadpdb %s\n' %input_pdb)
  f.write('set x box {20.000   20.000   20.000}\n')
  f.write('set default nocenter on\n')
  f.write('set default reorder_residues %s\n' %reorder_residues)
  f.write('saveAmberParm x %s.prmtop %s.rst7\n' %(output_base, output_base))
  f.write('quit\n')   
  f.close()
  cmd = 'tleap -f tleap.in >%s' %logfile
  print cmd
  ero=easy_run.fully_buffered(cmd)
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
  print cmd
  ero=easy_run.fully_buffered(cmd)
  return 0

# make pdb 
def run_ambpdb(base):
  cmd='ambpdb -p  %s.prmtop <%s.rst7 >new.pdb' %(base, base)
  print cmd
  ero=easy_run.fully_buffered(cmd)
  ero.show_stdout()
  ero.show_stderr()
  return 0

#add cryst1 and sscale
def finalizePdb(pdb_filename,cryst,base):
  pdb_inp = pdb.input(pdb_filename)
  pdb_hierarchy=pdb_inp.construct_hierarchy()
  pdbstring=pdb_hierarchy.as_pdb_string(crystal_symmetry=cryst)
  f=open('4phenix_'+base+'.pdb','w')
  f.write(pdbstring)
  f.close()
  #~ import code; code.interact(local=locals())
  return 0

def run_UnitCell(input_file,output_file):
  cmd='UnitCell -p %s -o %s' %(input_file, output_file)
  print cmd
  ero=easy_run.fully_buffered(cmd)
  ero.show_stdout()
  ero.show_stderr()

def uc(pdb_filename,ns_names,cryst1, base, redq=False):
  #add SMTRY/CRYST1 to 4tleap.pdb -> 4UnitCell.pdb
  with open("4UnitCell.pdb","wb") as fout:
    with open(pdb_filename) as fin:
      lines = fin.readlines()
      smtry = [line for line in lines if "SMTRY" in line]
      smtry = ''.join(smtry)
      if not smtry:
        print '"%s"' % smtry
        raise Sorry("REMARK 290 SMTRY1,2,3 records required")
      fout.write(smtry)
      # import code; code.interact(local=locals())
      cryst1card = [line for line in lines if "CRYST1" in line]
      if len(cryst1card) <1:
        raise Sorry("CRYST1 record required")
      fout.write(cryst1card[0])
    with open("new.pdb") as fin:
      for line in fin:
        fout.write(line)
  run_UnitCell('4UnitCell.pdb', '4tleap_uc.pdb')
  run_tleap('4tleap_uc.pdb', 'uc', ns_names, reorder_residues='off', logfile='tleap_uc.log', redq=redq)
  run_ChBox('uc', cryst1)
  os.rename('uc.rst7', '4amber_%s.rst7' % base)
  os.rename('uc.prmtop', '4amber_%s.prmtop' % base)

def run_minimise(base, cryst1, option=None):
  if option is None: return

  if option=="amber_h":
    f=open('min_H.in', 'w')
    f.write("Initial minimization\n")
    f.write("&cntrl\n")
    f.write(" ntwx   = 0, ntb    = 1, cut    = 9.0,     nsnb   = 10,\n")
    f.write(" ntr    = 1, restraint_wt = 50.0, restraintmask ='!@H=',\n")
    f.write(" imin   = 1, maxcyc =1000, ncyc   = 200, ntmin  = 1, \n")
    f.write("&end\n")
    f.close()
    cmd='sander -O -i min_H.in -p asu.prmtop -c asu.rst7 -o min_H.out \
         -ref asu.rst7 -r asu_minH.rst7'
    print cmd
    ero=easy_run.fully_buffered(cmd)
    assert (ero.return_code == 0)
    run_ChBox('asu_minH',cryst1)
    cmd='ambpdb -p asu.prmtop <asu_minH.rst7 >new.pdb'
    print cmd
    ero=easy_run.fully_buffered(cmd)
    assert (ero.return_code == 0)
    ero.show_stdout()
    ero.show_stderr()
    fix_ambpdb.run('4tleap.pdb', 'new.pdb', 'new2.pdb' )
    finalizePdb('new2.pdb',cryst1, base)

  elif option=="amber_all":  
    f=open('min_all.in', 'w')
    f.write("Initial minimization\n")
    f.write("&cntrl\n")
    f.write(" ntwx   = 0, ntb    = 1, cut    = 9.0,     nsnb   = 10,\n")
    f.write(" imin   = 1, maxcyc = 1000, ncyc   = 200, ntmin  = 1, \n")
    f.write("&end\n")
    f.close()
    cmd='sander -O -i min_all.in -p asu.prmtop -c asu.rst7 -o min_H.out \
         -ref asu.rst7 -r asu_minall.rst7'
    print cmd
    ero=easy_run.fully_buffered(cmd)  
    assert (ero.return_code == 0)
    run_ChBox('asu_minall',cryst1)
    cmd='ambpdb -p asu.prmtop <asu_minall.rst7 >new.pdb'
    print cmd
    ero=easy_run.fully_buffered(cmd)
    assert (ero.return_code == 0)
    ero.show_stdout()
    ero.show_stderr()
    fix_ambpdb.run('4tleap.pdb', 'new.pdb', 'new2.pdb' )
    finalizePdb('new2.pdb',cryst1, base+'_minall')  

  elif option=="phenix_all":
    cmd='phenix.geometry_minimization 4phenix_%s.pdb amber.use_amber=True \
         amber.topology_file_name=4amber_%s.prmtop \
         amber.coordinate_file_name=4amber_%s.rst7  \
         output_file_name_prefix=4phenix_%s_minPhenix ' %(base,base,base,base)
    restraints = "%s.ligands.cif" % base
    if os.path.exists(restraints):
      cmd += " %s" % restraints
    print cmd
    ero=easy_run.fully_buffered(cmd).raise_if_errors() #.return_code
    assert (ero.return_code == 0)
    ero.show_stdout()
    ero.show_stderr()
    os.rename('4phenix_%s_minPhenix.pdb' %base, '4phenix_%s.pdb' % base)
  return 0

def run_clean():
  files_to_clean = """
    4tleap_nonprot.pdb
    4tleap.pdb
    4tleap_renum.txt
    4tleap_sslink
    4tleap_uc.pdb
    4UnitCell.pdb
    asu.prmtop
    asu.rst7
    init.pdb
    leap.log
    new2.pdb
    new.pdb
    test.eff
    tleap_asu.log
    tleap.in
    tleap_uc.log
    uc.prmtop
    uc.rst7
    sqm.pdb
    sqm.out
    sqm.in
    4antechamber*
  """

  import glob
  for filename in files_to_clean.strip().split():
    for file in glob.glob(filename):
      os.remove(file)

def run(rargs):
  working_params = setup_options_args(rargs)
  inputs = working_params.amber_prep.input
  actions = working_params.amber_prep.actions
  base = os.path.basename(inputs.pdb_file_name).split('.')[0]
  pdb_inp = pdb.input(inputs.pdb_file_name)
  pdb_hierarchy=pdb_inp.construct_hierarchy()
  invalid = validatePdb(pdb_hierarchy)
  if invalid: raise Sorry( 'PDB input is not "valid"' )
  cryst1=initializePdb(inputs.pdb_file_name)
  ns_names=run_pdb4amber('init.pdb')
  #
  run_elbow_antechamber(pdb_hierarchy, ns_names, nproc=inputs.nproc)
  #
  print >> sys.stderr, "\n=================================================="
  print >> sys.stderr, "Preparing asu files and 4phenix_%s.pdb" %base
  print >> sys.stderr, "=================================================="
  run_tleap('4tleap.pdb','asu',ns_names,reorder_residues='on', logfile='tleap_asu.log', redq=actions.redq)
  run_ChBox('asu',cryst1)
  run_ambpdb('asu')
  fix_ambpdb.run('4tleap.pdb', 'new.pdb', 'new2.pdb' )
  finalizePdb('new2.pdb', cryst1, base)
  print >> sys.stderr, "\n=================================================="
  print >> sys.stderr, "Preparing uc files: %s.prmtop and %s.rst7" %(base,base)
  print >> sys.stderr, "=================================================="
  uc(inputs.pdb_file_name, ns_names,cryst1, base, redq=actions.redq)
  if actions.minimise == "off":
    pass
  else:
    print >> sys.stderr, "\n=================================================="
    print >> sys.stderr, "Minimizing input coordinates."
    print >> sys.stderr, "=================================================="
    run_minimise(base, cryst1, actions.minimise)
  if actions.clean == "on":
    run_clean()
  return '4phenix_%s.pdb' % base

if __name__=="__main__":
  if 1:
    args = sys.argv[1:]
    del sys.argv[1:]
    run(args)
  else:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_file_name", help="name of pdb file")
    parser.add_argument("min", help="option to minimize", default=0)
    parser.add_argument("-c","--no_clean", help="don't remove "
                        "intermediate files", action='True', default=False )
    args = parser.parse_args()
    run(args.pdb_file_name, minimize=args.min, clean=args.no_clean)
