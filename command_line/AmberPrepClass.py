# LIBTBX_SET_DISPATCHER_NAME phenix.AmberPrep

import os, sys
import shutil
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
      clean = *on off
        .type = choice
      redq = False
        .type = bool
      save_cpp_traj_prmtop = True
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

def get_output_preamble(params):
  preamble = params.amber_prep.input.pdb_file_name
  if params.amber_prep.output.file_name:
    preamble = params.amber_prep.output.file_name
  preamble = os.path.basename(preamble)
  preamble = preamble.replace(".pdb","")
  return preamble

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
  preamble = get_output_preamble(working_params)
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

def check_required_output_filenames(filenames=None):
  # only tleap currently
  if filenames is None: return
  print 'Checking output filenames'
  for filename in filenames:
    print '  file : %s' % filename
    if not os.path.exists(filename):
      raise Sorry("  Output file not present : %s" % filename)
    if os.stat(filename).st_size==0:
      raise Sorry("  Output file is empty : %s" % filename)

def print_cmd(cmd):
  print "\n~> %s\n" % cmd

class amber_prep_run_class:
  def __init__(self, base_name):
    self.base = base_name

  def __repr__(self):
    outl = "AmberPrepRunner"
    outl += "\n  Base : %s" % self.base
    return outl
    
    
  # get box info & space group info
  def initializePdb(self, pdb_filename):
    self.pdb_filename  = pdb_filename
    self.pdb_inp       = pdb.input(pdb_filename)
    self.cryst1        = self.pdb_inp.crystal_symmetry_from_cryst1()
    self.pdb_hierarchy = self.pdb_inp.construct_hierarchy()

  def validatePdb(self):
    assert self.pdb_hierarchy
    from mmtbx import conformation_dependent_library
    gaps=[]
    for three in conformation_dependent_library.generate_protein_threes(
        hierarchy=self.pdb_hierarchy,
        geometry=None,
        include_non_linked=True,
        ):
      if not three.are_linked():
        print three.show()
        gaps.append(three)
    if gaps:
      return gaps
    return False

  # run pdb4amber
  def run_pdb4amber(self):
    assert self.pdb_hierarchy
    pdbstring=self.pdb_hierarchy.as_pdb_string()
    init = "%s_init.pdb" % self.base
    f=open(init,'w')
    f.write(pdbstring)
    f.close()
    self.tleap_input_pdb = "%s_4tleap.pdb" % self.base
    self.ns_names=pdb4amber.run(self.tleap_input_pdb,
                                self.pdb_filename,
                                arg_elbow=True,
      )
    # init is identical to old method
    if os.path.exists(init): os.remove(init)

  def run_elbow_antechamber(self, nproc=1, debug=False):
    assert self.pdb_hierarchy
    assert hasattr(self, "ns_names")
    if nproc>1:
      print "\n\tParallel processes not implemented\n"
    for residue_name in self.ns_names:
      if amber_library_server.is_in_components_lib(residue_name):
        print "%s already in amber monomer library. Skipping elbow/antechamber run for this residue.\n" %residue_name
        continue
      else:
        _run_elbow_antechamber(self.pdb_hierarchy, residue_name, debug=debug)
    return 0

  # prepare tleap input (set box) or run pytleap?
  def run_tleap(self,
                input_pdb,
                output_base,
                #ns_names,
                reorder_residues,
                #logfile,
                redq=False,
                verbose=False,
                ):
    assert self.tleap_input_pdb
    assert hasattr(self, "ns_names")
    tleap_input_file = "%s_%s_tleap_input_run" % (self.base, output_base)
    f=file(tleap_input_file, "wb")
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
    for res in self.ns_names:
      if amber_library_server.is_in_components_lib(res):
        res_path=amber_library_server.path_in_components_lib(res)
        f.write('%s = loadmol2 %s\n' %(res,res_path[1]))
        f.write('loadAmberParams %s\n' %(res_path[0]))
      else:
        f.write('%s = loadmol2 %s.mol2\n' %(res,res))
        f.write('loadAmberParams %s.frcmod\n' %res)
    #
    # input PDB file
    #
    f.write('x = loadpdb %s\n' % input_pdb)
    f.write('set x box {20.000   20.000   20.000}\n')
    f.write('set default nocenter on\n')
    f.write('set default reorder_residues %s\n' % reorder_residues)
    f.write('saveAmberParm x %s.prmtop %s.rst7\n' %(
      "%s_%s" % (self.base, output_base),
      "%s_%s" % (self.base, output_base),
      )
      )
    f.write('quit\n')   
    f.close()
    cmd = 'tleap -f %s' % tleap_input_file
    print_cmd(cmd)
    ero=easy_run.fully_buffered(cmd)
    if verbose:
      ero.show_stdout()
      ero.show_stderr()
    # check output
    check_required_output_filenames(["%s_%s.prmtop" % (self.base, output_base),
                                     "%s_%s.rst7" % (self.base, output_base),
                                     ])
    # remove files
    #if os.path.exists(tleap_input_file): os.remove(tleap_input_file)
    return 0

  #change box size in rst7 file
  def run_ChBox(self, output_base):
    assert self.cryst1
    uc = self.cryst1.unit_cell().parameters()
    cmd="ChBox -c %s_%s.rst7 -o %s_%s.rst7" % (self.base,
                                               output_base,
                                               self.base,
                                               output_base,
                                               )
    cmd += " -X %s -Y %s -Z %s -al %s -bt %s -gm %s " % tuple(uc)
    print_cmd(cmd)
    ero=easy_run.fully_buffered(cmd)
    ero.show_stdout()
    ero.show_stderr()
    return 0

  # make pdb 
  def run_ambpdb(self, save_cpp_traj_prmtop=False):
    assert self.base
    cmd='ambpdb -p %s_asu.prmtop < %s_asu.rst7 > %s_new.pdb' % tuple(
      [self.base]*3
      )
    print_cmd(cmd)
    ero=easy_run.fully_buffered(cmd)
    ero.show_stdout()
    ero.show_stderr()
    fix_ambpdb.run('%s_4tleap.pdb' % self.base,
                   '%s_new.pdb' % self.base,
                   '%s_new2.pdb' % self.base,
      )
    if save_cpp_traj_prmtop:
      shutil.copyfile("%s_asu.prmtop" % self.base,
                      "4cpptraj_%s.prmtop" % self.base,
                      )
    return 0

  #add cryst1 and sscale
  def finalizePdb(self,
                  pdb_filename=None, #cryst,base):
                  ):
    assert self.pdb_hierarchy
    assert self.cryst1
    assert self.base
    if pdb_filename:
      pdb_inp       = pdb.input(pdb_filename)
      pdb_hierarchy = pdb_inp.construct_hierarchy()
    else:
      # this returns the altloc to the model so not good!!!
      pdb_hierarchy = self.pdb_hierarchy
      assert 0
    pdbstring=pdb_hierarchy.as_pdb_string(crystal_symmetry=self.cryst1)
    new_file = '4phenix_'+self.base+'.pdb'
    print 'Writing 4phenix file',new_file
    f=open(new_file,'w')
    f.write(pdbstring)
    f.close()
    #~ import code; code.interact(local=locals())
    return 0

  def uc(self, redq=False):
    #add SMTRY/CRYST1 to 4tleap.pdb -> 4UnitCell.pdb
    assert self.base
    uc_pdb_file="%s_4UnitCell.pdb" % self.base
    with open(uc_pdb_file, "wb") as fout:
      with open(self.pdb_filename) as fin:
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
      with open("%s_new.pdb" % self.base) as fin:
        for line in fin:
          fout.write(line)
    tleap_pdb_file = "%s_4tleap_uc.pdb" % self.base
    run_UnitCell(uc_pdb_file, tleap_pdb_file)
    self.run_tleap(tleap_pdb_file,
                   output_base='uc',
                   reorder_residues='off',
                   #logfile='tleap_uc.log',
                   redq=redq,
                   )
    self.run_ChBox('uc')
    os.rename('%s_uc.rst7'   % self.base, '4amber_%s.rst7'   % self.base)
    os.rename('%s_uc.prmtop' % self.base, '4amber_%s.prmtop' % self.base)

  def run_minimise(self, option=None):
    assert self.base
    if option is None: return

    inputs = {"amber_h" : """Initial minimization
      &cntrl
       ntwx   = 0, ntb    = 1, cut    = 9.0,     nsnb   = 10,
       ntr    = 1, restraint_wt = 50.0, restraintmask ='!@H=',
       imin   = 1, maxcyc =1000, ncyc   = 200, ntmin  = 1,
      &end
      """,
              "amber_all" : """Initial minimization
      &cntrl
       ntwx   = 0, ntb    = 1, cut    = 9.0,     nsnb   = 10,
       imin   = 1, maxcyc = 1000, ncyc   = 200, ntmin  = 1,
      &end
      """
    }

    if option in ["amber_h", "amber_all"]:
    
      input_file = '%s_%s.in' % (self.base, option)
      f=open('%s_%s.in' % (self.base, option), 'wb')
      f.write(inputs[option])
      f.close()
      cmd='sander -O -i %s -p %s_asu.prmtop -c %s_asu.rst7 -o %s_%s.out \
           -ref %s_asu.rst7 -r %s_asu_%s.rst7' % (
             input_file,
             self.base,
             self.base,
             self.base,
             option,
             self.base,
             self.base,
             option,
             )
      print_cmd(cmd)
      ero=easy_run.fully_buffered(cmd)
      assert (ero.return_code == 0)
      self.run_ChBox('asu_%s' % option)
      cmd='ambpdb -p %s_asu.prmtop < %s_asu_%s.rst7 > %s_new.pdb' % (
        self.base,
        self.base,
        option,
        self.base,
        )
      print_cmd(cmd)
      ero=easy_run.fully_buffered(cmd)
      assert (ero.return_code == 0)
      ero.show_stdout()
      ero.show_stderr()
      fix_ambpdb.run('%s_4tleap.pdb' % self.base,
                     '%s_new.pdb' % self.base,
                     '%s_new2.pdb' % self.base,
        )
      self.finalizePdb(pdb_filename='%s_new2.pdb' % self.base)

    elif option=="phenix_all":
      cmd='phenix.geometry_minimization 4phenix_%s.pdb amber.use_amber=True \
           amber.topology_file_name=4amber_%s.prmtop \
           amber.coordinate_file_name=4amber_%s.rst7  \
           output_file_name_prefix=4phenix_%s_minPhenix ' % tuple([self.base]*4)
      restraints = "%s.ligands.cif" % self.base
      if os.path.exists(restraints):
        cmd += " %s" % restraints
      print_cmd(cmd)
      ero=easy_run.fully_buffered(cmd).raise_if_errors() #.return_code
      assert (ero.return_code == 0)
      ero.show_stdout()
      ero.show_stderr()
      os.rename('4phenix_%s_minPhenix.pdb' % self.base,
                '4phenix_%s.pdb' % self.base)
    return 0

  def check_special_positions(self):
    pdb_file = '4phenix_%s.pdb' % self.base
    print 'checking special positions in %s' % pdb_file
    xrs = iotbx.pdb.input(file_name=pdb_file).xray_structure_simple()
    site_symmetry_table = xrs.site_symmetry_table()
    if site_symmetry_table.n_special_positions() > 0:
      print  "WARNING: The following atoms occupy special positions."
      for i in site_symmetry_table.special_position_indices():
        print  "  Atom %d" %i
      print  "It may be a good idea to inspect manually and remove"
      print  "atoms from special positions or rerun minimization."
    # the following high gradient checker does not work yet
    # rst = Rst7.open('4amber_%s.rst7' %base)     #box
    # sander.setup('4amber_%s.prmtop' %base, rst.coordinates, rst.box, sander.pme_input())
    # crystal_symmetry = xrs.crystal_symmetry()
    # sites_cart = xrs.sites_cart()
    # sites_cart_uc=expand_coord_to_unit_cell(sites_cart, crystal_symmetry)
    # sander.set_positions( list(sites_cart_uc.as_double()) )
    # ene, frc = sander.energy_forces()
    # sander.cleanup()
    # big_frcs = [(i, frc_i) for i, frc_i in enumerate(frc) if frc_i>100000]
    # if len(big_frcs) >0:
    #   print >> sys.stderr, "WARNING: The following atoms have very high gradients."
    #   for i in big_frcs:
    #     print >> sys.stderr, "%d %f" %(i[0]/3, i[1])

  def run_clean(self):
    files_to_clean = """
      4tleap_nonprot.pdb
      4tleap.pdb
      4tleap_renum.txt
      4tleap_sslink
      4tleap_uc.pdb
      4UnitCell.pdb
      init.pdb
      leap.log
      new2.pdb
      new.pdb
      test.eff
      tleap_asu.log
      tleap.in
      tleap_uc.log
      asu.prmtop
      asu.rst7
      uc.prmtop
      uc.rst7
      sqm.pdb
      sqm.out
      sqm.in
      4antechamber
      amber_all.out
      amber_all.in
      asu.prmtop
      asu.rst7
      """
    import glob

    print "\nChecking for files to remove"
    for filename in files_to_clean.strip().split():
      for pre in ["", "%s_" % self.base]:
        for filename in glob.glob("%s%s" % (pre, filename)):
          print '  removing' , filename
          os.remove(filename)

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
  
# run elbow and antechamber
def _run_elbow_antechamber(pdb_hierarchy,
                           residue_name,
                           use_am1_and_maxcyc_zero=False,
                           debug=False):
  pdb_mol = get_molecule_from_hierarchy(pdb_hierarchy, residue_name)
  names = []
  for atom in pdb_mol:
    names.append(atom.name.strip())
  pdb_set = set(names)
  print >> sys.stderr, "\n=================================================="
  print >> sys.stderr, "Running elbow/antechamber for %s " %residue_name
  print >> sys.stderr, "=================================================="
  if debug and 0:
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
    if 0:
      filename = "4elbow.pdb"
      pdb_mol.WritePDB(filename)
      mol = builder.run(filename=filename,
                        id=residue_name,
                        no_output=True,
                        silent=True,
        )
    else:
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

  cmds = []
  cmd='antechamber -i 4antechamber_%s.pdb -fi pdb -o %s.mol2 -fo mol2 \
      -nc %d -m %d -s 2 -pf y -c bcc -at gaff' \
      %(residue_name, residue_name, mol.charge, mol.multiplicity)
  if use_am1_and_maxcyc_zero:
    cmd += " -ek 'qm_theory='AM1', grms_tol=0.0005, scfconv=1.d-10, maxcyc=0, ndiis_attempts=700,'"
  cmds.append(cmd)
  if not use_am1_and_maxcyc_zero:
    cmd='antechamber -i sqm.pdb -fi pdb -o %s.mol2 -fo mol2 \
      -nc %s -m %d -s 2 -pf y -c bcc -at gaff' \
      %(residue_name, mol.charge, mol.multiplicity)    
    cmds.append(cmd)
  
  for cmd in cmds:
    print_cmd(cmd)
    ero=easy_run.fully_buffered(cmd)
    stdo = StringIO.StringIO()
    ero.show_stdout(out=stdo)
    for line in stdo.getvalue().splitlines():
      if line.find('APS')>-1:
        print line
      if line.find('Error')>-1:
        raise Sorry(line)
      
  cmd='parmchk2 -i %s.mol2 -f mol2 -o %s.frcmod' %(residue_name, residue_name)
  print_cmd(cmd)
  easy_run.fully_buffered(cmd)

def run_UnitCell(input_file,output_file):
  cmd='UnitCell -p %s -o %s' %(input_file, output_file)
  print_cmd(cmd)
  ero=easy_run.fully_buffered(cmd)
  ero.show_stdout()
  ero.show_stderr()

def run(rargs):
  working_params = setup_options_args(rargs)
  inputs = working_params.amber_prep.input
  actions = working_params.amber_prep.actions
  base = get_output_preamble(working_params)
  amber_prep_runner = amber_prep_run_class(base)
  amber_prep_runner.initializePdb(inputs.pdb_file_name)
  invalid = amber_prep_runner.validatePdb()
  if invalid: raise Sorry( 'PDB input is not "valid"' )
  amber_prep_runner.run_pdb4amber()
  #
  amber_prep_runner.run_elbow_antechamber(nproc=inputs.nproc)
  #
  print >> sys.stderr, "\n=================================================="
  print >> sys.stderr, "Preparing asu files and 4phenix_%s.pdb" % base
  print >> sys.stderr, "=================================================="
  amber_prep_runner.run_tleap('%s_4tleap.pdb' % amber_prep_runner.base,
                              'asu',
                              #ns_names,
                              reorder_residues='on',
                              #logfile='tleap_asu.log',
                              redq=actions.redq,
    )
  amber_prep_runner.run_ChBox("asu")
  amber_prep_runner.run_ambpdb(save_cpp_traj_prmtop=actions.save_cpp_traj_prmtop)
  amber_prep_runner.finalizePdb(
    pdb_filename='%s_new2.pdb' % amber_prep_runner.base)
  
  print >> sys.stderr, "\n=================================================="
  print >> sys.stderr, "Preparing uc files: %s.prmtop and %s.rst7" %(base,base)
  print >> sys.stderr, "=================================================="
  amber_prep_runner.uc(redq=actions.redq)
  if actions.minimise == "off":
    pass
  else:
    print >> sys.stderr, "\n=================================================="
    print >> sys.stderr, "Minimizing input coordinates."
    print >> sys.stderr, "=================================================="
    amber_prep_runner.run_minimise(actions.minimise)
  amber_prep_runner.check_special_positions()
  if actions.clean == "on": amber_prep_runner.run_clean()
  outl = "\n\nExample\n\n  phenix.geometry_minimization"
  outl += " 4phenix_%s.pdb use_amber=True" % (
    amber_prep_runner.base,
    )
  outl += " amber.topology_file_name=4amber_%s.prmtop" % amber_prep_runner.base
  outl += " amber.coordinate_file_name=4amber_%s.rst7" % amber_prep_runner.base
  outl += "\n\n\n"
  print outl
  return '4phenix_%s.pdb' % amber_prep_runner.base

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
