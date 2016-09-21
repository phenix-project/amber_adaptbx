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
from amber_adaptbx import amber_library_server
from amber_adaptbx.scripts.les_build import LESBuilder
from libtbx.utils import Sorry
import libtbx.load_env

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
      antechamber
        .caption = Options for Amber program antechamber
      {
        prefer_input_method = chemical_component elbow *Auto
          .type = choice
          .caption = Control to first input chosen for antechamber
      }
    }
    actions
    {
      minimise = amber_all amber_h phenix_all *off
        .type = choice
      amber_min_options = ''
        .type = str
        .caption = add more amber options for minimization for amber minimization. If not specified, use default values
        .help = add more amber options for minimization for amber minimization. If not specified, use default values
      clean = True
        .type = bool
      redq = False
        .type = bool
      LES = False
        .type = bool
      use_reduce = True
        .type = bool
        .caption = Run reduce on the input pdb file to place hydrogens
        .help = Run reduce on the input pdb file to place hydrogens
    }
    output
    {
      file_name = None
        .type = path
    }
  }
  """
#      save_cpp_traj_prmtop = True
#        .type = bool

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
  #working_phil.show()
  working_params = working_phil.extract()
  working_params.amber_prep.input.pdb_file_name = pdbs[0]
  #check_working_params(working_params)
  preamble = get_output_preamble(working_params)
  #print "  Writing effective parameters to %s.eff\n" % preamble
  #working_phil.format(python_object=working_params).show()
  #print "#phil __ON__"
  #master_phil.fetch_diff(source=master_phil.format(
  #    python_object=working_params)).show()
  #print "#phil __OFF__\n\n"
  f=file("%s.eff" % preamble, "wb")
  f.write(working_phil.format(python_object=working_params).as_str())
  f.close()
  return working_params

def test_files_exist(filenames):
  #return # doesn't work for minimi=amber_all
  for filename in filenames:
    if not os.path.exists(filename):
      raise Sorry("Filename %s not found" % filename)

def check_required_output_filenames(filenames=None,
                                    error_display_file=None,
                                   ):
  # only tleap currently
  if filenames is None: return
  print 'Checking output filenames'
  for filename in filenames:
    print '  file : %s' % filename
    if not os.path.exists(filename):
      s = "  Output file not present : %s" % filename
      if error_display_file:
        s += '\n  Check contents of "%s"' % error_display_file
      raise Sorry(s)
    if os.stat(filename).st_size==0:
      s = "  Output file is empty : %s" % filename
      if error_display_file:
        s += '\n  Check contents of "%s"' % error_display_file
      raise Sorry("  Output file is empty : %s" % filename)

def print_cmd(cmd, verbose=False):
  print "\n~> %s\n" % cmd
  if verbose:
    print '\nAMBERHOME: %s' % os.environ.get("AMBERHOME", None)
    path = os.environ.get("PATH", None)
    print '\nPATH: %s' % ('\n    : '.join(path.split(":"))) 

def get_chemical_components_file_name(code):
  cc_dir = libtbx.env.dist_path("chem_data", default=None)
  if not cc_dir: return None
  cc_file_name = os.path.join(cc_dir,
                                 'chemical_components',
                                 code[0].lower(),
                                 'data_%s.cif' % code.upper(),
                                 )
  if os.path.exists(cc_file_name):
    return cc_file_name
  return None

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

  def curate_model(self, remove_alt_confs=False):
    assert self.pdb_hierarchy
    from mmtbx import pdbtools
    if remove_alt_confs:
      # performs removal in place
      pdbtools.remove_alt_confs(self.pdb_hierarchy,
                                always_keep_one_conformer=True,
      )

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

  def run_elbow_antechamber(self,
                            ns_names=[],
                            nproc=1,
                            prefer_input_method=None,
                            debug=False):
    assert self.pdb_hierarchy
    if nproc>1:
      print "\n\tParallel processes not implemented\n"
    for residue_name in ns_names:
      if prefer_input_method:
        if prefer_input_method=="chemical_component":
          _run_antechamber_ccif( residue_name )
          continue
        elif prefer_input_method=="elbow":
          _run_elbow_antechamber(self.pdb_hierarchy, residue_name, debug=debug)
          continue
      if amber_library_server.is_in_components_lib(residue_name):
        print """
  Residue "%s" already in amber monomer library. Skipping elbow/antechamber 
    run for this residue.
        """ % residue_name
        continue
      elif os.path.isfile('%s.mol2' % residue_name):
        print "%s.mol2 is present. Skipping elbow/antechamber run for this residue.\n" %residue_name
        continue
      elif get_chemical_components_file_name(residue_name):
        _run_antechamber_ccif( residue_name )
      else:
        _run_elbow_antechamber(self.pdb_hierarchy, residue_name, debug=debug)
    return 0

  #=================================
  # prepare tleap input:
  #=================================
  def run_tleap(self,
                input_pdb,
                output_base,
                ns_names,
                gaplist,
                sslist,
                reorder_residues,
                logfile="leap.log",
                redq=False,
                verbose=False,
                ):
    def _parse_tleap_logfile(logfile):
      errors = []
      warnings = []
      fatals=[]
      f=file(logfile, "rb")
      lines=f.read()
      f.close()
      warnings = {
        #"Failed to generate parameters": None,
        "WARNING: The unperturbed charge of the unit:" : \
        "Note: Highly charged molecules with no indication of how they are to be neutralized might cause problems",
        }
      for line in lines.splitlines():
        if line.find("FATAL:")>-1:
          fatals.append(line)
        for warning, msg in warnings.items():
          if line.find(warning)>-1:
            errors.append(line)
            if msg: errors.append(" > %s" % msg)
      if errors:
        print "Errors and warnings in tleap"
        # print "  check logfile %s for explaination\n" % logfile
        for line in errors:
          print "  %s\n" % line
      if fatals:
        print "Fatal errors in tleap"
        for line in fatals:
          print line
        raise Sorry("Fatal errors in tleap")
      return fatals
    def _check_tleap_output(lines):
      stop = ["Could not open file"]
      for line in lines:
        for s in stop:
          if line.find(s)>-1:
            print_cmd('" --- Testing environment ---"', verbose=True)
            raise Sorry('tleap error : "%s"' % line)
    tleap_input_file = "%s_%s_tleap_input_run" % (self.base, output_base)
    f=file(tleap_input_file, "wb")
    f.write('logFile %s\n' % logfile)

    # Following should be true in AmberTools14/15:
    amber_dir = libtbx.env.dist_path("amber")
    if(os.path.isfile(os.path.join(amber_dir,
                                   'dat',
                                   'leap',
                                   'cmd',
                                   'leaprc.ff14SB',
                                   ))):
      raise Sorry('Amber environment appears to be older than AmberTools16; quitting')

    # Now we can assume that we are dealing with AmberTools16:
    f.write('source leaprc.protein.ff14SB\n')
    f.write('source leaprc.DNA.OL15\n')
    f.write('source leaprc.RNA.OL3\n')
    f.write('source leaprc.water.tip3p\n')
    f.write('source leaprc.gaff2\n')
    #  (for the future: have some mechanism for modifying the above list)

    for res in ns_names:
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
    #
    #  process gaplist
    #
    if gaplist:
      for d,res1,resid1,res2,resid2 in gaplist:
        f.write('deleteBond x.%d.C x.%d.N\n' % ( resid1, resid2 ) )
    #
    #  process sslist
    #
    if sslist:
      for resid1,resid2 in sslist:
        f.write('bond x.%d.SG x.%d.SG\n' % ( resid1, resid2 ) )

    f.write('saveAmberParm x %s.prmtop %s.rst7\n' %(
        "%s_%s" % (self.base, output_base),
        "%s_%s" % (self.base, output_base),
        )
      )
    f.write('quit\n')
    f.close()
    #
    # strangely tleap appends to the logfile so must delete first
    #
    if os.path.exists(logfile): os.remove(logfile)
    cmd = 'tleap -f %s' % tleap_input_file
    print_cmd(cmd)
    ero=easy_run.fully_buffered(cmd)
    _check_tleap_output(ero.stdout_lines)
    if verbose:
      ero.show_stdout()
      ero.show_stderr()
    # check output
    fatals = _parse_tleap_logfile(logfile)
    check_required_output_filenames(["%s_%s.prmtop" % (self.base, output_base),
                                     "%s_%s.rst7" % (self.base, output_base),
                                     ],
                                    error_display_file="leap.log",
                                    )
    # remove files
    #if os.path.exists(tleap_input_file): os.remove(tleap_input_file)
    return 0

  #change box size in rst7 file
  def run_ChBox(self, output_base):
    assert self.cryst1
    uc = self.cryst1.unit_cell().parameters()
    # note: dangerous to use same file for input and output here?
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

  def _pdb_hierarchy_and_rename_wat(self, filename):
    pdb_inp=iotbx.pdb.input(file_name=filename)
    pdb_hierarchy = pdb_inp.construct_hierarchy(sort_atoms=False)
    # the -bres option in ambpdb does not (yet) change "WAT" to "HOH"
    for atom_group in pdb_hierarchy.atom_groups():
      if atom_group.resname in ["WAT"]:
        atom_group.resname = "HOH"
    return pdb_inp, pdb_hierarchy

  def _match_hierarchies_and_transfer_to(self,
                                         hierachy1,
                                         hierachy2,
                                         transfer_b=False,
                                         transfer_occ=False,
                                         transfer_chain_id=False,
                                         transfer_xyz=False,
                                         ):
    #match residues based on resseq and resname
    #match atoms based on name an i_seq
    for chain_post in hierachy2.chains():
      for resi_post in chain_post.conformers()[0].residues():
        for atom_post in resi_post.atoms():
          for chain_pre in hierachy1.chains():
            for resi_pre in chain_pre.conformers()[0].residues():
              if ( resi_pre.resseq==resi_post.resseq and
                   resi_pre.resname.strip()==resi_post.resname.strip()
                   ):
                for atom_pre in resi_pre.atoms():
                  if atom_pre.name == atom_post.name:
                    if transfer_b: atom_post.b=atom_pre.b
                    if transfer_occ: atom_post.occ=atom_pre.occ
                    if transfer_chain_id: chain_post.id=chain_pre.id
                    if transfer_xyz:
                      atom_post.xyz=(atom_pre.xyz[0],
                                     atom_pre.xyz[1],
                                     atom_pre.xyz[2])

  # make pdb
  def run_ambpdb(self):
    assert self.base
    cmd='ambpdb -bres -p %s_asu.prmtop < %s_asu.rst7 > %s_new.pdb' % tuple(
      [self.base]*3
      )
    print_cmd(cmd)
    ero=easy_run.fully_buffered(cmd)
    ero.show_stdout()
    ero.show_stderr()

    pdb_pre, pdb_h_pre = self._pdb_hierarchy_and_rename_wat('%s_4tleap.pdb' % self.base)
    pdb_post, pdb_h_post = self._pdb_hierarchy_and_rename_wat('%s_new.pdb' % self.base)

    self._match_hierarchies_and_transfer_to(pdb_h_pre,  # from
                                            pdb_h_post, # to
                                            transfer_b=True,
                                            transfer_occ=True,
                                            transfer_chain_id=True,
                                            )
                                          
    pdb_h_post.write_pdb_file(file_name='%s_new2.pdb' %self.base,
                              append_end=True,
                              crystal_symmetry=pdb_pre.crystal_symmetry(),
                              )
    return 0

  #add cryst1 and sscale; creates 4phenix_base_type_.pdb file
  def finalizePdb(self,
                  pdb_filename=None, sort_atoms=True, type='',
                  ):
    assert self.pdb_hierarchy
    assert self.cryst1
    assert self.base
    if pdb_filename:
      pdb_inp       = pdb.input(pdb_filename)
      pdb_hierarchy = pdb_inp.construct_hierarchy(sort_atoms=sort_atoms)
    else:
      # this returns the altloc to the model so not good!!!
      pdb_hierarchy = self.pdb_hierarchy
      assert 0
    pdbstring=pdb_hierarchy.as_pdb_string(crystal_symmetry=self.cryst1)
    new_file = '4phenix_'+self.base+type+'.pdb'
    print 'Writing 4phenix file',new_file
    f=open(new_file,'w')
    f.write(pdbstring)
    f.close()
    #~ import code; code.interact(local=locals())
    return 0

  def write_remark_290(self):
    assert self.pdb_filename
    with open(self.pdb_filename) as fin:
      lines = fin.readlines()
      smtry = [line for line in lines if "SMTRY" in line]
      smtry = ''.join(smtry)
      rem = 'remark_290.txt'
      if not smtry:
        print '"%s"' % smtry
        raise Sorry("REMARK 290 SMTRY1,2,3 records required")
      else:
        # print smtry
        file(rem, "wb").write(smtry)

  def uc(self, redq=False):
    #add SMTRY/CRYST1 to 4tleap.pdb -> 4UnitCell.pdb
    assert self.base
    uc_pdb_file="%s_4UnitCell.pdb" % self.base
    with open(uc_pdb_file, "wb") as fout:
      with open(self.pdb_filename) as fin:
        lines = fin.readlines()
        smtry = [line for line in lines if "SMTRY" in line]
        smtry = ''.join(smtry)
        rem = 'remark_290.txt'
        if not smtry:
          if os.path.exists(rem):
            smtry=file(rem, "rb").read()
          else:
            print '"%s"' % smtry
            raise Sorry("REMARK 290 SMTRY1,2,3 records required")
        else:
          # print smtry
          file(rem, "wb").write(smtry)
        fout.write(smtry)
        # import code; code.interact(local=locals())
        cryst1card = [line for line in lines if "CRYST1" in line]
        if len(cryst1card) <1:
          raise Sorry("CRYST1 record required")
        fout.write(cryst1card[0])
      with open("4phenix_%s.pdb" % self.base) as fin:
        for line in fin:
          if not "CRYST1" in line:
            fout.write(line)
    tleap_pdb_file = "%s_4tleap_uc.pdb" % self.base

    # temporary file: will be output from UnitCell, input to pdb4amber:
    tleap_pdb_file1 = "%s_4tleap_uc1.pdb" % self.base
    run_UnitCell(uc_pdb_file, tleap_pdb_file1)

    #run back through pdb4amber to get new CONECT records for SS bonds:
    ns_names,gaplist,sslist=pdb4amber.run(
          tleap_pdb_file, tleap_pdb_file1, arg_elbow=True,
          )
    self.run_tleap(tleap_pdb_file,
                   output_base='uc',
                   ns_names=ns_names,
                   gaplist=gaplist,
                   sslist=sslist,
                   reorder_residues='off',
                   #logfile='tleap_uc.log',
                   redq=redq,
                   )
    self.run_ChBox('uc')
    os.rename('%s_uc.rst7'   % self.base, '4amber_%s.rst7'   % self.base)
    os.rename('%s_uc.prmtop' % self.base, '4amber_%s.prmtop' % self.base)

  def run_minimise(self, option=None, is_LES=False, amber_min_options=''):
    assert self.base
    if option is None: return

    inputs = {"amber_h" : """Initial minimization
      &cntrl
       ntwx   = 0, ntb    = 1, cut    = 9.0,     nsnb   = 10,
       ntr    = 1, restraint_wt = 50.0, restraintmask ='!@H=',
       imin   = 1, maxcyc = 1000, ncyc   = 200, ntmin  = 1, ntxo = 1,
      /
      """,
              "amber_all" : """Initial minimization
      &cntrl
       ntwx   = 0, ntb    = 1, cut    = 9.0,     nsnb   = 10,
       imin   = 1, maxcyc = 500, ncyc   = 200, ntmin  = 1, ntxo = 1,
       ntpr=50, ntr=1, restraint_wt=2.0, restraintmask=':*',
      /
      """
    }
    if amber_min_options:
      inputs['amber_h'] = inputs['amber_h'].replace('/', amber_min_options + '\n /')
      inputs['amber_all'] = inputs['amber_all'].replace('/', amber_min_options + '\n /')
    prmtop_file_name = '4amber_%s.prmtop' % self.base
    rst7_file_name = '4amber_%s.rst7' % self.base
    output_rst7_file_name = ' %s_%s.rst7' % (self.base, option)
    if is_LES:
      prmtop_file_name = prmtop_file_name.replace('.prmtop', '.LES.prmtop')
      rst7_file_name = rst7_file_name.replace('.rst7', '.LES.rst7')
      output_rst7_file_name = output_rst7_file_name.replace('.rst7', '.LES.rst7')

    if option in ["amber_h", "amber_all"]:

      input_file = '%s_%s.in' % (self.base, option)
      f=open('%s_%s.in' % (self.base, option), 'wb')
      f.write(inputs[option])
      f.close()
      cmd='sander -O -i %s -p %s -c %s -o %s_%s.out \
           -ref %s -r %s' % (
             input_file,
             prmtop_file_name,
             rst7_file_name,
             self.base,
             option,
             rst7_file_name,
             output_rst7_file_name
             )
      if is_LES:
        cmd = cmd.replace('sander', 'sander.LES')
      print_cmd(cmd)
      # test function that may be useful...
      test_files_exist([input_file,
                        prmtop_file_name,
                        rst7_file_name,
                      ])
      ero=easy_run.fully_buffered(cmd)
      assert (ero.return_code == 0)
      test_files_exist([input_file,
                        prmtop_file_name,
                        rst7_file_name,
                      ])

      cmd='ambpdb -bres -p %s < %s > %s_new.pdb' % (
        prmtop_file_name,
        rst7_file_name,
        self.base,
        )
      print_cmd(cmd)
      ero=easy_run.fully_buffered(cmd)
      assert (ero.return_code == 0)
      ero.show_stdout()
      ero.show_stderr()

      pdb_pre, pdb_h_pre = self._pdb_hierarchy_and_rename_wat('4phenix_%s.pdb' % self.base)
      pdb_post, pdb_h_post = self._pdb_hierarchy_and_rename_wat('%s_new.pdb' % self.base)

      # there is a function that will transfer the coordinates from one PDB
      # hierarchy to another but the atoms have to be the same number & order
      self._match_hierarchies_and_transfer_to(pdb_h_post, # from
                                              pdb_h_pre,  # to
                                              transfer_xyz=True,
                                             )

      pdb_h_pre.write_pdb_file(file_name='%s_new2.pdb' % self.base,
                               append_end=True,
                               crystal_symmetry=pdb_pre.crystal_symmetry(),
                               )

      self.finalizePdb(pdb_filename='%s_new2.pdb' % self.base,
              sort_atoms=False, type='.min')

    elif option=="phenix_all":
      if is_LES:
        cmd='phenix.geometry_minimization 4phenix_%s.LES.pdb amber.use_amber=True \
             amber.topology_file_name=%s \
             amber.coordinate_file_name=%s  \
             output_file_name_prefix=4phenix_%s_minPhenix' % (self.base, prmtop_file_name, rst7_file_name, self.base)
      else:
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

  def check_special_positions(self,is_LES):
    if is_LES:
      pdb_file = '4phenix_%s.LES.pdb' % self.base
    else: 
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

  def run_clean(self):
    files_to_clean = """
      4tleap_nonprot.pdb
      4tleap_uc_nonprot.pdb
      4tleap.pdb
      4tleap_uc1.pdb
      4tleap_renum.txt
      4tleap_uc_renum.txt
      4tleap_uc_sslink
      4tleap_sslink
      4tleap_uc.pdb
      uc_tleap_input_run
      asu_tleap_input_run
      4UnitCell.pdb
      leap.log
      new2.pdb
      new.pdb
      tleap_asu.log
      tleap.in
      tleap_uc.log
      sqm.pdb
      sqm.out
      sqm.in
      4antechamber
      amber_all.in
      mdinfo
      asu.prmtop
      asu.rst7
      remark_290.txt
      addles.in
      reduce_info.log
      """
    import glob

    # print "\nChecking for files to remove"
    for filename in files_to_clean.strip().split():
      for pre in ["", "%s_" % self.base]:
        for filename in glob.glob("%s%s" % (pre, filename)):
          # print '  removing' , filename
          os.remove(filename)
    for s in ['%s.eff',
              '%s_curated.pdb',
              '%s_uc.pdb',
              '%s_uc_H.pdb',
              '%sab.rst7',
              '4amber_%s.pdb',
              ]:
      if os.path.isfile( s % self.base ):
        os.remove( s % self.base)

  def write_pdb_hierarchy(self,pdb_file_name):
    assert self.pdb_hierarchy
    self.pdb_hierarchy.write_pdb_file(
      file_name=pdb_file_name,
      append_end=True,
      crystal_symmetry=self.pdb_inp.crystal_symmetry(),
    )

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

# run antechamber from a components.cif file:
def _run_antechamber_ccif( residue_name,
                           use_am1_and_maxcyc_zero=True,
                           debug=False):
  print "\n=================================================="
  print "Running antechamber_ccif for %s " %residue_name
  print "=================================================="

  ccif = get_chemical_components_file_name(residue_name)
  cmds = []
  cmd='antechamber -i %s -fi ccif -bk %s -o %s.mol2 -fo mol2 \
      -s 2 -pf y -c bcc -at gaff2'  %(ccif, residue_name, residue_name)
  if use_am1_and_maxcyc_zero:
    cmd += ' -ek "qm_theory=\'AM1\', grms_tol=0.0005, scfconv=1.d-10, maxcyc=0, ndiis_attempts=700,"'
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

  cmd='parmchk2 -s 2 -i %s.mol2 -f mol2 -o %s.frcmod' %(residue_name,
                                                        residue_name)
  print_cmd(cmd)
  easy_run.fully_buffered(cmd)

  # should there be a check for output???

# run elbow and antechamber
def _run_elbow_antechamber(pdb_hierarchy,
                           residue_name,
                           use_am1_and_maxcyc_zero=True,
                           debug=False):
  pdb_mol = get_molecule_from_hierarchy(pdb_hierarchy, residue_name)
  names = []
  for atom in pdb_mol:
    names.append(atom.name.strip())
  pdb_set = set(names)
  print "\n=================================================="
  print "Running elbow/antechamber for %s " %residue_name
  print "=================================================="
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
                        #no_output=True,
                        #pH=8, # special Amber number...
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

  mol.WritePDB('4antechamber_%s.pdb' %residue_name,
               pymol_pdb_bond_order=False,
              )
  mol.Multiplicitise()
  print mol.DisplayBrief()

  cmds = []
  cmd='antechamber -i 4antechamber_%s.pdb -fi pdb -o %s.mol2 -fo mol2 \
      -nc %d -m %d -s 2 -pf y -c bcc -at gaff2' \
      %(residue_name, residue_name, mol.charge, mol.multiplicity)
  if use_am1_and_maxcyc_zero:
    cmd += ' -ek "qm_theory=\'AM1\', grms_tol=0.0005, scfconv=1.d-10, maxcyc=0, ndiis_attempts=700,"'
  cmds.append(cmd)
  if not use_am1_and_maxcyc_zero:
    cmd='antechamber -i sqm.pdb -fi pdb -o %s.mol2 -fo mol2 \
      -nc %s -m %d -s 2 -pf y -c bcc -at gaff2' \
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

  cmd='parmchk2 -s 2 -i %s.mol2 -f mol2 -o %s.frcmod' %(residue_name, residue_name)
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
  amber_prep_runner.write_remark_290()
  amber_prep_runner.curate_model(remove_alt_confs=(not actions.LES))
  # need to write PDB for some of the other methods
  current_pdb_file_name = inputs.pdb_file_name.replace(
    '.pdb',
    '_curated.pdb',
    )
  amber_prep_runner.write_pdb_hierarchy(current_pdb_file_name)

  print "\n=================================================="
  print "Running pdb4amber on %s" % inputs.pdb_file_name
  print "=================================================="

  tleap_input_pdb = "%s_4tleap.pdb" % base
  log = []
  ns_names,gaplist,sslist=pdb4amber.run(tleap_input_pdb,
                                        current_pdb_file_name,
                                        arg_elbow=True,
                                        arg_reduce=actions.use_reduce,
                                        #log=log,
      )

  print "\n=================================================="
  print "Setting up library files for non-standard residues"
  print "=================================================="

  amber_prep_runner.run_elbow_antechamber(
    ns_names, 
    nproc=inputs.nproc, 
    prefer_input_method=inputs.antechamber.prefer_input_method,
  )
  #
  print "\n=================================================="
  print "Preparing asu files and 4phenix_%s.pdb" % base
  print "=================================================="
  
  # at this point, we can ignore any gaps found by pdb4amber, since
  #   we are just creating pdb files for phenix (which automatically
  #   finds gaps itself), and for making the unit cell pdb file; in
  #   the latter case, the second run of pdb4amber (with the unit cell
  #   pdb file) will find the gaps needed for the construction of the
  #   4amber_base.prmtop file.
  #
  # N.B.: this means that the base_asu.prmtop file should never be used!
  #   we might want to make user that this file is always removed.

  dummy_gaplist=[]
  amber_prep_runner.run_tleap('%s_4tleap.pdb' % amber_prep_runner.base,
                              'asu',
                              ns_names,
                              dummy_gaplist,
                              sslist,
                              reorder_residues='on',
                              #logfile='tleap_asu.log',
                              redq=actions.redq,
    )
  amber_prep_runner.run_ChBox("asu")
  amber_prep_runner.run_ambpdb()   # (note: only called once)
  amber_prep_runner.finalizePdb(
      pdb_filename='%s_new2.pdb' % amber_prep_runner.base, sort_atoms=False)

  print "\n============================================================"
  print "Preparing unit cell files: 4amber_%s.prmtop and 4amber_%s.rst7" %(base,base)
  print "============================================================"

  amber_prep_runner.uc(redq=actions.redq)
  if actions.LES:
   pdb_file_name = working_params.amber_prep.input.pdb_file_name
   non_les_prmtop_file_name = '4amber_%s.prmtop' % amber_prep_runner.base
   non_les_rst7_file_name = '4amber_%s.rst7' % amber_prep_runner.base
   les_builder = LESBuilder(
           pdb_file_name,
           prmtop=non_les_prmtop_file_name,
           rst7_file=non_les_rst7_file_name)
   les_builder.run()

  if actions.minimise == "off":
    pass
  else:
    print "\n=================================================="
    print "Minimizing input coordinates."
    print "=================================================="
    amber_prep_runner.run_minimise(actions.minimise, actions.LES, 
        amber_min_options=actions.amber_min_options)

  amber_prep_runner.check_special_positions(actions.LES)

  if actions.clean: amber_prep_runner.run_clean()

  outl = "\n==================================================\n"
  outl += "Done.  Three new files have been made:\n"
  outl += "      4phenix_%s.pdb\n" % amber_prep_runner.base
  outl += "      4amber_%s.prmtop\n" % amber_prep_runner.base
  outl += "      4amber_%s.rst7\n" % amber_prep_runner.base
  outl += "==================================================\n\n"
  outl += "Example\n\n  phenix.refine"
  outl += " 4phenix_%s.pdb use_amber=True \\\n" % (
    amber_prep_runner.base,
    )
  outl += "    amber.topology_file_name=4amber_%s.prmtop \\\n" % amber_prep_runner.base
  outl += "    amber.coordinate_file_name=4amber_%s.rst7 \\\n" % amber_prep_runner.base
  outl += "    ....(other refinement keywords here)....."
  outl += "\n\n\n"
  if actions.LES:
    outl = (outl.replace('.pdb', '.LES.pdb')
                .replace('.prmtop', '.LES.prmtop')
                .replace('.rst7', '.LES.rst7'))
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
    amber_prep_runner.run_minimise(actions.minimise, actions.LES)
