# LIBTBX_SET_DISPATCHER_NAME phenix.AmberPrep

import os
import sys
import iotbx.pdb
import StringIO
from libtbx import phil
from libtbx.utils import Sorry
import libtbx.load_env
import libtbx.phil.command_line
from libtbx import easy_run
from elbow.command_line import builder

# try: import pdb4amber
# except ImportError:
#   raise Sorry('  Import error - Please check that AMBERHOME is set correctly: %s' % (
#     os.environ.get('AMBERHOME', None)))

from amber_adaptbx import pdb4amber

from amber_adaptbx import amber_library_server
from amber_adaptbx.utils import build_unitcell, \
  get_indices_convert_dict_from_array
from amber_adaptbx.les_builder.build import LESBuilder
import parmed

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
        prefer_input_method = chemical_component *elbow Auto
          .type = choice
          .caption = Control to first input chosen for antechamber
      }
    }
    actions
    {
      minimise = amber_all amber_h *off
        .type = choice
      minimization_options = ''
        .type = str
        .caption = add more amber options for minimization for amber minimization. If not specified, use default values
        .help = add more amber options for minimization for amber minimization. If not specified, use default values
      clean = True
        .type = bool
      redq = False
        .type = bool
        .caption = Use reduced-charge Amber force field
        .help = Use reduced-charge Amber force field
      LES = False
        .type = bool
      use_reduce = True
        .type = bool
        .caption = Run reduce on the input pdb file to place hydrogens
        .help = Run reduce on the input pdb file to place hydrogens
      use_glycam = False
        .type = bool
        .caption = Load GLYCAM carbohydrate force field
        .help = Load GLYCAM carbohydrate force field
      addles_input = ''
        .type = str
        .caption = User specify addles input filename. Optional.
        .help = User specify addles input filename. Optional.
      use_amber_unitcell = False
        .type = bool
        .help = 'If True, use UnitCell from Amber; else use expand_to_p1()'
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

master_params = master_phil_string  # need for auto documentation
master_phil = phil.parse(master_phil_string,
                         process_includes=True,
                         )


def setup_parser():
  from libtbx.option_parser import OptionParser
  usage = """
  phenix.AmberPrep 3a37.pdb minimise=amber_all
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
  preamble = preamble.replace(".pdb", "")
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
  if len(args) == 0:
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
    if os.path.isfile(arg):
      if iotbx.pdb.is_pdb_file(arg) or arg.endswith('.cif'):
        pdbs.append(arg)
      else:
        try:
          file_phil = phil.parse(file_name=arg)
        except RuntimeError:
          pass
        else:
          phils.append(file_phil)
    else:
      phil_args.append(arg)
      phils.append(argument_interpreter.process(arg))
  working_phil = master_phil.fetch(sources=phils)
  assert len(pdbs) == 1
  # working_phil.show()
  working_params = working_phil.extract()
  working_params.amber_prep.input.pdb_file_name = pdbs[0]
  # check_working_params(working_params)
  preamble = get_output_preamble(working_params)
  # print "  Writing effective parameters to %s.eff\n" % preamble
  # working_phil.format(python_object=working_params).show()
  # print "#phil __ON__"
  # master_phil.fetch_diff(source=master_phil.format(
  #    python_object=working_params)).show()
  # print "#phil __OFF__\n\n"
  with open("%s.eff" % preamble, "wb") as f:
    f.write(working_phil.format(python_object=working_params).as_str())
  return working_params


def test_files_exist(filenames):
  # return # doesn't work for minimi=amber_all
  for filename in filenames:
    if not os.path.exists(filename):
      raise Sorry("Filename %s not found" % filename)


def check_required_output_filenames(filenames=None,
                                    error_display_file=None,
                                    ):
  # only tleap currently
  if filenames is None:
    return
  print 'Checking output filenames'
  for filename in filenames:
    print '  file : %s' % filename
    if not os.path.exists(filename):
      s = "  Output file not present : %s" % filename
      if error_display_file:
        s += '\n  Check contents of "%s"' % error_display_file
      raise Sorry(s)
    if os.stat(filename).st_size == 0:
      s = "  Output file is empty : %s" % filename
      if error_display_file:
        s += '\n  Check contents of "%s"' % error_display_file
      raise Sorry(s)


def print_cmd(cmd, verbose=False):
  print "\n| ~> %s\n" % cmd
  if verbose:
    print '\nAMBERHOME: %s' % os.environ.get("AMBERHOME", None)
    path = os.environ.get("PATH", None)
    print '\nPATH: %s' % ('\n    : '.join(path.split(":")))


def get_chemical_components_file_name(code):
  cc_dir = libtbx.env.dist_path("chem_data", default=None)
  if not cc_dir:
    return None
  cc_file_name = os.path.join(cc_dir,
                              'chemical_components',
                              code[0].lower(),
                              'data_%s.cif' % code.upper(),
                              )
  if os.path.exists(cc_file_name):
    return cc_file_name
  return None


class AmberPrepRunner:

  def __init__(self, base_name, LES=False):
    self.base = base_name
    self.LES = LES
    self.final_prmtop_file = '4amber_%s.prmtop' % self.base
    self.final_rst7_file = '4amber_%s.rst7' % self.base
    self.final_order_file = '4amber_%s.order' % self.base
    self.final_pdb_file_4phenix = '4phenix_%s.pdb' % self.base
    self.non_les_prmtop_file_name = '4amber_%s.prmtop' % self.base
    self.non_les_rst7_file_name = '4amber_%s.rst7' % self.base

  def __repr__(self):
    outl = "AmberPrepRunner"
    outl += "\n  Base : %s" % self.base
    return outl

  # get box info & space group info
  def initialize_pdb(self, pdb_filename):
    self.pdb_filename = pdb_filename
    self.pdb_inp = iotbx.pdb.input(pdb_filename)
    self.cryst1 = self.pdb_inp.crystal_symmetry_from_cryst1()
    self.pdb_hierarchy = self.pdb_inp.construct_hierarchy()

  def asu_parm7_to_4phenix_pdb(self, pdb_filename):
    '''
    Combine information in xxxx_asu.{prmtop,rst7} with that in the input
       pdb file to create 4phenix_xxxx.pdb; (uses parmed)

    We want each atom in the 4phenix_xxxx.pdb file to have an occupancy
    that is the sum of the occupancies of the all the alternate conformers
    in the original pdb file.  This should make the total number of
    electrons in the 4phenix_xxxx.pdb file the same as in the original pdb
    file.

    Note:  the atom order will be the "Amber" atom order; use phenix
    routines after this to conver to phenix atom order.
    '''

    #  Note: here self.pdb_filename is typically the input pdb file
    template_parm = parmed.load_file(self.pdb_filename)
    parm = parmed.load_file("%s_asu.prmtop" % self.base, 
                         xyz="%s_asu.rst7" % self.base )

    # following kludge is needed to get got atom numbers; not needed if
    #   conversion to phenix order is done right after this
    # parm.write_pdb( pdb_filename )
    # parm = parmed.load_file( pdb_filename )

    # update occupancy to take into account alternate locations:
    for atom in template_parm.atoms:
      if atom.other_locations:
        # if atom has alternative locations, ParmEd saves them to 
        # other_locations (dict)
        atom.occupancy += sum(a.occupancy for _, 
                              a in atom.other_locations.items())

    # transfer occupancy and bfactor for heavy atoms;
    for template_residue, residue in zip(template_parm.residues, parm.residues):

      # transfer residue numbers, chainId's and insertion codes:
      residue.number = template_residue.number
      residue.chain = template_residue.chain
      residue.insertion_code = template_residue.insertion_code

      # for each residue, sort atom by name to make sure we get the same 
      #   atom order between the two parms
      template_heavy_atoms = sorted((atom for atom in template_residue.atoms 
                                   if atom.atomic_number > 1),
                                   key=lambda x : x.name)
      heavy_atoms = sorted((atom for atom in residue.atoms 
                                   if atom.atomic_number > 1),
                                   key=lambda x : x.name)
      for template_atom, atom in zip(template_heavy_atoms, heavy_atoms):
        atom.occupancy = template_atom.occupancy
        atom.bfactor = template_atom.bfactor
        atom.anisou = template_atom.anisou

    # update hydrogens by copying occupancy and bfactor from its bond partner
    for atom in parm.atoms:
      if atom.atomic_number == 1:
        atom.occupancy = atom.bond_partners[0].occupancy
        atom.bfactor = atom.bond_partners[0].bfactor

    parm.write_pdb(pdb_filename, standard_resnames=True, renumber=False,
        write_anisou=True)

  def uc_parm7_to_4phenix_pdb(self, parm7_file, rst7_file, 
                              template_pdb, outpdb):
    '''
    Combine information in {parm7,rst7} with that in the template
       pdb file to create 4phenix_xxxx.pdb; (uses parmed)
    Original application was that parm7/rst7 came from initial minimization;
       but they could come from any type of Amber run

    Note:  the atom order in outpdb will be the "Amber" atom order; use phenix
    routines after this to convert to phenix atom order.

    Note also: there is a lot of duplication between this routine and code
    in asu_parm7_to_4phenix_pdb()
    '''

    template_parm = parmed.load_file(template_pdb)
    parm = parmed.load_file( parm7_file, xyz=rst7_file )

    # truncate to unit cell:
    n_asu_residue = len(template_parm.residues)
    selection = ':1-{}'.format(n_asu_residue)
    asu_parm = parm[selection]
    asu_parm.box = parm.box
    asu_parm.space_group = parm.space_group
    asu_parm.symmetry = parm.symmetry

    # following kludge is needed to get got atom numbers; not needed if
    #   conversion to phenix order is done right after this
    #asu_parm.write_pdb( pdb_filename )
    #asu_parm = parmed.load_file( pdb_filename )

    # transfer occupancy and bfactor for heavy atoms;
    for template_residue, residue in zip(template_parm.residues, 
                                         asu_parm.residues):

      # transfer residue numbers, chainId's, and insertion codes:
      residue.number = template_residue.number
      residue.chain = template_residue.chain
      residue.insertion_code = template_residue.insertion_code

      # for each residue, sort atom by name to make sure we get the same 
      #   atom order between the two parms
      #   TODO: handle extra atoms for LES (may not be needed?)
      template_atoms = sorted((atom for atom in template_residue.atoms),
                                   key=lambda x : x.name)
      atoms = sorted((atom for atom in residue.atoms),
                                   key=lambda x : x.name)
      for template_atom, atom in zip(template_atoms, atoms):
        atom.occupancy = template_atom.occupancy
        atom.bfactor = template_atom.bfactor
        atom.anisou = template_atom.anisou

    asu_parm.write_pdb(outpdb, standard_resnames=True, renumber=False,
        write_anisou=True)

  def return_protein_chain_gaps(self):
    assert self.pdb_hierarchy
    from mmtbx import conformation_dependent_library
    gaps = []
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

  def process_ligands(self,
                      ns_names=[],
                      nproc=1,
                      prefer_input_method=None,
                      debug=False):
    assert self.pdb_hierarchy
    if nproc > 1:
      print "\n\tParallel processes not implemented\n"

    for residue_name in ns_names:

      #  use existing mol2 or lib files if present
      if os.path.isfile('%s.lib' % residue_name):
        print "%s.lib is present. Skipping elbow/antechamber run for this residue.\n" % residue_name
      elif os.path.isfile('%s.mol2' % residue_name):
        print "%s.mol2 is present. Skipping elbow/antechamber run for this residue.\n" % residue_name

      # else check if this has already been entered in the amber library:
      elif amber_library_server.is_in_components_lib(residue_name):
        print """
  Residue "%s" already in amber monomer library. Skipping elbow/antechamber
    run for this residue.
        """ % residue_name

      #  else use chemical component dictionary if requested (should be rare):
      elif prefer_input_method == "chemical_component":
        _run_antechamber_ccif(residue_name)

      #  default is to use elbow/antechamber (this will be the most common):
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
                use_glycam=False,
                verbose=False,
                ):

    def _parse_tleap_logfile(logfile):
      errors = []
      warnings = []
      fatals = []
      f = file(logfile, "rb")
      lines = f.read()
      f.close()
      warnings = {
          #"Failed to generate parameters": None,
          "WARNING: The unperturbed charge of the unit:": \
          "Note: Highly charged molecules with no indication of how they are to be neutralized might cause problems",
      }
      for line in lines.splitlines():
        if line.find("FATAL:") > -1:
          fatals.append(line)
        for warning, msg in warnings.items():
          if line.find(warning) > -1:
            errors.append(line)
            if msg:
              errors.append(" > %s" % msg)
      if errors:
        print "Errors and warnings in tleap"
        # print "  check logfile %s for explaination\n" % logfile
        for line in errors:
          print "  %s\n" % line
      if fatals:
        print "Fatal errors in tleap"
        for line in fatals:
          print line
        raise Sorry("Fatal errors in tleap: %s" % os.path.abspath(logfile))
      return fatals

    def _check_tleap_output(lines):
      stop = ["Could not open file"]
      for line in lines:
        for s in stop:
          if line.find(s) > -1:
            print_cmd('" --- Testing environment ---"', verbose=True)
            raise Sorry('tleap error : "%s"' % line)

    tleap_input_file = "%s_%s_tleap_input_run" % (self.base, output_base)
    f = file(tleap_input_file, "wb")
    f.write('logFile %s\n' % logfile)

    amber_dir = os.environ["AMBERHOME"]

    # Following should be true in AmberTools14/15:
    if(os.path.isfile(os.path.join(amber_dir,
                                   'dat',
                                   'leap',
                                   'cmd',
                                   'leaprc.ff14SB',
                                   ))):
      raise Sorry('Amber environment appears to be older than AmberTools16; quitting')

    # Now we can assume that we are dealing with AmberTools16 or later:

    if(os.path.isfile(os.path.join(amber_dir, 'dat', 'leap', 'cmd',
                                   'leaprc.phenix',))):
      f.write('source leaprc.phenix\n')
    if( redq ):
      f.write('source leaprc.ff14SB.redq\n')
    else:
      f.write('source leaprc.protein.ff14SB\n')
      f.write('source leaprc.DNA.OL15\n')
      f.write('source leaprc.RNA.OL3\n')
    if( use_glycam ):
       f.write('source leaprc.GLYCAM_06j-1\n')
    f.write('source leaprc.water.tip3p\n')
    f.write('source leaprc.gaff2\n')
    f.write('set default nocenter on\n')
    f.write('set default reorder_residues %s\n' % reorder_residues)

    #  mechanism for user modifications:
    if os.path.isfile('myleaprc'):
       f.write('source myleaprc\n')

    for res in ns_names:
      if os.path.isfile('%s.lib' % res):
        f.write('loadOff %s.lib\n' % res)
        f.write('loadAmberParams %s.frcmod\n' % res)
      elif amber_library_server.is_in_components_lib(res):
        res_path = amber_library_server.path_in_components_lib(res)
        if res_path[1].find('.lib') > 0:
           f.write('loadOff %s\n' % res_path[1])
        else:
           f.write('%s = loadMol2 %s\n' % (res, res_path[1]))
        f.write('loadAmberParams %s\n' % (res_path[0]))
      else:
        f.write('%s = loadmol2 %s.mol2\n' % (res, res))
        f.write('loadAmberParams %s.frcmod\n' % res)

    #
    # input PDB file
    #
    f.write('x = loadpdb %s\n' % input_pdb)

    assert self.cryst1
    uc = self.cryst1.unit_cell().parameters()
    f.write('set x box { %s  %s  %s }\n' % (uc[0], uc[1], uc[2]))

    #
    #  process gaplist
    #
    if gaplist:
      for d, res1, resid1, res2, resid2 in gaplist:
        f.write('deleteBond x.%d.C x.%d.N\n' % (resid1+1, resid2+1))
    #
    #  process sslist
    #
    if sslist:
      for resid1, resid2 in sslist:
        # convert from 0-based to 1-based index:
        f.write('bond x.%d.SG x.%d.SG\n' % (resid1+1, resid2+1))

    #  second mechanism for user modifications:
    if output_base == 'uc' and os.path.isfile('myuclinks'):
       f.write('source myuclinks\n')

    f.write('saveAmberParm x %s.prmtop %s.rst7\n' % (
        "%s_%s" % (self.base, output_base),
        "%s_%s" % (self.base, output_base) ) )
    f.write('quit\n')
    f.close()
 
    #
    # strangely tleap appends to the logfile so must delete first
    #
    if os.path.exists(logfile):
      os.remove(logfile)
    cmd = os.path.join( os.environ["AMBERHOME"],'bin','tleap' )
    cmd += ' -f %s' % tleap_input_file
    print_cmd(cmd)
    ero = easy_run.fully_buffered(cmd)
    _check_tleap_output(ero.stdout_lines)
    if verbose:
      ero.show_stdout()
      ero.show_stderr()
    # check output
    # fatals = _parse_tleap_logfile(logfile)
    check_required_output_filenames(["%s_%s.prmtop" % (self.base, output_base),
                                     "%s_%s.rst7" % (self.base, output_base),
                                     ],
                                    error_display_file="leap.log",
                                    )
    # remove files
    #if os.path.exists(tleap_input_file): os.remove(tleap_input_file)
    return 0

  # change box size in rst7 file
  def update_rst7_box(self, output_base):
    assert self.cryst1
    uc = self.cryst1.unit_cell().parameters()
    # note: dangerous to use same file for input and output here?
    cmd = os.path.join( os.environ["AMBERHOME"],'bin','ChBox' )
    cmd += " -c %s_%s.rst7 -o %s_%s.rst7" % (self.base,
                                                 output_base,
                                                 self.base,
                                                 output_base,
                                                 )
    cmd += " -X %s -Y %s -Z %s -al %s -bt %s -gm %s " % tuple(uc)
    print_cmd(cmd)
    ero = easy_run.fully_buffered(cmd)
    ero.show_stdout()
    ero.show_stderr()
    return 0

  def build_unitcell_prmtop_and_rst7_files(self, redq=False, use_glycam=False,
          use_amber_unitcell=False):

    #-----------------------------------------------------------------
    # Step 1: add SYMTRY/CRYST1 to 4phenix_xxxx.pdb -> xxxx_4UnitCell.pdb
    #-----------------------------------------------------------------

    uc_pdb_file = "%s_4UnitCell.pdb" % self.base
    with open(uc_pdb_file, "wb") as fout:
      with open(self.pdb_filename) as fin:

        # lines = fin.readlines()
        # cryst1card = [line for line in lines if "CRYST1" in line or "SMTRY" in line]
        # if len(cryst1card) < 1:
        #   raise Sorry("CRYST1 record required in input pdb file")
        # fout.write(cryst1card)

        for line in fin:
          if "CRYST1" in line or "SMTRY" in line:
            fout.write(line)

      with open("4phenix_%s.pdb" % self.base) as fin:
        for line in fin:
          if not "CRYST1" in line and not "ANISOU" in line:
            fout.write(line)

    #-----------------------------------------------------------------
    # Step 2: invoke either the phenix build_unitcell() method,
    #         or Amber's UnitCell program, to convert
    #         xxxx_4UnitCell.pdb to xxxx_4tleap_uc1.pdb
    #-----------------------------------------------------------------

    tleap_pdb_file1 = "%s_4tleap_uc1.pdb" % self.base
    build_unitcell(uc_pdb_file, tleap_pdb_file1,
                   use_amber_unitcell=use_amber_unitcell)

    #-----------------------------------------------------------------
    # Step 3: run xxxx_4tleap_uc1.pdb back through pdb4amber to get new 
    #         lists that describe gaps and SS bonds.  Output will be
    #         xxxx_4tleap_uc.pdb
    #         Note: don't need to call reduce this time around.
    #-----------------------------------------------------------------

    tleap_pdb_file = "%s_4tleap_uc.pdb" % self.base
    ns_names, gaplist, sslist = pdb4amber.run(
        tleap_pdb_file, tleap_pdb_file1, arg_elbow=True,
        arg_logfile=sys.stderr,
        arg_conect=False,
    )

    #-----------------------------------------------------------------
    # Step 4:  feed xxxx_4tleap_uc.pdb to tleap
    #-----------------------------------------------------------------

    # dac note: next line is currently needed for modified residues
    #   that also have modified connectivities, such as IAS in 1dy5.
    #   For now, we cannot both have modified connectivities and also
    #   check for gaps.
    #  dummy_gaplist = []
    self.run_tleap(tleap_pdb_file,
                   output_base='uc',
                   ns_names=ns_names,
                   gaplist=gaplist,
                   sslist=sslist,
                   reorder_residues='off',
                   # logfile='tleap_uc.log',
                   redq=redq,
                   use_glycam=use_glycam,
                   )
    self.update_rst7_box('uc')

    #-----------------------------------------------------------------
    #  Step 5:  rename files to "4amber_xxxx.{prmtop,rst7}
    #-----------------------------------------------------------------

    os.rename('%s_uc.rst7' % self.base, '4amber_%s.rst7' % self.base)
    os.rename('%s_uc.prmtop' % self.base, '4amber_%s.prmtop' % self.base)

  def run_minimise(self, mintype=None, minimization_options=''):
    assert self.base
    if mintype is None:
      return

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
       imin   = 1, maxcyc = 50, ncyc   = 200, ntmin  = 1, ntxo = 1,
       ntpr=10, ntr=1, restraint_wt=2.0, restraintmask='!@H=',
      /
      """
              }

    if self.LES:
      LEStype=".LES"
    else:
      LEStype=""

    if minimization_options and mintype in ["amber_h", "amber_all"]:
      inputs['amber_h'] = inputs['amber_h'].replace('/', minimization_options + '\n /')
      inputs['amber_all'] = inputs['amber_all'].replace('/', minimization_options + '\n /')

    output_rst7_file = '4amber_%s.min.rst7' % self.base
    output_pdb_file = '4phenix_%s.min.pdb' % self.base

    if mintype in ["amber_h", "amber_all"]:

      input_file = '%s_%s.in' % (self.base, mintype)
      f = open('%s_%s.in' % (self.base, mintype), 'wb')
      f.write(inputs[mintype])
      f.close()
      cmd = os.path.join( os.environ["AMBERHOME"],'bin','sander' )
      cmd += '%s -O -i %s -p %s -c %s -o %s.min.out \
           -ref %s -r %s' % (
          LEStype,
          input_file,
          self.final_prmtop_file,
          self.final_rst7_file,
          self.base, 
          self.final_rst7_file,
          output_rst7_file
      )
      print_cmd(cmd)
      test_files_exist([input_file,
                        self.final_prmtop_file,
                        self.final_rst7_file,
                        ])
      ero = easy_run.fully_buffered(cmd)
      assert (ero.return_code == 0)

      self.uc_parm7_to_4phenix_pdb( self.final_prmtop_file, output_rst7_file,
                    "4phenix_%s.pdb" % self.base ,
                    output_pdb_file )

      pdb_inp = iotbx.pdb.input(file_name=output_pdb_file )
      pdb_h = pdb_inp.construct_hierarchy(sort_atoms=True)
      pdb_h.write_pdb_file(file_name=output_pdb_file,
                           append_end=True,
                           crystal_symmetry=self.pdb_inp.crystal_symmetry(),
                           )

    return 0

  def check_special_positions(self):
    pdb_file = '4phenix_%s.pdb' % self.base
    print 'checking special positions in %s' % pdb_file
    xrs = iotbx.pdb.input(file_name=pdb_file).xray_structure_simple()
    site_symmetry_table = xrs.site_symmetry_table()
    if site_symmetry_table.n_special_positions() > 0:
      print "WARNING: The following atoms occupy special positions."
      for i in site_symmetry_table.special_position_indices():
        print "  Atom %d" % i
      print "It may be a good idea to inspect manually and remove"
      print "atoms from special positions or rerun minimization."

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
      asu.prmtop
      asu.rst7
      uc_tleap_input_run
      asu_tleap_input_run
      4UnitCell.pdb
      leap.log
      tleap_asu.log
      tleap.in
      tleap_uc.log
      sqm.pdb
      sqm.out
      sqm.in
      amber_all.in
      mdinfo
      reduce_info.log
      reduce_lesbuilder.log
      reduce.log
      addles.log
      """
    import glob

    # print "\nChecking for files to remove"
    for filename in files_to_clean.strip().split():
      for pre in ["", "%s_" % self.base]:
        for filename in glob.glob("%s%s" % (pre, filename)):
          # print '  removing' , filename
          os.remove(filename)
    for filename in glob.glob("./4antechamber*"):
      # print '  removing' , filename
      os.remove(filename)
    for filename in glob.glob("./*.pickle"):
      # print '  removing' , filename
      os.remove(filename)
    for s in ['%s.eff',
              '%s_curated.pdb',
              '%s_uc.pdb',
              '%s_uc_H.pdb',
              '%sab.rst7',
              '4amber_%s.pdb',
              '4amber_%s.LES.pdb',
              ]:
      if os.path.isfile(s % self.base):
        # print '  removing' , s % self.base
        os.remove(s % self.base)

#  end of the AmberPrepRunner class block

def get_molecule_from_hierarchy(hierarchy, resname):
  # only works for non altloc files
  # note: only called once
  from elbow.chemistry.MoleculeClass import MoleculeClass
  mol = MoleculeClass()
  for residue_group in hierarchy.residue_groups():
    for atom_group in residue_group.atom_groups():
      if atom_group.resname == resname:
        for atom in atom_group.atoms():
          mol.AddAtom(atom.element, xyz=atom.xyz)
          mol[-1].name = atom.name
        break
    if mol: break
  return mol

def _run_antechamber_ccif(residue_name,
                          use_am1_and_maxcyc_zero=True,
                          debug=False):
  '''
  run antechamber from a components.cif file:
  '''
  print "\n=================================================="
  print "Running antechamber_ccif for %s " % residue_name
  print "=================================================="

  ccif = get_chemical_components_file_name(residue_name)
  cmds = []
  cmd = os.path.join( os.environ["AMBERHOME"],'bin','antechamber' )
  cmd += ' -i %s -fi ccif -bk %s -o %s.mol2 -fo mol2 \
      -s 2 -pf y -c bcc -at gaff2' % (ccif, residue_name, residue_name)
  if use_am1_and_maxcyc_zero:
    cmd += ' -ek "qm_theory=\'AM1\',grms_tol=0.0005,scfconv=1.d-10,maxcyc=0,ndiis_attempts=700,"'
  cmds.append(cmd)

  for cmd in cmds:
    print_cmd(cmd)
    ero = easy_run.fully_buffered(cmd)
    stdo = StringIO.StringIO()
    ero.show_stdout(out=stdo)
    for line in stdo.getvalue().splitlines():
      if line.find('APS') > -1:
        print line
      if line.find('Error') > -1:
        raise Sorry(line)

  cmd = os.path.join( os.environ["AMBERHOME"],'bin','parmchk2' )
  cmd += ' -s 2 -i %s.mol2 -f mol2 -o %s.frcmod' % (residue_name,
                                                           residue_name)
  print_cmd(cmd)
  easy_run.fully_buffered(cmd)

  # should there be a check for output???

def _run_elbow_antechamber(pdb_hierarchy,
                           residue_name,
                           use_am1_and_maxcyc_zero=True,
                           debug=False):
  '''
  run elbow and antechamber
  '''
  pdb_mol = get_molecule_from_hierarchy(pdb_hierarchy, residue_name)
  names = []
  for atom in pdb_mol:
    names.append(atom.name.strip())
  pdb_set = set(names)
  print "\n=================================================="
  print "Running elbow/antechamber for %s " % residue_name
  print "=================================================="
  if debug and 0:
    import pickle
    pf = "%s.pickle" % residue_name
    if os.path.exists(pf):
      f = file(pf, "rb")
      mol = pickle.load(f)
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
                        # no_output=True,
                        # pH=8, # special Amber number...
                        silent=True)
  names = []
  for atom in mol:
    names.append(atom.name.strip())
  cc_set_h = set(names)
  names = []
  for atom in mol:
    if atom.isH():
      continue
    names.append(atom.name.strip())
  cc_set = set(names)
  if (pdb_set.symmetric_difference(cc_set) and
      pdb_set.symmetric_difference(cc_set_h)
      ):
    from elbow.utilities import molecule_superposition
    pdb_mol.Bondise(add_bond_order=False)
    rc = molecule_superposition.run(mol,
                                    pdb_mol,
                                    use_hydrogens=False,
                                    seed_with_unique_atoms=True,
                                    )
    # molecule_superposition.print_return_list(rc)
    for atom1, atom2 in rc:
      atom1.name = atom2.name

  mol.WritePDB('4antechamber_%s.pdb' % residue_name,
               pymol_pdb_bond_order=False,
               )
  mol.Multiplicitise()
  print mol.DisplayBrief()

  cmds = []
  cmd = os.path.join( os.environ["AMBERHOME"],'bin','antechamber' )
  cmd += ' -i 4antechamber_%s.pdb -fi pdb -o %s.mol2 -fo mol2 \
      -nc %d -m %d -s 2 -pf y -c bcc -at gaff2' \
      % (residue_name, residue_name, mol.charge, mol.multiplicity)
  if use_am1_and_maxcyc_zero:
    cmd += ' -ek "qm_theory=\'AM1\',grms_tol=0.0005,scfconv=1.d-10,maxcyc=0,ndiis_attempts=700,"'
  cmds.append(cmd)
  if not use_am1_and_maxcyc_zero:
    cmd = os.path.join( os.environ["AMBERHOME"],'bin','antechamber' )
    cmd += ' -i sqm.pdb -fi pdb -o %s.mol2 -fo mol2 \
      -nc %s -m %d -s 2 -pf y -c bcc -at gaff2' \
      % (residue_name, mol.charge, mol.multiplicity)
    cmds.append(cmd)

  for cmd in cmds:
    print_cmd(cmd)
    ero = easy_run.fully_buffered(cmd)
    stdo = StringIO.StringIO()
    ero.show_stdout(out=stdo)
    for line in stdo.getvalue().splitlines():
      if line.find('APS') > -1:
        print line
      if line.find('Error') > -1:
        raise Sorry(line)

  cmd = os.path.join( os.environ["AMBERHOME"],'bin','parmchk2' )
  cmd += ' -s 2 -i %s.mol2 -f mol2 -o %s.frcmod' % (residue_name, residue_name)
  print_cmd(cmd)
  easy_run.fully_buffered(cmd)

def run(rargs):
  working_params = setup_options_args(rargs)
  inputs = working_params.amber_prep.input
  actions = working_params.amber_prep.actions
  base = get_output_preamble(working_params)
  amber_prep_runner = AmberPrepRunner(base, LES=actions.LES)
  amber_prep_runner.initialize_pdb(inputs.pdb_file_name)

  # basename = os.path.basename(inputs.pdb_file_name)
  # amber_prep_runner.curate_model(remove_alt_confs=True)
  #current_pdb_file_name = basename.replace(
  #    '.pdb',
  #    '_curated.pdb',
  #)
  #amber_prep_runner.write_pdb_hierarchy(current_pdb_file_name)

  print "\n=================================================="
  print "Running pdb4amber on %s" % inputs.pdb_file_name
  print "=================================================="

  tleap_input_pdb = "%s_4tleap.pdb" % base
  # log = []
  ns_names, gaplist, sslist = pdb4amber.run(tleap_input_pdb,
                                            inputs.pdb_file_name,
                                            arg_elbow=True,
                                            arg_reduce=actions.use_reduce,
                                            arg_logfile=sys.stderr,
                                            arg_conect=False,
                                            arg_no_reduce_db=True,
                                            )

  print "\n=================================================="
  print "Setting up library files for non-standard residues"
  print "=================================================="

  amber_prep_runner.process_ligands(
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
  #   4amber_xxxx.prmtop file.
  #
  # N.B.: this means that the base_asu.prmtop file should never be used!
  #   we might want to make sure that this file is always removed.

  dummy_gaplist = []
  amber_prep_runner.run_tleap('%s_4tleap.pdb' % amber_prep_runner.base,
                              output_base='asu',
                              ns_names=ns_names,
                              gaplist=dummy_gaplist,
                              sslist=sslist,
                              reorder_residues='off',
                              # logfile='tleap_asu.log',
                              redq=actions.redq,
                              use_glycam=actions.use_glycam,
                              )
  amber_prep_runner.update_rst7_box("asu")

  phenix_file = '4phenix_%s.pdb' % amber_prep_runner.base
  amber_prep_runner.asu_parm7_to_4phenix_pdb(phenix_file)

  # N.B: above 4phenix file does not have the phenix atom order inside
  #    residues: pass this through a phenix hierarchy function to get that.
  pdb_post = iotbx.pdb.input(file_name=phenix_file)
  pdb_h_post = pdb_post.construct_hierarchy(sort_atoms=True)
  pdb_h_post.write_pdb_file(file_name=phenix_file,
            append_end=True,
            crystal_symmetry=amber_prep_runner.pdb_inp.crystal_symmetry() )

  print "\n============================================================"
  print "Preparing unit cell files: 4amber_%s.prmtop and 4amber_%s.rst7" % (base, base)
  print "============================================================"

  amber_prep_runner.build_unitcell_prmtop_and_rst7_files(redq=actions.redq,
          use_glycam=actions.use_glycam,
          use_amber_unitcell=actions.use_amber_unitcell)

  if actions.LES:
    print "\n=================================================="
    print "Building the LES prmtop and rst7 files"
    print "=================================================="
    pdb_file_name = working_params.amber_prep.input.pdb_file_name
    les_builder = LESBuilder(
        pdb_file_name,
        prmtop=amber_prep_runner.non_les_prmtop_file_name,
        rst7_file=amber_prep_runner.non_les_rst7_file_name,
        addles_input_file=actions.addles_input)
    les_builder.run(use_amber_unitcell=actions.use_amber_unitcell,
                    use_reduce=actions.use_reduce)

    #  rename the files to the canonical three we want.  (Do this here,
    #    so that these lines can be commented out if debugging is 
    #    required.)
    os.rename('4amber_%s.LES.prmtop' % base,'4amber_%s.prmtop' % base)
    os.rename('4amber_%s.LES.rst7' % base, '4amber_%s.rst7' % base)
    os.rename('4phenix_%s.LES.pdb' % base, '4phenix_%s.pdb' % base)

  #
  # get atom order file
  #
  from iotbx import pdb
  pdb_inp =  pdb.input('4phenix_%s.pdb' % base)
  hierarchy = pdb_inp.construct_hierarchy()
  sites_cart = hierarchy.atoms().extract_xyz()
  parm = parmed.load_file("4amber_%s.prmtop" % base, 
                          xyz="4amber_%s.rst7" % base )
  print parm
  asu_n_atoms = len(sites_cart)
  n_models = int(len(parm.coordinates) / asu_n_atoms)
  order_converter = get_indices_convert_dict_from_array(
    sites_cart,
    parm.coordinates[:asu_n_atoms])
  import numpy as np
  # asu
  asu_p2a_indices = order_converter['p2a']
  asu_a2p_indices = order_converter['a2p']

  # unitcell
  uc_p2a_indices = []
  uc_a2p_indices = []

  for index in range(n_models):
    # need to increase atom indices
    uc_p2a_indices.extend((asu_p2a_indices + index * asu_n_atoms).tolist())
    uc_a2p_indices.extend((asu_a2p_indices + index * asu_n_atoms).tolist())

  # extend array for unitcell
  order_converter['p2a'] = np.array(uc_p2a_indices)
  order_converter['a2p'] = np.array(uc_a2p_indices)

  # save to disk for debugging
  saved_arr = np.array([order_converter['a2p'],
                        order_converter['p2a']],
                        dtype='i4')
  # 1st column: amber -> phenix
  # 2nd column: phenix -> amber
  filename = amber_prep_runner.final_order_file
  np.savetxt(filename, saved_arr.transpose(), fmt='%5d')

  if actions.minimise != "off":
    print "\n=================================================="
    print "Minimizing input coordinates."
    print "=================================================="
    amber_prep_runner.run_minimise(actions.minimise,
                      minimization_options=actions.minimization_options)

    #  rename the files to the canonical three we want.  (Do this here,
    #    so that these lines can be commented out if debugging is 
    #    required.)
    os.rename('4amber_%s.min.rst7' % base, '4amber_%s.rst7' % base)
    os.rename('4phenix_%s.min.pdb' % base, '4phenix_%s.pdb' % base)

  amber_prep_runner.check_special_positions()

  if actions.clean:
    amber_prep_runner.run_clean()

  outl = "\n==================================================\n"
  outl += "Done.  Four new files have been made:\n"
  outl += "      %s\n" % amber_prep_runner.final_pdb_file_4phenix
  outl += "      %s\n" % amber_prep_runner.final_rst7_file
  outl += "      %s\n" % amber_prep_runner.final_prmtop_file
  outl += "      %s\n" % amber_prep_runner.final_order_file
  outl += "==================================================\n\n"
  outl += "Example\n\n  phenix.refine"
  outl += " %s use_amber=True \\\n" % (
      amber_prep_runner.final_pdb_file_4phenix,
  )
  outl += "    amber.topology_file_name=%s \\\n" % amber_prep_runner.final_prmtop_file
  outl += "    amber.coordinate_file_name=%s \\\n" % amber_prep_runner.final_rst7_file
  outl += "    amber.order_file_name=%s \\\n" % amber_prep_runner.final_order_file
  outl += "    ....(other refinement keywords here)....."
  outl += "\n\n\n"
  print outl
  return amber_prep_runner.final_pdb_file_4phenix

if __name__ == "__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(args)
