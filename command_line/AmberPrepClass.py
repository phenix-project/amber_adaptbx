# LIBTBX_SET_DISPATCHER_NAME phenix.AmberPrep

import os
import sys
from datetime import datetime
import iotbx.pdb
import StringIO
from libtbx import phil
from libtbx.utils import Sorry
import libtbx.load_env
import libtbx.phil.command_line
from libtbx import easy_run
from elbow.command_line import builder
from amber_adaptbx import pdb4amber
from amber_adaptbx import amber_library_server
from amber_adaptbx.utils import build_unitcell, write_standard_pdb
from amber_adaptbx.les_builder.build import LESBuilder
import parmed as pmd

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
      minimization_options = ''
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
      addles_input = ''
        .type = str
        .caption = User specify addles input filename. Optional.
        .help = User specify addles input filename. Optional.
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
      if iotbx.pdb.is_pdb_file(arg):
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
  print "\n~> %s\n" % cmd
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
    # TODO: remove unused variables
    self.base = base_name
    self.LES = LES
    self.final_prmtop_file = '4amber_%s.prmtop' % self.base
    self.final_rst7_file = '4amber_%s.rst7' % self.base
    self.final_pdb_file_4phenix = '4phenix_%s.pdb' % self.base
    if self.LES:
      self.final_prmtop_file = self.final_prmtop_file.replace('.prmtop', '.LES.prmtop')
      self.final_rst7_file = self.final_rst7_file.replace('.rst7', '.LES.rst7')
      self.final_pdb_file_4phenix = self.final_pdb_file_4phenix.replace('.pdb', '.LES.pdb')
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

  def curate_model(self, remove_alt_confs=False):
    assert self.pdb_hierarchy
    from mmtbx import pdbtools
    if remove_alt_confs:
      # performs removal in place
      pdbtools.remove_alt_confs(self.pdb_hierarchy,
                                always_keep_one_conformer=True,
                                )
  def correct_atom_occupancy_and_bfactor(self, pdb_filename):
    # from Dave:
    # Here's the goal: We want each atom in the 4phenix_xxxx.pdb file to have an
    # occupancy that is the sum of the occupancies of the all the alternate
    # conformers in the original pdb file.
    # 
    # If there was only one conformer in the original, we should keep that
    # occupancy, whatever it is.  If there was more than one conformer, sum the
    # occupancies of the various conformers.  We should do this on an atom-by-atom
    # basis.  This should make the total number of electrons in the 4phenix_xxxx.pdb
    # file the same as in the original pdb file.
    # 
    # I think/hope(?) parmed keeps information about the occupancies of multiple
    # conformers, so all one needs to do is to sum them up.
    # end Dave commend

    # reduce program does not add H for water
    # tleap will do that, so we need to update this information.
    # also update H-occupancy for other residues too.
    # b-factor will be also updated
    # this is used for non-LES case
    # load very original pdb file as info source

    template_parm = pmd.load_file(self.pdb_filename)
    parm = pmd.load_file(pdb_filename)

    # update occupancy for alt-atom of template_parm
    # so we can copy to parm
    for atom in template_parm.atoms:
      if atom.other_locations:
        # if atom has alternative locations, ParmEd save them to other_locations (dict)
        atom.occupancy += sum(a.occupancy for _, a in atom.other_locations.items())
    # copy occupancy and bfactor for heavy atoms (original pdb file does not have H)
    for template_residue, residue in zip(template_parm.residues, parm.residues):
      # for each residue, sort atom by name to make sure we get the same atom order between
      # two parms
      template_heay_atoms = sorted((atom for atom in template_residue.atoms if atom.atomic_number > 1),
                                   key=lambda x : x.name)
      heay_atoms = sorted((atom for atom in residue.atoms if atom.atomic_number > 1),
                                   key=lambda x : x.name)
      for template_atom, atom in zip(template_heay_atoms, heay_atoms):
        atom.occupancy = template_atom.occupancy
        atom.bfactor = template_atom.bfactor

    # now updating hydrogens by copying number from its bond partner
    for atom in parm.atoms:
      if atom.atomic_number == 1:
        atom.occupancy = atom.bond_partners[0].occupancy
        atom.bfactor = atom.bond_partners[0].bfactor
    write_standard_pdb(parm, pdb_filename)

  def validate_pdb(self):
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

  def run_elbow_antechamber(self,
                            ns_names=[],
                            nproc=1,
                            prefer_input_method=None,
                            debug=False):
    assert self.pdb_hierarchy
    if nproc > 1:
      print "\n\tParallel processes not implemented\n"
    for residue_name in ns_names:
      if prefer_input_method:
        if prefer_input_method == "chemical_component":
          _run_antechamber_ccif(residue_name)
          continue
        elif prefer_input_method == "elbow":
          _run_elbow_antechamber(self.pdb_hierarchy, residue_name, debug=debug)
          continue
      if amber_library_server.is_in_components_lib(residue_name):
        print """
  Residue "%s" already in amber monomer library. Skipping elbow/antechamber
    run for this residue.
        """ % residue_name
        continue
      elif os.path.isfile('%s.mol2' % residue_name):
        print "%s.mol2 is present. Skipping elbow/antechamber run for this residue.\n" % residue_name
        continue
      elif get_chemical_components_file_name(residue_name):
        _run_antechamber_ccif(residue_name)
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

    #amber_dir = libtbx.env.dist_path("amber")
    #if os.environ["AMBERHOME"]!=amber_dir:
    #  raise Sorry("$AMBERHOME %s\nnot pointing to Phenix module %s" % (
    #    os.environ["AMBERHOME"],
    #    amber_dir,
    #    ))
    amber_dir = os.environ["AMBERHOME"]

    # Following should be true in AmberTools14/15:
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
    # f.write('source leaprc.GLYCAM_06j-1\n') #un-comment for glycoproteins
    f.write('source leaprc.water.tip3p\n')
    f.write('source leaprc.gaff2\n')
    #  (for the future: have some mechanism for modifying the above list)
    f.write('set default nocenter on\n')
    f.write('set default reorder_residues %s\n' % reorder_residues)

    for res in ns_names:
      if amber_library_server.is_in_components_lib(res):
        res_path = amber_library_server.path_in_components_lib(res)
        f.write('%s = loadmol2 %s\n' % (res, res_path[1]))
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
        f.write('deleteBond x.%d.C x.%d.N\n' % (resid1, resid2))
    #
    #  process sslist
    #
    if sslist:
      for resid1, resid2 in sslist:
        f.write('bond x.%d.SG x.%d.SG\n' % (resid1, resid2))

    f.write('saveAmberParm x %s.prmtop %s.rst7\n' % (
        "%s_%s" % (self.base, output_base),
        "%s_%s" % (self.base, output_base),
    )
    )
    f.write('quit\n')
    f.close()
    #
    # strangely tleap appends to the logfile so must delete first
    #
    if os.path.exists(logfile):
      os.remove(logfile)
    cmd = 'tleap -f %s' % tleap_input_file
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
    cmd = "ChBox -c %s_%s.rst7 -o %s_%s.rst7" % (self.base,
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

  def _pdb_hierarchy_and_rename_wat(self, filename):
    pdb_inp = iotbx.pdb.input(file_name=filename)
    pdb_hierarchy = pdb_inp.construct_hierarchy(sort_atoms=True)
    # the -bres option in ambpdb does not (yet) change "WAT" to "HOH"
    for atom_group in pdb_hierarchy.atom_groups():
      if atom_group.resname in ["WAT"]:
        atom_group.resname = "HOH"
    return pdb_inp, pdb_hierarchy

  def _match_hierarchies_and_transfer_to(self,
                                         hierachy1,  # pre
                                         hierachy2,  # post
                                         transfer_b=False,
                                         transfer_occ=False,
                                         transfer_chain_id=False,
                                         transfer_xyz=False,
                                         ):
    # match residues based on resseq and resname
    # match atoms based on name an i_seq

    # dac note, 12/16: this routine is impossibly slow!  I have re-ordered
    #   the loops for some speedup, but we need to re-think this: I think
    #   we can assume that the chains and residues (but not the atoms) are
    #   in the same order in both hierarchies.

    for chain_pre in hierachy1.chains():
      for chain_post in hierachy2.chains():
        if transfer_chain_id:
          chain_post.id = chain_pre.id

        for resi_pre in chain_pre.conformers()[0].residues():

          # TODO: put this into a function?:
          # convert Amber residue names to Brookhaven standards:
          #  (only needed here for the "pre" hierarchy)
          pre_resname = resi_pre.resname.strip()
          if pre_resname == "CYX": pre_resname = "CYS"
          if pre_resname == "HID": pre_resname = "HIS"
          if pre_resname == "HIE": pre_resname = "HIS"
          if pre_resname == "HIP": pre_resname = "HIS"
          if pre_resname == "C3":  pre_resname = "C"
          if pre_resname == "U3":  pre_resname = "U"
          if pre_resname == "G3":  pre_resname = "G"
          if pre_resname == "A3":  pre_resname = "A"
          if pre_resname == "C5":  pre_resname = "C"
          if pre_resname == "U5":  pre_resname = "U"
          if pre_resname == "G5":  pre_resname = "G"
          if pre_resname == "A5":  pre_resname = "A"
          if pre_resname == "DC3":  pre_resname = "DC"
          if pre_resname == "DT3":  pre_resname = "DT"
          if pre_resname == "DG3":  pre_resname = "DG"
          if pre_resname == "DA3":  pre_resname = "DA"
          if pre_resname == "DC5":  pre_resname = "DC"
          if pre_resname == "DT5":  pre_resname = "DT"
          if pre_resname == "DG5":  pre_resname = "DG"
          if pre_resname == "DA5":  pre_resname = "DA"

          for resi_post in chain_post.conformers()[0].residues():
            if (resi_pre.resseq == resi_post.resseq and
                    pre_resname == resi_post.resname.strip()
                ):
              for atom_pre in resi_pre.atoms():
                for atom_post in resi_post.atoms():
                  if atom_pre.name == atom_post.name:
                    if transfer_b:
                      atom_post.b = atom_pre.b
                    if transfer_occ:
                      atom_post.occ = atom_pre.occ
                    if transfer_xyz:
                      atom_post.xyz = (atom_pre.xyz[0],
                                       atom_pre.xyz[1],
                                       atom_pre.xyz[2])

  #--------------------------------------------------------------------------
  # make an Amber-compatible asu-only pdb file, and transfer occupancies,
  #    b-factors and chain-ids to it
  # (Only called once: takes xxxx_asu.{prmtop,rst7} as inputs, writes
  #    xxxx_new2.pdb.  Also creates intermediate file xxxx_new.pdb
  #--------------------------------------------------------------------------
  def run_ambpdb_and_transfer(self):
    assert self.base
    cmd = 'ambpdb -bres -p %s_asu.prmtop < %s_asu.rst7 > %s_new.pdb' % tuple(
        [self.base] * 3
    )
    print_cmd(cmd)
    ero = easy_run.fully_buffered(cmd)
    ero.show_stdout()
    ero.show_stderr()

    pdb_pre, pdb_h_pre = self._pdb_hierarchy_and_rename_wat('%s_4tleap.pdb' % self.base)
    pdb_post, pdb_h_post = self._pdb_hierarchy_and_rename_wat('%s_new.pdb' % self.base)

    print "starting match_hierarchies: %s\n" % str(datetime.now())
    self._match_hierarchies_and_transfer_to(pdb_h_pre,  # from
                                            pdb_h_post,  # to
                                            transfer_b=True,
                                            transfer_occ=True,
                                            transfer_chain_id=True,
                                            )

    print "starting write_pdb: %s\n" % str(datetime.now())
    pdb_h_post.write_pdb_file(file_name='%s_new2.pdb' % self.base,
                              append_end=True,
                              crystal_symmetry=pdb_pre.crystal_symmetry(),
                              )
    return 0

  # add cryst1 and sscale; creates 4phenix_base_type_.pdb file
  def finalize_pdb(self,
                   pdb_filename=None, sort_atoms=True, type='',
                   ):
    # should we change 'type' to something else. This is Python keyword
    # Please update doc for this method.
    assert self.pdb_hierarchy
    assert self.cryst1
    assert self.base
    if pdb_filename:
      pdb_inp = iotbx.pdb.input(pdb_filename)
      pdb_hierarchy = pdb_inp.construct_hierarchy(sort_atoms=sort_atoms)
    else:
      # this returns the altloc to the model so not good!!!
      pdb_hierarchy = self.pdb_hierarchy
      assert 0
    pdbstring = pdb_hierarchy.as_pdb_string(crystal_symmetry=self.cryst1)
    print('--> final_pdb_file_4phenix', self.final_pdb_file_4phenix)
    print 'Writing 4phenix file', self.final_pdb_file_4phenix
    output_filename = '4phenix_%s.pdb' % self.base
    with open(output_filename, 'w') as f:
      f.write(pdbstring)
    if type:
      # update filename after doing minimization
      # need to make a copy after doing minimization
      # TODO: should do this in another place?
      self.final_pdb_file_4phenix = '4phenix_' + self.base + type + '.pdb'
      easy_run.fully_buffered('cp %s %s' % (output_filename, self.final_pdb_file_4phenix))
    #~ import code; code.interact(local=locals())
    return 0

  def build_unitcell_prmtop_and_rst7_files(self, redq=False):

    #-----------------------------------------------------------------
    # Step 1: add SYMTRY/CRYST1 to 4phenix_xxxx.pdb -> xxxx_4UnitCell.pdb
    #-----------------------------------------------------------------

    assert self.base, 'must provide base name'
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
          if not "CRYST1" in line:
            fout.write(line)

    #-----------------------------------------------------------------
    # Step 2: invoke the phenix build_unitcell() method to convert
    #         xxxx_4UnitCell.pdb to xxxx_4tleap_uc1.pdb
    #-----------------------------------------------------------------

    tleap_pdb_file1 = "%s_4tleap_uc1.pdb" % self.base
    build_unitcell(uc_pdb_file, tleap_pdb_file1)

    #-----------------------------------------------------------------
    # Step 3: run xxxx_4leap_uc1.pdb back through pdb4amber to get new 
    #         sslist that describes SS bonds.  Output will be
    #         xxxx_4tleap_uc.pdb
    #-----------------------------------------------------------------

    tleap_pdb_file = "%s_4tleap_uc.pdb" % self.base
    ns_names, gaplist, sslist = pdb4amber.run(
        tleap_pdb_file, tleap_pdb_file1, arg_elbow=True,
    )

    #-----------------------------------------------------------------
    # Step 4:  feed xxxx_4tleap_uc.pdb to tleap
    #-----------------------------------------------------------------

    self.run_tleap(tleap_pdb_file,
                   output_base='uc',
                   ns_names=ns_names,
                   gaplist=gaplist,
                   sslist=sslist,
                   reorder_residues='off',
                   # logfile='tleap_uc.log',
                   redq=redq,
                   )
    self.update_rst7_box('uc')

    #-----------------------------------------------------------------
    #  Step 5:  rename files to "4amber_xxxx.{prmtop,rst7}
    #-----------------------------------------------------------------

    os.rename('%s_uc.rst7' % self.base, '4amber_%s.rst7' % self.base)
    os.rename('%s_uc.prmtop' % self.base, '4amber_%s.prmtop' % self.base)

  @classmethod
  def correct_resid(cls, template_pdb_file, output_file):
    ''' ensure output_file has the same resnum, chain as `template_pdb_file`

    `output_file` will be overwriten. Make sure that two pdb files
    have the same residue order
    '''
    template_parm = pmd.load_file(template_pdb_file)
    target_parm = pmd.load_file(output_file)

    for template_residue, target_residue in zip(template_parm.residues, target_parm.residues):
      target_residue.number = template_residue.number
      target_residue.chain = template_residue.chain
    write_standard_pdb(target_parm, output_file)

  def _write_LES_pdb_4phenix(self):
    # TODO: update ocupancy
    asu_parm = pmd.load_file(self.pdb_filename)
    n_asu_residue = len(asu_parm.residues)
    selection = ':1-{}'.format(n_asu_residue)

    parm = pmd.load_file(self.final_prmtop_file, xyz=self.final_rst7_file)
    # strip atoms
    asu_new_parm = parm[selection]
    asu_new_parm.box = parm.box
    asu_new_parm.space_group = asu_parm.space_group
    asu_new_parm.symmetry = asu_parm.symmetry
    # update occupancy
    # we can use original occupancy from ASU pdb file too.
    for atom in asu_new_parm.atoms:
      atom.occupancy = 1.0
    print('--> final_pdb_file_4phenix', self.final_pdb_file_4phenix)
    write_standard_pdb(asu_new_parm, self.final_pdb_file_4phenix)

  def run_minimise(self, minimization_type=None, minimization_options=''):
    assert self.base
    if minimization_type is None:
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
    if minimization_options and minimization_type in ["amber_h", "amber_all"]:
      inputs['amber_h'] = inputs['amber_h'].replace('/', minimization_options + '\n /')
      inputs['amber_all'] = inputs['amber_all'].replace('/', minimization_options + '\n /')
    output_rst7_file_name = ' %s_%s.rst7' % (self.base, minimization_type)
    if self.LES:
      output_rst7_file_name = output_rst7_file_name.replace('.rst7', '.LES.rst7')

    if minimization_type in ["amber_h", "amber_all"]:

      input_file = '%s_%s.in' % (self.base, minimization_type)
      f = open('%s_%s.in' % (self.base, minimization_type), 'wb')
      f.write(inputs[minimization_type])
      f.close()
      cmd = 'sander -O -i %s -p %s -c %s -o %s_%s.out \
           -ref %s -r %s' % (
          input_file,
          self.final_prmtop_file,
          self.final_rst7_file,
          self.base,
          minimization_type,
          self.final_rst7_file,
          output_rst7_file_name
      )
      if self.LES:
        cmd = cmd.replace('sander', 'sander.LES')
      print_cmd(cmd)
      # test function that may be useful...
      test_files_exist([input_file,
                        self.final_prmtop_file,
                        self.final_rst7_file,
                        ])
      ero = easy_run.fully_buffered(cmd)
      assert (ero.return_code == 0)
      test_files_exist([input_file,
                        self.final_prmtop_file,
                        self.final_rst7_file,
                        ])

      cmd = 'ambpdb -bres -p %s < %s > %s_new.pdb' % (
          self.final_prmtop_file,
          output_rst7_file_name,
          self.base,
      )
      print_cmd(cmd)
      ero = easy_run.fully_buffered(cmd)
      assert (ero.return_code == 0)
      ero.show_stdout()
      ero.show_stderr()
      # rename
      if self.LES:
        self.final_rst7_file = '4amber_' + self.base + '.LES.min.{}.rst7'.format(minimization_type)
      else:
        self.final_rst7_file = '4amber_' + self.base + '.min.{}.rst7'.format(minimization_type)
      # TODO: Why I can not use os.rename? (Got OSError about file not found. Weird)
      # save minimized rst7
      easy_run.fully_buffered('cp {} {}'.format(output_rst7_file_name, self.final_rst7_file))

      if self.LES:
        # TODO: remove this and use Nigel's version.
        #    dac: what do you mean by "Nigel's version"???
        # Why using this right now? Seems too me that the output pdb from Nigel's code
        # does is not a reordered version of minimized rst7 file. Or may be I made a bug.
        # we should avoid writing too many files to disk.
        self.final_pdb_file_4phenix = '4phenix_%s.LES.min.%s.pdb' % (self.base, minimization_type)
        self._write_LES_pdb_4phenix()
      else:
        pdb_pre, pdb_h_pre = self._pdb_hierarchy_and_rename_wat('4phenix_%s.pdb' % self.base)
        pdb_post, pdb_h_post = self._pdb_hierarchy_and_rename_wat('%s_new.pdb' % self.base)

        # there is a function that will transfer the coordinates from one PDB
        # hierarchy to another but the atoms have to be the same number & order
        self._match_hierarchies_and_transfer_to(pdb_h_post,  # from
                                                pdb_h_pre,  # to
                                                transfer_xyz=True,
                                                )

        pdb_h_pre.write_pdb_file(file_name='%s_new2.pdb' % self.base,
                                 append_end=True,
                                 crystal_symmetry=pdb_pre.crystal_symmetry(),
                                 )

        type_ = '.min.%s' % minimization_type
        self.finalize_pdb(pdb_filename='%s_new2.pdb' % self.base,
                          sort_atoms=True, type=type_)
    elif minimization_type == "phenix_all":
      if self.LES:
        cmd = 'phenix.geometry_minimization 4phenix_%s.LES.pdb amber.use_amber=True \
             amber.topology_file_name=%s \
             amber.coordinate_file_name=%s  \
             output_file_name_prefix=4phenix_%s_minPhenix' % (self.base, self.final_prmtop_file, self.final_rst7_file, self.base)
      else:
        cmd = 'phenix.geometry_minimization 4phenix_%s.pdb amber.use_amber=True \
             amber.topology_file_name=4amber_%s.prmtop \
             amber.coordinate_file_name=4amber_%s.rst7  \
             output_file_name_prefix=4phenix_%s_minPhenix ' % tuple([self.base] * 4)

      cmd = cmd + ' ' + minimization_options
      restraints = "%s.ligands.cif" % self.base
      if os.path.exists(restraints):
        cmd += " %s" % restraints
      print_cmd(cmd)
      ero = easy_run.fully_buffered(cmd).raise_if_errors()  # .return_code
      assert (ero.return_code == 0)
      ero.show_stdout()
      ero.show_stderr()
      if self.LES:
        self.final_pdb_file_4phenix = '4phenix_%s.LES.min.%s.pdb' % (self.base, minimization_type)
      else:
        self.final_pdb_file_4phenix = '4phenix_%s.min.%s.pdb' % (self.base, minimization_type)
      os.rename('4phenix_%s_minPhenix.pdb' % self.base,
                self.final_pdb_file_4phenix)
    return 0

  def check_special_positions(self):
    if self.LES:
      pdb_file = '4phenix_%s.LES.pdb' % self.base
    else:
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
      addles.in
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
    for s in ['%s.eff',
              '%s_curated.pdb',
              '%s_uc.pdb',
              '%s_uc_H.pdb',
              '%sab.rst7',
              '4amber_%s.pdb',
              '4amber_%s.LES.pdb',
              ]:
      if os.path.isfile(s % self.base):
        os.remove(s % self.base)

  def write_pdb_hierarchy(self, pdb_file_name):
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
      if atom_group.resname == resname:
        for atom in atom_group.atoms():
          mol.AddAtom(atom.element, xyz=atom.xyz)
          mol[-1].name = atom.name
        break
  return mol

# run antechamber from a components.cif file:


def _run_antechamber_ccif(residue_name,
                          use_am1_and_maxcyc_zero=True,
                          debug=False):
  print "\n=================================================="
  print "Running antechamber_ccif for %s " % residue_name
  print "=================================================="

  ccif = get_chemical_components_file_name(residue_name)
  cmds = []
  cmd = 'antechamber -i %s -fi ccif -bk %s -o %s.mol2 -fo mol2 \
      -s 2 -pf y -c bcc -at gaff2' % (ccif, residue_name, residue_name)
  if use_am1_and_maxcyc_zero:
    cmd += ' -ek "qm_theory=\'AM1\', grms_tol=0.0005, scfconv=1.d-10, maxcyc=0, ndiis_attempts=700,"'
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

  cmd = 'parmchk2 -s 2 -i %s.mol2 -f mol2 -o %s.frcmod' % (residue_name,
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
  cmd = 'antechamber -i 4antechamber_%s.pdb -fi pdb -o %s.mol2 -fo mol2 \
      -nc %d -m %d -s 2 -pf y -c bcc -at gaff2' \
      % (residue_name, residue_name, mol.charge, mol.multiplicity)
  if use_am1_and_maxcyc_zero:
    cmd += ' -ek "qm_theory=\'AM1\', grms_tol=0.0005, scfconv=1.d-10, maxcyc=0, ndiis_attempts=700,"'
  cmds.append(cmd)
  if not use_am1_and_maxcyc_zero:
    cmd = 'antechamber -i sqm.pdb -fi pdb -o %s.mol2 -fo mol2 \
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

  cmd = 'parmchk2 -s 2 -i %s.mol2 -f mol2 -o %s.frcmod' % (residue_name, residue_name)
  print_cmd(cmd)
  easy_run.fully_buffered(cmd)


def run(rargs):
  working_params = setup_options_args(rargs)
  inputs = working_params.amber_prep.input
  actions = working_params.amber_prep.actions
  base = get_output_preamble(working_params)
  amber_prep_runner = AmberPrepRunner(base, LES=actions.LES)
  amber_prep_runner.initialize_pdb(inputs.pdb_file_name)
  invalid = amber_prep_runner.validate_pdb()
  if invalid:
    raise Sorry('PDB input is not "valid"')
  # amber_prep_runner.curate_model(remove_alt_confs=True)
  # need to write PDB for some of the other methods
  basename = os.path.basename(inputs.pdb_file_name)
  current_pdb_file_name = basename.replace(
      '.pdb',
      '_curated.pdb',
  )
  amber_prep_runner.write_pdb_hierarchy(current_pdb_file_name)

  print "\n=================================================="
  print "Running pdb4amber on %s" % inputs.pdb_file_name
  print "=================================================="

  tleap_input_pdb = "%s_4tleap.pdb" % base
  # log = []
  ns_names, gaplist, sslist = pdb4amber.run(tleap_input_pdb,
                                            current_pdb_file_name,
                                            arg_elbow=True,
                                            arg_reduce=actions.use_reduce,
                                            # log=log,
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
  #   we might want to make sure that this file is always removed.

  dummy_gaplist = []
  amber_prep_runner.run_tleap('%s_4tleap.pdb' % amber_prep_runner.base,
                              'asu',
                              ns_names,
                              dummy_gaplist,
                              sslist,
                              reorder_residues='off',
                              # logfile='tleap_asu.log',
                              redq=actions.redq,
                              )
  amber_prep_runner.update_rst7_box("asu")
  amber_prep_runner.run_ambpdb_and_transfer()   # (note: only called once)
  amber_prep_runner.finalize_pdb(
      pdb_filename='%s_new2.pdb' % amber_prep_runner.base, sort_atoms=True)

  print "\n============================================================"
  print "Preparing unit cell files: 4amber_%s.prmtop and 4amber_%s.rst7" % (base, base)
  print "============================================================"

  amber_prep_runner.build_unitcell_prmtop_and_rst7_files(redq=actions.redq)
  inout = '4phenix_%s.pdb' % amber_prep_runner.base
  amber_prep_runner.correct_resid(amber_prep_runner.pdb_filename, inout)
  amber_prep_runner.correct_atom_occupancy_and_bfactor(inout)
  if actions.LES:
    pdb_file_name = working_params.amber_prep.input.pdb_file_name
    les_builder = LESBuilder(
        pdb_file_name,
        prmtop=amber_prep_runner.non_les_prmtop_file_name,
        rst7_file=amber_prep_runner.non_les_rst7_file_name,
        addles_input_file=actions.addles_input)
    les_builder.run()

  if actions.minimise == "off":
    pass
  else:
    print "\n=================================================="
    print "Minimizing input coordinates."
    print "=================================================="
    amber_prep_runner.run_minimise(actions.minimise,
                                   minimization_options=actions.minimization_options)

  amber_prep_runner.check_special_positions()

  if actions.clean:
    amber_prep_runner.run_clean()

  outl = "\n==================================================\n"
  outl += "Done.  Three new files have been made:\n"
  outl += "      %s\n" % amber_prep_runner.final_pdb_file_4phenix
  outl += "      %s\n" % amber_prep_runner.final_prmtop_file
  outl += "      %s\n" % amber_prep_runner.final_rst7_file
  outl += "==================================================\n\n"
  outl += "Example\n\n  phenix.refine"
  outl += " %s use_amber=True \\\n" % (
      amber_prep_runner.final_pdb_file_4phenix,
  )
  outl += "    amber.topology_file_name=%s \\\n" % amber_prep_runner.final_prmtop_file
  outl += "    amber.coordinate_file_name=%s \\\n" % amber_prep_runner.final_rst7_file
  outl += "    ....(other refinement keywords here)....."
  outl += "\n\n\n"
  print outl
  return amber_prep_runner.final_pdb_file_4phenix

if __name__ == "__main__":
  if 1:
    args = sys.argv[1:]
    del sys.argv[1:]
    run(args)
  else:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_file_name", help="name of pdb file")
    parser.add_argument("min", help="option to minimize", default=0)
    parser.add_argument("-c", "--no_clean", help="don't remove "
                        "intermediate files", action='True', default=False)
    args = parser.parse_args()
    run(args.pdb_file_name, minimize=args.min, clean=args.no_clean)
    AmberPrepRunner.run_minimise(args.minimise, minimization_options=args.minimization_options)
