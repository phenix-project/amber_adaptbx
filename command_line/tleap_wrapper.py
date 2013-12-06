# LIBTBX_SET_DISPATCHER_NAME elbow.amber.tleap
import os, sys
import StringIO

import iotbx
from iotbx import pdb

from libtbx import easy_run
from libtbx.utils import Sorry

from libtbx import phil
import libtbx.phil.command_line
from libtbx.option_parser import OptionParser
#from libtbx import runtime_utils

his_names = ["HIS", "HID", "HIE", "HIP"]

changing_residues = []
changing_residues += his_names

master_phil_string = """

amber_tleap
  .caption = None
{
  input
  {
    pdb_file_name = None
      .type = path
      .short_caption = model
      .help = PDB filename
      .style = bold file_type:pdb
  }
  output
  {
    pdb_file_name = None
      .type = path
  }
}
"""
master_params = master_phil_string # need for auto documentation
if False:
  print '-'*80
  print master_phil_string
  print '-'*80
master_phil = phil.parse(master_phil_string)

def setup_parser():
  parser = OptionParser(
    prog="phenix.carbo_load",
    version="""
  up-to-date version
""",
    usage="""
  phenix.carbo_load protein_pdb_file_name=pdb3a37.ent \\ 
                    carbohydrate_file_name=man.txt \\
                    map_coeffs_file_name=pdb3a36_refine_001_map_coeffs.mtz \\
                    residue_selection="resname LG1"
""",
    )
  # Input options
  parser.add_option("",
                    "--show_defaults",
                    dest="show_defaults",
                    default=False,
                    action="store_true",
                    help="Display defaults",
                    )
  if 0:
    parser.add_option("",
                      "--verbose",
                      dest="verbose",
                      default=False,
                      action="store_true",
                      help="Verbose output",
                      )
    parser.add_option("",
                      "--silent",
                      dest="silent",
                      default=False,
                      action="store_true",
                      help="No output to screen",
                      )
    parser.add_option("",
                      "--dry-run",
                      dest="dry_run",
                      default=False,
                      action="store_true",
                      help="Display which residues will be processed",
                      )
  return parser

def write_tleap_cmd_file(pdb_filename, tleap_input="tleap.in"):
  cmd_template = """
source leaprc.ff12SB
loadamberparams frcmod.ionsjc_tip3p
loadAmberParams frcmod.tip3pf


x = loadPdb "%s"
"""
  loads = """
loadAmberPrep nitrate.prepin
loadAmberParams nitrate.frcmod
loadAmberPrep acetate.prepin
"""
  bond_line_template = "bond x.%s.%s x.%s.%s"
  set_box_template = "set x box {%0.3f %0.3f %0.3f}"
  save_and_quit = """
set default nocenter on
saveAmberParm x %s.prmtop %s.rst7
quit
"""
  cmd = cmd_template % (pdb_filename)
  sulfur_bonds = []
  for sb in sulfur_bonds:
    print sb
    cmd += bond_line_template(1,"SG", 2, "FE")
  cmd += set_box_template % (99, 100, 101)
  cmd += save_and_quit % (pdb_filename,
                          pdb_filename,
    )
  f=file(tleap_input, "wb")
  f.write(cmd)
  f.close()

def run_tleap_cmd_file(tleap_input):
  fatal_errors = [
    "Failed to generate parameters",
    "Unknown residue",
    ]
  non_fatal_errors = [
    "Could not open file",
    "duplicate [",
    ]
  cmd = "tleap -f %s" % tleap_input
  ero = easy_run.fully_buffered(command=cmd)
  err = StringIO.StringIO()
  ero.show_stdout(out=err)
  outl = ""
  for line in err.getvalue().split("\n"):
    for fe in non_fatal_errors:
      if line.find(fe)>-1:
        print line
    for fe in fatal_errors:
      if line.find(fe)>-1:
        raise Sorry(line)

def run_tleap(pdb_filename):
  tleap_input = "%s.in" % pdb_filename
  write_tleap_cmd_file(pdb_filename, tleap_input)
  run_tleap_cmd_file(tleap_input)

def adjust_his(atom_group):
  hydrogens = ["ND1", "NE2"]
  found = []
  for atom in atom_group.atoms():
    if atom.name.strip() in hydrogens:
      found.append(atom.name.strip())
  if len(found)==2:
    new_residue = "HIP"
  elif len(found)==1:
    if "ND1" in found:
      new_residue = "HID"
    elif "NE2" in found:
      new_residue = "HIE"
    else: assert 0
  else: assert 0
  atom_group.resname = new_residue

def convert_to_amber(hierarchy):
  # assert hydrogens
  for atom_group in hierarchy.atom_groups():
    if atom_group.resname in changing_residues:
      if atom_group.resname in his_names:
        adjust_his(atom_group)

def generate_altloc_filenames(hierarchy):
  conformer_hierarchy = hierarchy.deep_copy()
  for model in conformer_hierarchy.models():
    for chain in model.chains():
      conformers = chain.conformers()
      break
  for i in range(len(conformers)):
    altloc_hierarchy = hierarchy.deep_copy()
    for residue_group in altloc_hierarchy.residue_groups():
      blank = 0
      remove = []
      for j, atom_group in enumerate(residue_group.atom_groups()):
        if j==0 and not atom_group.altloc.strip():
          blank=1
          continue
        if i+blank!=j:
          remove.append(atom_group)
      if remove:
        for r in remove:
          residue_group.remove_atom_group(r)
    output_filename = "tleap_input_%d.pdb" % (i+1)
    f=file(output_filename, "wb")
    f.write(altloc_hierarchy.as_pdb_string(
      #crystal_symmetry=pdb_inp.crystal_symmetry()
      )
      )
    f.close()
    yield output_filename

def convert_tleap_to_phenix(filename):
  pdb_inp = pdb.input(filename)
  hierarchy = pdb_inp.construct_hierarchy()
  for atom_group in hierarchy.atom_groups():
    if atom_group.resname in his_names:
      atom_group.resname = "HIS"
    elif atom_group.resname == "CYX":
      atom_group.resname = "CYS"
    elif atom_group.resname in ["WAT"]:
      atom_group.resname = "HOH"
  return hierarchy

def adjust_occupany(hierarchy):
  for atom in hierarchy.atoms():
    atom.occ = 1
  return hierarchy

def adjust_b_factors(hierarchy, tleap_hierarchy):
  return tleap_hierarchy

def run(rargs):
  print """
  Convert a standard model to AMBER input using "tleap"
  """
  rargs = list(rargs)
  parser = setup_parser()
  (options, args) = parser.parse_args(args=rargs)
  if options.show_defaults:
    for line in master_phil_string.splitlines():
      if line.strip().find(".")==0: continue
      print line
    sys.exit()
  if len(args)==0:
    parser.print_help()
    sys.exit()  
  #
  argument_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_phil,
    home_scope="amber_tleap")
  #
  phils = []
  phil_args = []
  pdbs = []
  for arg in args:
    if os.path.isfile(arg):
      if iotbx.pdb.is_pdb_file(arg):
        pdbs.append(arg)
        continue
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
  #working_phil.show()
  working_params = working_phil.extract()

  assert len(pdbs)==1, "only one model"

  params = working_params.amber_tleap
  params.input.pdb_file_name = pdbs[0]
  if params.output.pdb_file_name is None:
    preamble = params.input.pdb_file_name.split(".")[0]
    params.output.pdb_file_name = "%s_amber.pdb" % preamble

  working_phil.format(python_object=working_params).show()

  pdb_inp = pdb.input(params.input.pdb_file_name)
  hierarchy = pdb_inp.construct_hierarchy()
  convert_to_amber(hierarchy)
  for i, tleap_input_filename in enumerate(generate_altloc_filenames(hierarchy)):
    print tleap_input_filename
    # run tleap
    tleap_output_filename = "%s_tleap_output.pdb" % preamble
    if 1:
      pass #os.system("cp 1exr_tleap_output.pdb %s" % tleap_output_filename)
    else:
      run_tleap(tleap_input_filename,
                tleap_output_filename,
        )
    # convert tleap output model to phenix model
    tleap_hierarchy = convert_tleap_to_phenix(tleap_output_filename)
    # adjust b-factors from input
    tleap_hierarchy = adjust_occupany(tleap_hierarchy)
    tleap_hierarchy = adjust_b_factors(hierarchy, tleap_hierarchy)

    break

  f=file(params.output.pdb_file_name, "wb")
  f.write(tleap_hierarchy.as_pdb_string(
    crystal_symmetry=pdb_inp.crystal_symmetry()),
          )
  f.close()
  print "\n  Output model written to",params.output.pdb_file_name


if __name__=="__main__":
  run(sys.argv[1:])
  
