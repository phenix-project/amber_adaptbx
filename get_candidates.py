import os, sys
import time
import copy
import pickle
import gzip
from libtbx import easy_run
from iotbx import pdb
import StringIO

from cctbx.array_family import flex
from mmtbx import model_vs_data
from libtbx import easy_pickle
import libtbx.load_env
    
def generate_pdb_codes_01():
  done = []
  #for item in generate_pdb_codes_01_additional_codes():
  #  done.append(item)
  #  yield item
  for item in generate_pdb_codes_01_pavel():
    if item in done: continue
    done.append(item)
    yield item
  for item in generate_pdb_codes_01_high_res():
    if item in done:
      continue
    done.append(item)
    yield item

def generate_pdb_codes_01_additional_codes():
  for apc in additional_pdb_codes:
    yield apc

def generate_pdb_codes_01_high_res():
  for apc in generate_pdb_codes_01_pavel(
      d_min_max=1.05,
      diff_max=0.005,
      ):
    yield apc

def generate_pdb_codes_01_pavel(d_min_max      = 3.549,
                                r_work_pdb_max = 0.30,
                                r_free_pdb_max = 0.35,
                                diff_max       = 0.015,
                                ):
  file = libtbx.env.find_in_repositories(relative_path=
    "chem_data/polygon_data/all_mvd.pickle", test=os.path.isfile)
  database_dict = easy_pickle.load(file)
  #
  pdb_code           = database_dict["pdb_code"]
  r_work_pdb         = database_dict["pdb_header_r_work"]
  r_free_pdb         = database_dict["pdb_header_r_free"]
  r_work_re_computed = database_dict["r_work"]
  r_free_re_computed = database_dict["r_free"]
  d_min              = database_dict["high_resolution"]
  twinned            = database_dict["twinned"]
  cmpl_in_range      = database_dict["completeness_in_range"]
  resname_classes    = database_dict["resname_classes"]
  #print '1'*80
  #for i,(a,b) in enumerate(zip(pdb_code, d_min)):
  #  print 'pdb_code',a,'d_min',b
  #  if i>10: break
  #
  sel = pdb_code != "none"
  sel &= r_work_pdb != "none"
  sel &= r_free_pdb != "none"
  sel &= r_work_re_computed != "none"
  sel &= r_free_re_computed != "none"
  sel &= d_min != "none"
  sel &= twinned != "none"
  sel &= cmpl_in_range != "none"
  sel &= resname_classes != "none"
  #
  pdb_code           = pdb_code.select(sel)
  r_work_pdb         = r_work_pdb.select(sel)
  r_free_pdb         = r_free_pdb.select(sel)
  r_work_re_computed = r_work_re_computed.select(sel)
  r_free_re_computed = r_free_re_computed.select(sel)
  d_min              = d_min.select(sel) 
  twinned            = twinned.select(sel)
  cmpl_in_range      = cmpl_in_range.select(sel)
  resname_classes    = resname_classes.select(sel)
  #
  def str_to_float(x):
    tmp = flex.double()
    for x_ in x:
      tmp.append(float(x_))
    return tmp
  #
  d_min              = str_to_float(d_min)
  r_work_pdb         = str_to_float(r_work_pdb)
  r_free_pdb         = str_to_float(r_free_pdb)
  r_work_re_computed = str_to_float(r_work_re_computed)
  r_free_re_computed = str_to_float(r_free_re_computed)
  cmpl_in_range      = str_to_float(cmpl_in_range) 
  diff               = r_free_pdb - r_work_pdb
  #
  #
  sel  = d_min < d_min_max
  sel &= r_work_pdb < r_work_pdb_max
  sel &= r_free_pdb < r_free_pdb_max
  sel &= diff > diff_max #0.015 # this probably excludes 0.7 and better
  sel &= twinned == "false"
  sel &= cmpl_in_range > 0.9 
  #
  r_work_pdb_ = r_work_pdb.select(sel)
  r_free_pdb_ = r_free_pdb.select(sel)
  diff_       = diff.select(sel)
  name_       = pdb_code.select(sel)
  d_min_      = d_min.select(sel)
  resname_classes_ = resname_classes.select(sel)
  
  for n, rw, rf, dm, rc in zip(name_,r_work_pdb_,r_free_pdb_, d_min_, resname_classes_):
    if rc.find("amino_acid")==-1: continue
    #print "code:%s r_work:%f r_free:%f d_min:%f contents:%s" % (n, rw, rf, dm, rc)
    yield n
  
  print "Rwork:"
  show_histogram(data = r_work_pdb_, n_slots=5)
  print "Rfree:"
  show_histogram(data = r_free_pdb_, n_slots=5)
  print "Rfree - Rwork:"
  show_histogram(data = diff_, n_slots=5)
  print "d min:"
  show_histogram(data = d_min_, n_slots=20)

def generate_pdb_codes_amber(d_min_max      = 3.549,
                             r_work_pdb_max = 0.30,
                             r_free_pdb_max = 0.35,
                             diff_max       = 0.015,
                             exclude_resname_classes = [],
                             verbose        = False,
                             ):
  file = libtbx.env.find_in_repositories(relative_path=
    "chem_data/polygon_data/all_mvd.pickle", test=os.path.isfile)
  database_dict = easy_pickle.load(file)
  if 0:
    for key in sorted(database_dict.keys()):
      print key, list(database_dict[key][:9])
    assert 0
  count = 0
  for i, (sg, rc, na) in enumerate(zip(database_dict["space_group"],
                                   database_dict["resname_classes"],
                                   database_dict["number_of_atoms"],
                                   )):
    #if sg.find("p 1 (no. 1)")==-1: continue
    if int(na)>5000: continue
    for e in exclude_resname_classes:
      if rc.find(e)>-1: break
    else:
      if verbose:
        print "  %5d %s %s %s" % (i+1,database_dict["pdb_code"][i],sg,rc)
      yield database_dict["pdb_code"][i]
      count+=1
    continue

def run(only_i=None,
        only_code=None,
        ):
  list_only=False
  if only_i=="list":
    list_only=True
  try: only_i=int(only_i)
  except: only_i=None
  if only_i==0: only_i=None
  dry_run=False
  if only_i is not None and only_i<0:
    only_i=abs(only_i)
    dry_run=True

  codes = []
  for i, code in enumerate(generate_pdb_codes_amber(
      exclude_resname_classes=[
        "other",
        "rna_dna",
        ])):
    print '...',i, code
    codes.append(code)

    if only_i is not None: break
    if only_code is not None: break

  if 0:
    f=file("90_EDS_pairs.csv", "rb")
    lines = f.readlines()
    f.close()

    hi_lo = []
    for i, line in enumerate(lines):
      if not i: continue
      print line
      tmp = line.split(",")
      print tmp
      if tmp[0] in codes:
        hi_lo.append(tmp)
      if tmp[0] in codes and tmp[3] in codes:
        assert 0
    print hi_lo
    for h in hi_lo:
      print h


if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
