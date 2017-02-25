import os
import subprocess
import parmed
from amber_adaptbx.tests.utils import get_fn, tempfolder, equal_files
from glob import glob


def test_tleap_input():
    code = '4lzt'
    pdb_fn = get_fn('{code}/{code}.pdb'.format(code=code))

    regression_folder = os.path.dirname(pdb_fn) + '/regression'

    excluded_files = ['README.md']

    all_expected_files = [os.path.basename(fn) for fn in glob(regression_folder + '/*')]
    for fn in excluded_files:
        all_expected_files.remove(fn)

    print(all_expected_files)

    with tempfolder():
        command_non_les = [
            'phenix.AmberPrep',
            pdb_fn,
            'clean=False',
        ]
        subprocess.check_call(command_non_les)
        for fn in all_expected_files:
            expected_fn = os.path.join(regression_folder, fn)
            equal_files(expected_fn, fn)

def test_gap_in_tleap():
    code = '4lzt'
    pdb_fn = get_fn('{code}/{code}.pdb'.format(code=code))
    parm = parmed.load_file(pdb_fn)

    # with tempfolder():
    if True:
        # create a new protein with only 4 residues
        # there are 2 gaps 
        new_parm = parm[':1,3,9,10']
        gap_pdb_fn = 'new0.pdb'
        new_parm.save(gap_pdb_fn, overwrite=True)
        subprocess.check_call('cat {} | sed "/TER/d" > new1.pdb'.format(gap_pdb_fn), shell=True)
        for residue in new_parm.residues:
            residue.ter = False
        command_non_les = [
            'phenix.AmberPrep',
            'new1.pdb',
            'clean=False',
        ]
        subprocess.check_call(command_non_les)
