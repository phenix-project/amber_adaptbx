from __future__ import print_function
import os
import subprocess
import parmed
from amber_adaptbx.tests.utils import get_fn, tempfolder, equal_files
from glob import glob
import pytest

template = '''
phenix.AmberPrep {test_folder}/{pdb_code}.pdb minimise=amber_all

{test_folder}/dacdif {test_folder}/non-les/4phenix_{pdb_code}.pdb 4phenix_{pdb_code}.pdb
{test_folder}/dacdif {test_folder}/non-les/4amber_{pdb_code}.prmtop 4amber_{pdb_code}.prmtop
{test_folder}/dacdif {test_folder}/non-les/4amber_{pdb_code}.rst7 4amber_{pdb_code}.rst7
{test_folder}/dacdif {test_folder}/non-les/{pdb_code}.min.out {pdb_code}.min.out

echo "Running LES test on {pdb_code} with a gap:"
phenix.AmberPrep {test_folder}/{pdb_code}.pdb minimise=amber_all LES=true

{test_folder}/dacdif {test_folder}/les/4phenix_{pdb_code}.pdb 4phenix_{pdb_code}.pdb
{test_folder}/dacdif {test_folder}/les/4amber_{pdb_code}.prmtop 4amber_{pdb_code}.prmtop
{test_folder}/dacdif {test_folder}/les/4amber_{pdb_code}.rst7 4amber_{pdb_code}.rst7
{test_folder}/dacdif {test_folder}/les/{pdb_code}.min.out {pdb_code}.min.out
'''


# FIXME: turn on 3cfb test?
# @pytest.mark.parametrize('pdb_code', ['3cfb', '1aho'])
@pytest.mark.parametrize('pdb_code', ['1aho'])
def test_gap(pdb_code):
    test_AmberPrep_fn = get_fn(
        '{pdb_code}/test_AmberPrep'.format(pdb_code=pdb_code))
    test_folder = os.path.dirname(test_AmberPrep_fn)

    outputs = []
    lines = [
        line
        for line in template.format(
            test_folder=test_folder, pdb_code=pdb_code).split('\n')
        if line.startswith('phenix') or 'dacdif' in line
    ]

    for line in lines:
        print(('command = ', line))
        outputs.append(subprocess.check_output(line, shell=True).decode())
        print(outputs[-1])
    assert 'possible FAILURE:' not in ''.join(outputs)
