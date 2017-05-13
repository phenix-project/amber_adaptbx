import os
import subprocess
import pytest
import numpy as np

from numpy.testing import assert_almost_equal as aa_eq
from amber_adaptbx.tests.utils import tempfolder, get_fn
import parmed as pmd
from parmed.amber.netcdffiles import NetCDFTraj


@pytest.mark.xfail
def test_netcdf_trajectory():
    command_refine = [
        'phenix.refine', get_fn('vAla3/vAla3.pdb'), get_fn('vAla3/vAla3.cif'),
        get_fn('vAla3/vAla3.mtz'),
        'topology_file_name={}'.format(get_fn('vAla3/vAla3.prmtop')),
        'coordinate_file_name={}'.format(get_fn('vAla3/vAla3.rst7')),
        'use_amber=True', 'wxc_scale=0.025', '--overwrite',
        'refinement.main.number_of_macro_cycles=1',
        'amber.netcdf_trajectory_file_name=hello.nc'
    ]

    n_frames = 59
    parm = pmd.load_file(get_fn('vAla3/vAla3.prmtop'))
    expected_boxes = np.array(
        n_frames * [
            [30., 30., 30., 90., 90., 90.],
        ], dtype='f8')

    with tempfolder():
        output = subprocess.check_output(command_refine)
        traj = NetCDFTraj.open_old('hello.nc')
        print(traj.coordinates.shape)
        aa_eq(traj.box, expected_boxes)
        assert traj.coordinates.shape == (n_frames, len(parm.atoms), 3)
