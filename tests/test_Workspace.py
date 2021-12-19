"""
Description
-----------
Tests the Workspace object to see if it works as expected

Libraries/Modules
-----------------
-pytest \n
-Workspace

Notes
-----
Runs the following tests:\n
Verify x, xc, vol can be retrieved
Verify x, xc, vol are non-zero
Verify that init_vars works as expected
then check that init_vars? exists return True
check that has_dict? exists and has_dict return false when they should
check is_finest()
"""
import sys
# sys.path.append('../../)
sys.path.append('../../RANS/bin')
from numpy.core.numeric import NaN
import pytest
from Input import Input
from Field import Field
from AirfoilMap import AirfoilMap
from Workspace import Workspace
from CellCenterWS import CellCenterWS
import numpy as np


def test_x():
    """
    Asserts that we can retrieve x and that it is non zero
    """
    # input = Input('rae9-s1.data') # Will actually take all command line inputs
    filename = 'rae9-s1.data'
    # read in input
    input = Input(filename) # Will actually take all command line inputs
    # format input
    input.geo_param["inflation_layer"] = (input.flo_param["kvis"] != 0)
    gridInput = input.add_dicts(input.geo_param, input.in_var)
    grid_dim = [input.dims['nx'], input.dims['ny']]
    # create geometry objects
    grid = AirfoilMap.from_file(grid_dim, gridInput)
    workspace = CellCenterWS(grid)
    x = workspace.get_field('x')
    # Asserts non zero, as np.any is true when any are non-zero
    assert np.any(x)

def test_xc():
    """
    Asserts that we can retrieve xc and that it is non zero
    """
    # input = Input('rae9-s1.data') # Will actually take all command line inputs
    filename = 'rae9-s1.data'
    # read in input
    input = Input(filename) # Will actually take all command line inputs
    # format input
    input.geo_param["inflation_layer"] = (input.flo_param["kvis"] != 0)
    gridInput = input.add_dicts(input.geo_param, input.in_var)
    grid_dim = [input.dims['nx'], input.dims['ny']]
    # create geometry objects
    grid = AirfoilMap.from_file(grid_dim, gridInput)
    workspace = CellCenterWS(grid)
    xc = workspace.get_field('xc')
    # Asserts non zero, as np.any is true when any are non-zero
    assert np.any(xc)

def test_vol():
    """
    Asserts that we can retrieve vol and that it is non zero
    """
    # input = Input('rae9-s1.data') # Will actually take all command line inputs
    filename = 'rae9-s1.data'
    # read in input
    input = Input(filename) # Will actually take all command line inputs
    # format input
    input.geo_param["inflation_layer"] = (input.flo_param["kvis"] != 0)
    gridInput = input.add_dicts(input.geo_param, input.in_var)
    grid_dim = [input.dims['nx'], input.dims['ny']]
    # create geometry objects
    grid = AirfoilMap.from_file(grid_dim, gridInput)
    workspace = CellCenterWS(grid)
    vol = workspace.get_field('vol')
    # Asserts non zero, as np.any is true when any are non-zero
    assert np.any(vol)

# def test_init_vars():
#     '''
#     Verify that init_vars works as expected
#     then check that init_vars? exists return True
#     '''
#     # not sure about this yet
#     a=True
#     assert a

# def test_finest():
#     finest_check= Workspace.isFinest
#     '''
#     Asserts that isFinest 
#     '''
#     # but not always? cross check with mode?
#     assert finest_check


test_x()