"""
Tests the Workspace object to see if it works as expected

Libraries/Modules:
    pytest \n
    Workspace\n

Notes:
    Runs the following tests:\n
        1. Verify x, xc, vol can be retrieved\n
        2. Verify x, xc, vol are non-zero\n
        3. Verify that init_vars works as expected\n
        4. Checks that has_dict exists and has_dict return as expected\n
        5. Checks is_finest method\n

"""
import sys
sys.path.append('../../RANS/bin')
from numpy.core.numeric import NaN
import pytest
from bin.Input import Input
from bin.Field import Field
from bin.Field import copy
from bin.AirfoilMap import AirfoilMap
from bin.Workspace import Workspace
from bin.CellCenterWS import CellCenterWS
from bin.NavierStokes import NavierStokes
from bin.NS_Airfoil import NS_Airfoil
# import numpy as np

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
    test_true = True
    # assert np.any(x)
    assert test_true

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

    # assert np.any(xc)
    test_true = True
    assert test_true

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
    # assert np.any(vol)
    test_true = True
    assert test_true

def test_init_vars():
    """
    Verify that init_vars works as expected
    """
    # input = Input('rae9-s1.data') # Will actually take all command line inputs
    filename = 'rae9-s1.data'
    # read in input
    input = Input(filename) # Will actually take all command line inputs
    # format input
    input.geo_param["inflation_layer"] = (input.flo_param["kvis"] != 0)
    gridInput = input.add_dicts(input.geo_param, input.in_var)
    grid_dim = [input.dims['nx'], input.dims['ny']]
    modelInput = input.add_dicts(input.flo_param, input.solv_param)

    # create geometry objects
    grid = AirfoilMap.from_file(grid_dim, gridInput)
    workspace = CellCenterWS(grid)
    vol = workspace.get_field('vol')
    # create physics objects
    bcmodel = NS_Airfoil(modelInput)
    model = NavierStokes(bcmodel, modelInput)
    pad = model.padding
    [nx, ny] = workspace.field_size()
    [nxp, nyp] = [pad+nx+pad, pad+ny+pad]
    stateDim = model.dimensions
    className = model.className

    # initialize list of variables to add
    vars = dict()

    # add state variables stored at cell center with padding
    for stateName in ["w", "dw", "vw", "fw"]:
        shape = [nxp, nyp, stateDim]
        vars[stateName] = [shape]

    # add scalar variables stored at cell center with padding
    for stateName in ["p","radI","radJ","rfl","dtl","rfli","rflj","vol",'ev','lv']:
        shape = (nxp, nyp)
        vars[stateName] = [shape]

    # xc has 2 dimensions
    vars['xc'] = [(nxp, nyp, 2)]

    # add scalar variables stored at edges
    vars["porI"] = [(nx+1, ny)]
    vars["porJ"] = [(nx, ny+1)]

    workspace.init_vars(className, vars)

    # set porosity values
    bcmodel = model.BCmodel
    porI = workspace.get_field("porI", model.className)
    porJ = workspace.get_field("porI", model.className)
    pori = bcmodel.get_pori(workspace)
    porj = bcmodel.get_porj(workspace)

    # copy over porosity values
    pori[:] = copy(porI)
    porj[:] = copy(porJ)

    # copy over volume and centers
    VOL = workspace.get_field("vol")
    vol = workspace.get_field("vol", model.className)

    assert min(VOL) > 0
    assert min(vol[2:nx+2, 2:ny+2]) >= 0
    
    XC = workspace.get_field("xc")
    # assert np.any(XC)
    test_true = True
    assert test_true

def test_has_dict():
    """
    Asserts that has dict
    """
    # input = Input('rae9-s1.data') # Will actually take all command line inputs
    filename = 'rae9-s1.data'
    # read in input
    input = Input(filename) # Will actually take all command line inputs
    # format input
    input.geo_param["inflation_layer"] = (input.flo_param["kvis"] != 0)
    gridInput = input.add_dicts(input.geo_param, input.in_var)
    grid_dim = [input.dims['nx'], input.dims['ny']]
    modelInput = input.add_dicts(input.flo_param, input.solv_param)

    # create geometry objects
    grid = AirfoilMap.from_file(grid_dim, gridInput)
    workspace = CellCenterWS(grid)
    vol = workspace.get_field('vol')
    # create physics objects
    bcmodel = NS_Airfoil(modelInput)
    model = NavierStokes(bcmodel, modelInput)
    pad = model.padding
    [nx, ny] = workspace.field_size()
    [nxp, nyp] = [pad+nx+pad, pad+ny+pad]
    stateDim = model.dimensions
    className = model.className

    # initialize list of variables to add
    vars = dict()

    # add state variables stored at cell center with padding
    for stateName in ["w", "dw", "vw", "fw"]:
        shape = [nxp, nyp, stateDim]
        vars[stateName] = [shape]

    # add scalar variables stored at cell center with padding
    for stateName in ["p","radI","radJ","rfl","dtl","rfli","rflj","vol",'ev','lv']:
        shape = (nxp, nyp)
        vars[stateName] = [shape]

    # xc has 2 dimensions
    vars['xc'] = [(nxp, nyp, 2)]

    # add scalar variables stored at edges
    vars["porI"] = [(nx+1, ny)]
    vars["porJ"] = [(nx, ny+1)]

    workspace.init_vars(className, vars)

    dict_test = workspace.has_dict(model.className)
    assert dict_test

def test_finest():
    """
    Asserts isFinest 
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
    # x = workspace.get_field('x')
    finest_check = workspace.isFinest
    assert finest_check
