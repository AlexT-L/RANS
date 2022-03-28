"""Used to set up an environment to run validation tests
    Libraries/Modules:
        pytest\n
        numpy\n
        Field\n
        Input\n
        AirfoilMap\n
        CellCenterWS\n
        NS_Airfoil\n
        NavierStokes\n
    
        """

import sys
sys.path.append("../../RANS/bin")
#sys.path.append("../../../RANS/bin")

import pytest
import numpy as np
from bin.Field import Field, isfinite, max, min, mean, abs
from bin.Input import Input
from bin.AirfoilMap import AirfoilMap
from bin.CellCenterWS import CellCenterWS
from bin.NS_Airfoil import NS_Airfoil
from bin.NavierStokes import NavierStokes


def get_environment(filename=None):
    # create input and grid
    if filename is None:
        filename = 'rae9-s1.data'

    # read in input
    input = Input(filename)

    # format input
    input.geo_param["inflation_layer"] = (input.flo_param["kvis"] != 0)
    gridInput = input.add_dicts(input.geo_param, input.in_var)
    grid_dim = [input.dims['nx'], input.dims['ny']]
    modelInput = input.add_dicts(input.flo_param, input.solv_param)

    # create geometry objects
    grid = AirfoilMap.from_file(grid_dim, gridInput)
    ws = CellCenterWS(grid)

    # create physics objects
    bcmodel = NS_Airfoil(modelInput)
    model = NavierStokes(bcmodel, modelInput)

    # initialize state
    [nx, ny] = [grid.dims['nx'], grid.dims['ny']]
    dim = model.dim()
    state = np.zeros((nx, ny, dim))
    model.init_state(ws,state)
    
    return [model, ws, state, [nx, ny]]