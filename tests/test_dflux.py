"""This module tests the python version of eflux against the fortran version

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
from tests.validation.save_fortran import save_fortan
from bin.Field import Field, isfinite, max, min, mean, abs
from bin.Input import Input
from bin.AirfoilMap import AirfoilMap
from bin.CellCenterWS import CellCenterWS
from bin.NS_Airfoil import NS_Airfoil
from bin.NavierStokes import NavierStokes

def test_dflux_validation():
    save_fortan()
    
    # create input and grid
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

    # create faux fields
    [nx, ny] = [grid.dims['nx'], grid.dims['ny']]
    dim = model.dim()
    state = np.zeros((nx, ny, dim))
    dw = np.zeros((nx, ny, dim))
    model.init_state(ws,state)
    model.test(ws,state,dw,'dflux')

    # compare with fortran
    TOL = 1e-5
    dw_fortan = np.load('tests/validation/dflux.npy', allow_pickle=False)

    print ("max(dw_fortran) = "+str(max(dw_fortan)))
    print ("max(dw) = "+str(max(dw)))
    print ("mean(dw_fortran-dw) = "+str(mean(dw_fortan-dw)))
    print ("min(dw_fortran-dw) = "+str(min(abs(dw_fortan-dw))))
    # assert max(abs(dw_fortan - dw)) < TOL
    