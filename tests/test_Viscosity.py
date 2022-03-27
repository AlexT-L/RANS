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
from bin.Field import Field, isfinite, max, min, mean, abs
from bin.Input import Input
from bin.AirfoilMap import AirfoilMap
from bin.CellCenterWS import CellCenterWS
from bin.NS_Airfoil import NS_Airfoil
from bin.NavierStokes import UPDATE_FORTRAN_DATA, NavierStokes


def test_ev_validation():    
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
    ev = np.zeros((nx,ny))
    model.init_state(ws,state)
    model.test(ws,state,ev,'ev')

    # compare with fortran
    if UPDATE_FORTRAN_DATA:
        ev_fortran = np.zeros(ev.shape)
        model.test(ws,state,ev_fortran,'ev','fortran')
        np.save('tests/validation/ev.npy', ev_fortran)
    TOL = 1e-8
    ev_fortan = np.load('tests/validation/ev.npy', allow_pickle=False)

    print ("max(ev_fortran) = "+str(max(abs(ev_fortan))))
    print ("max(ev) = "+str(max(abs(ev))))
    print ("mean(ev_fortran-ev) = "+str(mean(ev_fortan-ev)))
    print ("min(ev_fortran-ev) = "+str(min(abs(ev_fortan-ev))))
    assert max(abs(ev_fortan - ev)) < TOL
    