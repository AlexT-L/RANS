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


def test_ynot_validation():    
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
    ynot = np.zeros(nx)
    model.init_state(ws,state)
    model.test(ws,state,ynot,'ynot')

    # compare with fortran
    if UPDATE_FORTRAN_DATA:
        ynot_fortran = np.zeros(ynot.shape)
        model.test(ws,state,ynot_fortran,'ynot','fortran')
        np.save('tests/validation/ynot.npy', ynot_fortran)
    TOL = 1e-5
    ynot_fortan = np.load('tests/validation/ynot.npy', allow_pickle=False)

    print ("max(ynot_fortran) = "+str(max(ynot_fortan)))
    print ("max(ynot) = "+str(max(ynot)))
    print ("mean(ynot_fortran-dw) = "+str(mean(ynot_fortan-ynot)))
    print ("min(ynot_fortran-dw) = "+str(min(abs(ynot_fortan-ynot))))
    assert max(abs(ynot_fortan - ynot)) < TOL


def test_dsti_validation():    
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
    dsti = np.zeros(nx)
    model.init_state(ws,state)
    model.test(ws,state,dsti,'dsti')

    # compare with fortran
    if UPDATE_FORTRAN_DATA:
        dsti_fortran = np.zeros(dsti.shape)
        model.test(ws,state,dsti_fortran,'dsti','fortran')
        np.save('tests/validation/dsti.npy', dsti_fortran)
    TOL = 1e-8
    dsti_fortan = np.load('tests/validation/dsti.npy', allow_pickle=False)

    print ("max(dsti_fortran) = "+str(max(abs(dsti_fortan))))
    print ("max(dsti) = "+str(max(abs(dsti))))
    print ("mean(dsti_fortran-dw) = "+str(mean(dsti_fortan-dsti)))
    print ("min(dsti_fortran-dw) = "+str(min(abs(dsti_fortan-dsti))))
    assert max(abs(dsti_fortan - dsti)) < TOL
    