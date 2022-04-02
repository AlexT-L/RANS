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

import pytest
import numpy as np
from bin.Field import Field, isfinite, load, max, min, mean, abs, save, load
from bin.Input import Input
from bin.AirfoilMap import AirfoilMap
from bin.CellCenterWS import CellCenterWS
from bin.NS_Airfoil import NS_Airfoil
from bin.NavierStokes import NavierStokes, UPDATE_FORTRAN_DATA

TOLERANCE = 1.0e-6
RANDOM_STATE = True
UPDATE_STATE = False
DELTA = 0.5
FORCE_FAIL = False

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
    
    # randomize state
    if RANDOM_STATE:
        if UPDATE_STATE:
            sigma = np.random.rand(nx,ny,dim)*DELTA*2.0 - DELTA
            state = state*(1.0+sigma)
            save('tests/validation/state0', state)
        else:
            state = load('tests/validation/state0')
            
    
    return [model, ws, state, [nx, ny]]


# perform test of python script
def test_python(model, ws, state, testName):
    return model.test(ws,state,testName)
    
# compare with fortran
def test_fortran(model, ws, state, testName):
    if UPDATE_FORTRAN_DATA:
        out_fortran = model.test(ws,state,testName,'fortran')
        save('tests/validation/'+testName, out_fortran)
    return load('tests/validation/'+testName)


def analyze(output, output_fortran, outName="output"):
    
    print ("max("+outName+"_fortran) = "+str(max(output_fortran)))
    print ("max("+outName+") = "+str(max(output)))
    print ("mean("+outName+"_fortran - "+outName+") = "+str(mean(output_fortran-output)))
    print ("min("+outName+"_fortran - "+outName+") = "+str(min(abs(output_fortran-output))))
    print ("max("+outName+"_fortran - "+outName+") = "+str(max(abs(output_fortran-output))))
    assert max(abs((output_fortran - output)/max(abs(output)))) < TOLERANCE


def run_test(testName):
    # set up environment
    [model, ws, state, size] = get_environment()
    
    # perform test of python script
    out = test_python(model,ws,state,testName)

    # compare with fortran
    out_fort = test_fortran(model,ws,state,testName)

    analyze(out, out_fort, testName)
    if FORCE_FAIL: assert False