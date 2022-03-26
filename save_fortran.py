"""This module creates and saves numpy arrays to files to be used in validation tests

    Libraries/Modules:
        numpy\n
        typing\n
        Field\n
        Grid\n
        Input\n
        airfoil_map\n
    
        """

import sys
#sys.path.append("../../RANS/bin")
#sys.path.append("../../../RANS/bin")

import numpy as np
from bin.Field import Field, isfinite
from bin.Input import Input
from bin.AirfoilMap import AirfoilMap
from bin.CellCenterWS import CellCenterWS
from bin.NS_Airfoil import NS_Airfoil
from bin.NavierStokes import NavierStokes

# create input and grid
filename = 'rae9-s1.data'

# read in input
input = Input(filename) # Will actually take all command line inputs
print("Input Loaded")

# format input
input.geo_param["inflation_layer"] = (input.flo_param["kvis"] != 0)
gridInput = input.add_dicts(input.geo_param, input.in_var)
grid_dim = [input.dims['nx'], input.dims['ny']]
modelInput = input.add_dicts(input.flo_param, input.solv_param)

# create geometry objects
grid = AirfoilMap.from_file(grid_dim, gridInput)
print("Grid Created")
ws = CellCenterWS(grid)
print("Workspace Created")

# create physics objects
bcmodel = NS_Airfoil(modelInput)
model = NavierStokes(bcmodel, modelInput)
print("Model Created")

# create faux fields
[nx, ny] = [grid.dims['nx'], grid.dims['ny']]
dim = model.dim()
state = np.zeros((nx, ny, dim))
dw = np.zeros((nx, ny, dim))
model.init_state(ws,state)
model.test_eflux(ws,state,dw)
np.save('bin/validation/eflux.npy', dw)