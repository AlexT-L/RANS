# """
# Description
# 
# Tests the flux calculations

# Libraries/Modules
# 
# -pytest \n
# -Field
# -Model
# -Workspace
# -

# Notes
# 
# Runs the following tests:\n
# 1. Checks that a 2D Field can be created  \n
# 2. Checks that a 3D field can be created \n
# 3. Tests that we can set a whole field \n
# 4. Tests that we can set an individual elements \n
# 5. Tests store_product function \n
# 6. Tests store_difference function \n
# 7. Tests store_product function \n
# 8. Tests store_quotient function \n

# """



# import numpy as np
# from bin.Field import Field
# from bin.Input import Input
# from bin.AirfoilMap import AirfoilMap
# from bin.CellCenterWS import CellCenterWS
# from bin.NS_Airfoil import NS_Airfoil
# from bin.NavierStokes import NavierStokes

# # create input and grid
# filename = 'rae9-s1.data'

# # read in input
# input = Input(filename) # Will actually take all command line inputs
# print("Input Loaded")

# # format input
# input.geo_param["inflation_layer"] = (input.flo_param["kvis"] != 0)
# gridInput = input.add_dicts(input.geo_param, input.in_var)
# grid_dim = [input.dims['nx'], input.dims['ny']]
# modelInput = input.add_dicts(input.flo_param, input.solv_param)

# # create geometry objects
# grid = AirfoilMap.from_file(grid_dim, gridInput)
# print("Grid Created")
# ws = CellCenterWS(grid)
# print("Workspace Created")

# # create physics objects
# bcmodel = NS_Airfoil(modelInput)
# model = NavierStokes(bcmodel, modelInput)
# print("Model Created")

# # create faux fields
# np.random.seed(100)
# w_init = np.random.standard_normal([grid.dims['nx'],grid.dims['ny'],model.dim()])
# w = Field(0, w_init)
# res = np.zeros([grid.dims['nx'],grid.dims['ny'],model.dim()])
# dw = Field(0, res)

# #print(dw)
# model.init_state(ws,w)
# model.get_flux(ws,w,dw)
# print(dw)
