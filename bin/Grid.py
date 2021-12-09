
# oversimplified version of Grid class just for now
# so I can write other stuff

import numpy as np
from Field import Field

class Grid:
    
    def __init__(self, input):

        # let Grid contain the variables in dims.f
        dims = input.dims        
        nx = dims['nx']
        ny = dims['ny']
        
        # defining things to be consistent with gmesh.f
        self.nx = nx
        self.ny = ny
        self.il = nx + 1 # this is computational domain dimensions?
        self.jl = ny + 1
        self.ie = nx + 2
        self.je = ny + 2
        self.ib = nx + 3
        self.jb = ny + 3
        self.itu = np.NaN # unsure what these guys are yet
        self.itl = np.NaN
        
        # it will also contain the variables in mesh_var as members
        # these fields will require a lot more math but I barely 
        # even know how to spell conformal mapping so I leave that to someone else

        # rectangular cartesian grid for now
        x_vec = np.linspace(dims['x_bound'][0], dims['x_bound'][1], self.ib)
        y_vec = np.linspace(dims['y_bound'][0], dims['y_bound'][1], self.jb)
        xg, yg = np.meshgrid(x_vec, y_vec)
        dims = np.array([self.ib, self.jb, 2])
        
        # physical vertex locations
        self.X = Field(np.zeros(dims)) 
        self.X.vals[:,:,0] = xg.T
        self.X.vals[:,:,1] = yg.T

        # cell volumes, porosity, and far field mask
        self.Vol = Field( np.ones([self.ib, self.jb]) )
        self.PorJ = Field( np.ones([self.ib, self.jb]) )
        self.PorI = Field( np.ones([self.ib, self.jb]) )
        self.Fint = Field( np.ones([self.ib, self.jb]) )
