import numpy as np
from Field import Field

class Grid:
    
    def __init__(self, input):

        # let Grid contain the variables in input.dims and input.geo_param
        self.dims = input.dims        
        nx = self.dims['nx']
        ny = self.dims['ny']
        self.geo=input.geo_param
        xte=self.geo['xte']
        self.flo=input.flo_param
        kvis=self.flo["kvis"]

        self.nx = nx 
        self.ny = ny 
        # set mesh dimensions
        self.il = nx + 1 # number of points/edges in i dir of computational grid
        self.jl = ny + 1 # number of points/edges in j dir of computational grid
         # values below define number of points in
         # computational domain for a padded grid
        self.ie = nx + 2 
        self.je = ny + 2 
        self.ib = nx + 3
        self.jb = ny + 3

        #set the limits of the aerfoil profile
        self.ite       = .5000005*xte*nx #coordinate of trailing edge in computationa domain
        self.ile       = self.il/2  +1 #coordinate of leading edge in computationa domain
        self.itl       = self.ile - self.ite #not sure what these do
        self.itu       = self.ile + self.ite
        
        #set the limits of the outer mesh for a viscous simulation (kvis>1)
        if kvis > 1:
            nbl       = self.jl/2
            self.ny        = ny  -nbl
            self.jl        = self.jl  -nbl
        
        #define point distributions in each coordinate direction






        # it will also contain the variables in mesh_var as members
        # these fields will require a lot more math but I barely 
        # even know how to spell conformal mapping so I leave that to someone else

        # rectangular cartesian grid for now
        x_vec = np.linspace(dims['x_bound'][0], dims['x_bound'][1], self.ib)
        y_vec = np.linspace(dims['y_bound'][0], dims['y_bound'][1], self.jb)
        xg, yg = np.meshgrid(x_vec, y_vec)
        
        # physical vertex locations
        self.X = Field(self.get_size(), stateDim=2) 
        self.X.vals[:,:,0] = xg.T
        self.X.vals[:,:,1] = yg.T

        # cell volumes, porosity, and far field mask
        self.Vol = Field(self.get_size())
        self.PorJ = Field(self.get_size())
        self.PorI = Field(self.get_size())
        self.Fint = Field(self.get_size())

    def get_size(self):
        return [self.ib, self.jb]


