# libraries
import numpy as np

# append to path so we can access files
import sys
sys.path.append("../../../")

# project specific dependencies
from Field import Field
from Contractinator import conservative4way, simple, sum4way

# grid creation functions
from utils.coord_strch_func import coord_stretch
from utils.geom_func import geom 
from utils.mesh_func import mesh
from utils.sangho_func import sangho
from utils.plot_mesh_func import plot_mesh

def init_from_file(self, grid_dim, input):
    # read in number of divisions
    [nx, ny] = grid_dim
    nx = int(nx)
    ny = int(ny)
    self.divisions = (nx, ny)
    
    #append to geo_param dictionary
    self.geo=input
    geo =self.geo
    geo["scal"] = 0
    geo["trail"]=geo["trail"]*np.pi/180#convert trail angle to radians
    geo["nbl"] = 0
    geo["jlinv"] = 0
    geo["nyinv"] = 0

    xte=self.geo['xte']    

    # set mesh dimensions and update dims
    dim         = dict()
    self.dims   = dim     
    dim["nx"]   = nx
    dim["ny"]   = ny
    dim["il"]   = nx + 1 # number of points/edges in i dir of computational grid
    dim["jl"]   = ny + 1 # number of points/edges in j dir of computational grid
        # values below define number of points in
        # computational domain for a padded grid
    dim["ie"]   = nx + 2
    dim["je"]   = ny + 2
    dim["ib"]   = nx + 3
    dim["jb"]   = ny + 3
    self.il     = dim["il"]
    self.jl     = dim["jl"]
    il          = dim["il"]
    jl          = dim["jl"]
    ib          = dim["ib"]
    jb          = dim["jb"]
    

    # initialize x-y vertex, center, and volume fields
    self.x  = Field((il,jl),2)
    self.xc = Field((ib+1,jb+1),2)
    self.vol= Field((ib+1,jb+1),1)

    # store relevant fields
    fields = dict()
    fields['x'] = self.x
    fields['xc'] = self.xc
    fields['vol'] = self.vol
    self.fields = fields


    #set the limits of the aerfoil profile
    self.ite       = int(0.5*xte*nx) #coordinate of trailing edge in physical space
    self.ile       = np.floor(il/2  +1) #coordinate of leading edge in physical space
    self.itl       = int(self.ile - self.ite) #lower coordinate of trailing edge in computational space
    self.itu       = int(self.ile + self.ite) #upper coordinate of trailing edge in computational space
    
    self.a  = np.array([np.zeros(il),np.zeros(il)])
    self.b0 = np.zeros(jl)
    self.s0 = np.zeros(il)

    #define point distributions in each coordinate direction (coordinate stretching)
    coord_stretch(self)

    #sqrt root mapping of aerfoil profile to slit & interpolating 
    # to make sure points on aerfoil match those of the grid
    geom(self)
    
    #making mesh 
    mesh(self)

    #sanghos modification for non-dimensionalization(re-scaling the c-mesh)
    sangho(self)

    #plot mesh
    plot_mesh(self, self.x)
    

def init_from_grid(newGrid, grid):
    # get relavant values from old grid
    fields = grid.fields
    x = fields['x']
    xc = fields['xc']
    vol = fields['vol']
    [nx, ny] = grid.divisions
    [nxNew, nyNew] = newGrid.divisions

    # create new arrays
    xNew = Field(newGrid.get_size())
    xcNew = Field((nxNew + 4, nyNew + 4))
    volNew = Field((nxNew + 4, nyNew + 4))

    # condense mesh
    simple(x, xNew)
    conservative4way(xc, xcNew)
    sum4way(vol[2:nx+2, 2:ny+2], volNew[2:nxNew+2, 2:nyNew+2])

    # store fields
    newFields = dict()
    newFields['x'] = x
    newFields['xc'] = xc
    newFields['vol'] = vol
    newGrid.fields = newFields