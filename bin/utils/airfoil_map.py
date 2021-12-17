# libraries
import numpy as np

# append to path so we can access files
import sys

from bin.Field import isfinite
sys.path.append("../../../")

# project specific dependencies
from Field import Field
from Contractinator import conservative4way, simple, sum4way

# grid creation functions
from utils.dims_func import set_dims
from utils.coord_strch_func import coord_stretch
from utils.geom_func import geom 
from utils.mesh_func import mesh
from utils.sangho_func import sangho
from utils.metric_func import metric
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

    # set dimensions
    set_dims(self)

    # initialize x-y vertex, center, and volume fields
    il = nx+1
    jl = ny+1
    self.x  = Field((il,jl,2))
    self.xc = Field((nx,ny,2))
    self.vol= Field((nx,ny,1))

    # store relevant fields
    fields = dict()
    fields['x'] = self.x
    fields['xc'] = self.xc
    fields['vol'] = self.vol
    self.fields = fields

    # create arrays for mesh generation
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

    # calculate cell volumes and centers
    metric(self)

    #plot mesh
    plot_mesh(self)
    

def init_from_grid(newGrid, grid):
    # transfer geometry info and set new dimensions
    newGrid.geo = grid.geo
    set_dims(newGrid)

    # get relavant values from old grid
    fields = grid.fields
    x = fields['x']
    xc = fields['xc']
    vol = fields['vol']

    # create new arrays
    [nx, ny] = newGrid.get_size()
    xNew = Field((nx,ny,2))
    xcNew = Field((nx,ny,2))
    volNew = Field((nx,ny,1))

    # condense mesh
    simple(x, xNew)
    conservative4way(xc, xcNew)
    sum4way(vol, volNew)

    # store fields
    newFields = dict()
    newFields['x'] = xNew
    newFields['xc'] = xcNew
    newFields['vol'] = volNew
    newGrid.fields = newFields

    #plot mesh
    plot_mesh(newGrid)