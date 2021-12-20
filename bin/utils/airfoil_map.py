"""This module has functions that perform the conformal mapping and the coarser grid objects

    Libraries/Modules:
        numpy\n
        Field\n
        Contractinator\n
        dims_funs\n
        coord_strch_func\n
        geom_func\n
        sangho_func\n
        metric_func\n
        plot_mesh_func\n
        """
# libraries
import numpy as np
from bin.Field import isfinite
from bin.Field import Field
from bin.Contractinator import conservative4way, simple, sum4way

# grid creation functions
from bin.utils.dims_func import set_dims
from bin.utils.coord_strch_func import coord_stretch
from bin.utils.geom_func import geom 
from bin.utils.mesh_func import mesh
from bin.utils.sangho_func import sangho
from bin.utils.metric_func import metric
from bin.utils.plot_mesh_func import plot_mesh

def init_from_file(self, grid_dim, input):
    """Performs conformal mapping and finds x,xc and vol values in physical space.
       
       Also plots the c-mesh/grid in phyiscal space
        
        Args:
            grid_dim (list):Number of cells in the x and y directions.
            input (dict):Dictionary containing data-file values
        """
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
    self.vol= Field((nx,ny))

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
    #plot_mesh(self)
    

def init_from_grid(newGrid, grid):
    """Makes a grid coarser and finds new x,xc and vol values on the coarse grid
       
       Also plots the coarser grid in physical space
        
        Args:
            grid (obj):finer input AirfoilMap object
            newGrid (obj):coarser output new AirfoilMap object
        """
    # transfer geometry info and set new dimensions
    newGrid.geo = grid.geo
    set_dims(newGrid)

    # get relavant values from old grid
    fields = grid.fields
    x = fields['x']
    xc = fields['xc']
    vol = fields['vol']

    # create new arrays
    [il, jl] = newGrid.get_size()
    [nx, ny] = [il-1, jl-1]
    xNew = Field((il,jl,2))
    xcNew = Field((nx,ny,2))
    volNew = Field((nx,ny))

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
    #plot_mesh(newGrid)