import numpy as np
from numpy.core.numeric import isclose
from numpy.lib.function_base import iterable
from Field import Field
from Grid import Grid
from Input import Input

        
class AirfoilMap(Grid):
    
    def __init__(self, input):

        # let Grid contain the variables in input.dims, geo_param, flo_param and in_var
        self.dims = input.dims
        self.in_var=input.in_var        
        nx = self.dims['nx']
        ny = self.dims['ny']
        
        
        #modify input object by adding values to flo_param dictionary
        input.flo_param["scal"]=0
        
        #append to geo_param dictionary
        self.geo=input.geo_param
        geo =self.geo
        geo["trail"]=geo["trail"]*np.pi/180#convert trail angle to radians
        geo["nbl"] = 0
        geo["jlinv"] = 0
        geo["nyinv"] = 0

        xte=self.geo['xte']

        self.flo=input.flo_param
        kvis=self.flo["kvis"]
        

        # set mesh dimensions and update dims
        dim         = input.dims
        dim["il"]   = int(nx + 1) # number of points/edges in i dir of computational grid
        dim["jl"]   = int(ny + 1) # number of points/edges in j dir of computational grid
         # values below define number of points in
         # computational domain for a padded grid
        dim["ie"]   = int(nx + 2) 
        dim["je"]   = int(ny + 2) 
        dim["ib"]   = int(nx + 3)
        dim["jb"]   = int(ny + 3)
        self.il     = dim["il"]
        self.jl     = dim["jl"]
        il          = dim["il"]
        jl          = dim["jl"]
        ib          = dim["ib"]
        jb          = dim["jb"]

        

        # initialize x-y vertex, center,vol and porosity arrays
        self.x  = np.zeros((il,jl,2))
        self.xc = np.zeros((ib+1,jb+1,2))
        self.vol= np.zeros((ib+1,jb+1))


        #set the limits of the aerfoil profile
        self.ite       = int(0.5*xte*nx) #coordinate of trailing edge in physical space
        self.ile       = np.floor(il/2  +1) #coordinate of leading edge in physical space
        self.itl       = int(self.ile - self.ite) #lower coordinate of trailing edge in computational space
        self.itu       = int(self.ile + self.ite) #upper coordinate of trailing edge in computational space
        
        #param dictionaries needed
        self.dims = input.dims
        self.flo_param=input.flo_param
        self.geo_param=input.geo_param
        
        
        self.a  = np.array([np.zeros(il),np.zeros(il)])
        self.b0 = np.zeros(jl)
        self.s0 = np.zeros(il)

        #define point distributions in each coordinate direction (coordinate stretching)
        self.coord_stretch()

        #sqrt root mapping of aerfoil profile to slit & interpolating 
        # to make sure points on aerfoil match those of the grid
        self.geom()
        
        
        #making mesh 
        self.mesh()

        #sanghos modification for non-dimensionalization(re-scaling the c-mesh)
        self.sangho()

        #plot mesh
        self.plot_mesh(self.x)
        
        

    #class methods
    from coord_strch_func import coord_stretch
    from geom_func import geom 
    from mesh_func import mesh
    from sangho_func import sangho
    from plot_mesh_func import plot_mesh

   
    def get_size(self):
         return [self.ib, self.jb]




# input=Input("rae1-s1.data")
# print("INPUT")
# grid = AirfoilMap(input)
# print("SQRT")




