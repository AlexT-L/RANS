import numpy as np
from numpy.core.numeric import isclose
from numpy.lib.function_base import iterable
from Field import Field
from Grid import Grid
from Input import Input

        
class AirfoilMap():
    
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
        
        #set the limits of the outer mesh for a viscous simulation (kvis>1)
        if kvis > 1:
            geo["nbl"]= jl/2
            dim["ny"]   = ny  -geo["nbl"]
            dim["jl"]   = jl  -geo["nbl"]
        #geo_var (array of variables required for sqrt mapping)
        
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
        itl         = self.itl
        itu         = self.itu
        x           = self.x

        xmax        = self.x[itl,0,0]
        xmin        = self.x[itl,0,0]

        for i in range(itl-1,itu):
            xmin       = min(xmin,x[i,1,1])


        scal       = xmax  -xmin

        for  i in range(il):
            for j in range(jl):
                x[i,j,0]   = x[i,j,0]/scal
                x[i,j,1]   = x[i,j,1]/scal
        

    #class methods
    from coord_strch_func import coord_stretch
    from geom_func import geom 
    from mesh_func import mesh

        # #insert inner sub-layer for viscous simulations
        # if kvis > 0:
        #     jlinv     = jl
        #     nyinv     = ny
        #     ny        = ny  +nbl
        #     jl        = jl  +nbl
        #     call vmesh

      






    #     # it will also contain the variables in mesh_var as members
    #     # these fields will require a lot more math but I barely 
    #     # even know how to spell conformal mapping so I leave that to someone else

    #     # rectangular cartesian grid for now
    #     x_vec = np.linspace(dims['x_bound'][0], dims['x_bound'][1], self.ib)
    #     y_vec = np.linspace(dims['y_bound'][0], dims['y_bound'][1], self.jb)
    #     xg, yg = np.meshgrid(x_vec, y_vec)
        
    #     # physical vertex locations
    #     self.X = Field(self.get_size(), stateDim=2) 
    #     self.X.vals[:,:,0] = xg.T
    #     self.X.vals[:,:,1] = yg.T

    #     # cell volumes, porosity, and far field mask
    #     self.Vol = Field(self.get_size())
    #     self.PorJ = Field(self.get_size())
    #     self.PorI = Field(self.get_size())
    #     self.Fint = Field(self.get_size())

    # def get_size(self):
    #     return [self.ib, self.jb]


#trial
import matplotlib.pyplot as plt 
input=Input("rae9e-s3.data")
print("INPUT")
grid = AirfoilMap(input)
print("SQRT")
ver =grid.x
x=ver[:,:,0]
y=ver[:,:,1]
plt.plot(x,y)
plt.plot(input.in_var["xn"],input.in_var["yn"],"+")
plt.show()