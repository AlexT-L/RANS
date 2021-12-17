import numpy as np

# set dimesions
def set_dims(self):
    # get dimensions
    [nx, ny] = self.divisions

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
    
    
    #set the limits of the aerfoil profile
    xte            = self.geo['xte']
    self.ite       = int(0.5*xte*nx) #coordinate of trailing edge in physical space
    self.ile       = np.floor(il/2  +1) #coordinate of leading edge in physical space
    self.itl       = int(self.ile - self.ite) #lower coordinate of trailing edge in computational space
    self.itu       = int(self.ile + self.ite) #upper coordinate of trailing edge in computational space

    # store in dims
    dim['ite'] = self.ite + 1
    dim['ile'] = self.ile + 1
    dim['itl'] = self.itl + 1
    dim['itu'] = self.itu + 1
