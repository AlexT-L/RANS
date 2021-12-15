import numpy as np
import scipy.interpolate as sc

def geom(self):
    #
    pi      = np.pi
    #use parameters from geo_param in input
    geo     = self.geo
    scal    = geo["scal"]
    xlim    = geo["xlim"]
    nn      = geo["nn"]
    xsing   = geo["xsing"]
    ysing   = geo["ysing"]
    slopt   = geo["slopt"]
    trail   = geo["trail"]#in radians, converted in AirfoilMap constructor
    
    #geo_var
    a       = self.a
    s0      = self.s0

    #x and y coords of aerfoil in physical space
    xn      = geo["xn"]
    yn      = geo["yn"]

    #mesh dimensions
    dim    = self.dims
    il     = dim["il"]
    jl     = dim["jl"]

    #coordinates of lower trailing edge and upper trailing edge in comp space
    itl     = self.itl
    itu     = self.itu

    # instantiate xs and xi arrays: airfoil coords after conformal mapping in computational space
    xs      =np.zeros(nn)
    ys      =np.zeros(nn)

    #set values of angls in array
    scal   = .50001*xlim**2/(xn[nn-1]  -xsing)

    angl    = 2*pi 
    u       = 1.
    v       = 0.

    for i in range(nn):
        xa        = xn[i]  -xsing #distance in x and y from the origin 
        ya        = yn[i]  -ysing
        angl      = angl  + np.arctan2((u*ya  -v*xa),(u*xa  +v*ya))
        r         = scal*np.sqrt(xa**2  +ya**2)
        u         = xa
        v         = ya
        r         = np.sqrt(2*r)

        #coordinates of airfoil after mapping to computational spac
        xs[i]     = r*np.cos(0.5*angl)
        ys[i]     = r*np.sin(0.5*angl)
    
    
    geo["scal"]  = 1/scal #set new value in dcitionary
    
    

    #fitting a cubic spline to aerfoil coords (in computational domain: after mapping) 
    # xs and ys to interpolate 
    # and make sure points on the computational mesh in the x-dir a[0]
    # match up witht the points on the aerfoil in the compuatational domain

    #fitting a cubic spline to aerfoil coords
    #input_vector=np.vstack((xs,ys)).T
    tck =  sc.splrep(xs,ys)
    
    # interpolating to make sure aerfoil has points on the mesh
    #fills in aerfoil values for s0 from itl to itu-1
    for i in range(itl-1,itu):
        s0[i]=sc.splev(a[0][i],tck)
    
    # import matplotlib.pyplot as plt
    # plt.figure(1)
    
    # #plt.plot(a[0][itl-1:itu],s0[itl-1:itu])
    # plt.plot(a[0],s0,"x")
    
    
    

    x1        = 0.125*(xlim)**2
    angl      = 2*pi
    u         = 1
    v         = 0
    i1        = 0
    i2        = itl  -1
    bnd       = np.array([itl,itu])
    
    

    #This for loop fills up values before and after the aerfoil coordinates in s0(lenght il) array
    #For i = 0 it fills in values from 0 to itl i.e. before the aerfoil in s0
    #For i = 1 it fills in values from itu to il i.e. before the aerfoil in s0
    for i in range(2):
        r         = np.sqrt(a[0][bnd[i]]**2  +s0[bnd[i]]**2)
        ang       = 2*np.arctan2(s0[bnd[i]],a[0][bnd[i]])
        r         = 0.5*r**2
        x0        = r*np.cos(ang)
        y0        = r*np.sin(ang)
        g         = slopt*(x0  -x1)
        b         = 1./(x0  -x1)
        
        #print("I1=",i1,"I2=",i2)
        for j in range(i1,i2):
            xa        = 0.5*a[0][j]**2
            d         = b*(xa  -x1)
            ya        = y0  +g*np.log(d)/d
            angl      = angl  +np.arctan2((u*ya  -v*xa),(u*xa  +v*ya))
            r         = np.sqrt(xa**2  +ya**2)
            u         = xa
            v         = ya
            r         = np.sqrt(2*r)
            s0[j]    = r*np.sin(0.5*angl)
        
        angl      = 0
        u         = 1
        v         = 0
        i1        = itu
        i2        = il
        
    #print("GEOM")

    # plt.figure(2)
    # plt.plot(a[0],s0)
    # plt.show()
    return
    


    