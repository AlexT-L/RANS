"""This module creates the c-mesh in physical space after conformal mapping

    Libraries/Modules:
        numpy\n
        Field\n
        """
import numpy as np
from bin.Field import isfinite

def mesh(self):
    """This function maps back to physical space to create the c-mesh
       
       First a cubic spline interpolation is performed to make sure that points on airfoil geometry line up with a0 points
       
       Then a0,b0,xs and ys are mapped to an s0 array and this is used to creat the x 
       array which contains the vertices of the right hand side corners of all cells 
       in the physical domain.
    """
    #vertices x[i,j,1]=x vertex and x[i,j,2]=y vertex in physical space
    x       =self.x

    #mesh dimensions
    dim    = self.dims
    il     = dim["il"]
    jl     = dim["jl"]

    #geo_var
    a       = self.a
    b0      = self.b0
    s0      = self.s0

    #use param from geo_param
    geo     = self.geo
    scal    = geo["scal"]
    xsing   = geo["xsing"]
    ysing   = geo["ysing"]    

    x0        = xsing/scal
    y0        = ysing/scal

    for  j in range(jl):

      s1    = 1 - (j  -1)/(jl  -1)**2

      for i in range(il):
          yp        = a[1][i]*b0[j]  +s1*s0[i]
          x[i,j,0]  = x0  +.5*(a[0][i]*a[0][i]  -yp*yp)
          x[i,j,1]  = y0  +a[0][i]*yp
    
    print("MESH")
    assert(isfinite(x))
    # import matplotlib.pyplot as plt
    # plt.plot(x[:,:,0],x[:,:,1],'o')
    # plt.show()
    return