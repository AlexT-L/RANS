"""This module creates x and y array in computational space.
   
    Libraries/Modules:
        numpy\n
        """
import numpy as np

########################################################################################
# Coordinate stretching function: create array a0 and b0 in x and y dirn computationally 
#                                 with spacing such that when it is mapped back to physical
#                                 domain it points on the mesh are evenly spaced for a given
#                                 i/j direction
#########################################################################################

def coord_stretch(self):
    """It create array a0 and b0 in x and y dirn computationally\n
       with spacing such that when it is mapped back to physical\n
       domain it points on the mesh are evenly spaced for a given\n
       i/j direction."""
    #use parameters from geo_param in input
    geo    = self.geo
    xte    = geo["xte"]
    boundx = geo["boundx"]
    boundy = geo["boundy"]
    bunch  = geo["bunch"]
    ylim1  = geo["ylim1"]
    ylim2  = geo["ylim2"]
    ax     = geo["ax"]
    ay     = geo["ay"]
    sy     = geo["sy"]

    #geo_var
    a      = self.a
    b0     = self.b0
    s0     = self.s0

    #parameters from dims
    dim    = self.dims
    il     = dim["il"]
    jl     = dim["jl"]
    nx     = dim["nx"]
    ny     = dim["ny"]
    
    #coord stretching
    ile         = self.ile
    pi          = np.pi
    xlim        = xte*boundx
    geo["xlim"] = xlim #add xlim to geo_param dict to use in geom function
    dx          = 2.0*boundx/nx

    px     = pi/xlim
    bp     = bunch/px

    a2     = 3.0*ylim1 - 4.0*ylim2
    a3     = 2.0*ylim1 - 3.0*ylim2

    for i in range(il):
        d  = ((i+1) -ile)*dx
        if abs(d) <= xlim:
            d = d  +bp*np.sin(px*d)
        else:
            b         = 1.
            if d < 0.0:
                b = -1
            g  =1.0  -((d  -b*xlim)/(1.0  -xlim))**2 
            c  = g**ax
            d  = b*xlim  +(1.0  -bunch)*(d  -b*xlim)/c
        a[0][i] = d #define x-array of computational mesh 
        d     = abs(d/xlim)
        if d >= 1.0:
            a[1][i]     = ylim2*xlim/abs(a[0][i])
        else:
            a[1][i]     = ylim1  -d*d*(a2  -a3*d)

 

    dy        = boundy/(jl  -1)
    for i in range(jl):
        d         = ((i+1)  -1)*dy
        g         = 1.  -d*d
        c         = g**ay
        b0[i]     = sy*d/c
    

    return
      
      

    

