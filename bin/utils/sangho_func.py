"""This module scales the values of the vertices

    Libraries/Modules:
        numpy\n
"""
import numpy as np

def sangho(self):
    """ Performs sanghos non-dimensionalisation method
    
        Scales the vertices in x to wrap around the airfoil in physical space
    """
    il          = self.il
    jl          = self.jl
    itl         = self.itl
    itu         = self.itu
    x           = self.x

    xmax        = x[itl,0,0]
    xmin        = x[itl,0,0]

    for i in range(itl,itu+1):
        xmin       = min(xmin,x[i,1,0])


    chord       = xmax  -xmin

    for  i in range(il):
        for j in range(jl):
            x[i,j,0]   = x[i,j,0]/chord
            x[i,j,1]   = x[i,j,1]/chord
            
    self.geo['chord'] = 1.0
    self.geo['xm'] = xmin + 0.25
    self.geo['ym'] = x[itl,0,1]
    return
