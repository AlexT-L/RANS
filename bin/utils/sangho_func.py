import numpy as np

def sangho(self):
    il          = self.il
    jl          = self.jl
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
    return
