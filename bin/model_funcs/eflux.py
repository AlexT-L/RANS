# python implementation of eflux.f

from Field import Field
from Grid import Grid
from Workspace import Workspace
from NavierStokes import NavierStokes
import numpy as np

def eflux(ws,dw):
    # take a workspace ws and calculate convective fluxes

    G = ws.grd # grab grid
    w = ws.flds['w'] # state vector
    porJ = ws.flds['porJ'] # porosity
    P = ws.flds['P']   # pressure
    n = 4 # number of quantities being convected
    
    # i direction
    fs = np.zeros([G.ib, G.jb,n])
    for j in range(1,G.jl):
        for i in range(0,G.il):

            # normal vector 
            dxy = G.X[i,j,0] - G.X[i,j-1,0]
            dyy = G.X[i,j,1] - G.X[i,j-1,1]
            
            # pressure averaging
            Pa = P[i+1][j] + P[i][j]

            # flux operator
            qsp       = (dyy*w[i+1,j,1]  -dxy*w[i+1,j,2])/w[i+1,j,0]
            qsm       = (dyy*w[i,j,1]  - dxy*w[i,j,2])/w[i,j,0]

            # add up on faces
            fs[i,j,0] = qsp*w[i+1,j,0]   + qsm*w[i,j,0] # density
            fs[i,j,1] = qsp*w[i+1,j,1]  + qsm*w[i,j,1]  + dyy*Pa # x - momentum
            fs[i,j,2] = qsp*w[i+1,j,2]  + qsm*w[i,j,2]  - dxy*Pa # y - momentum
            fs[i,j,3] = qsp*(w[i+1,j,3] + P[+1,j]) + qsm*(w[i,j,3] + P[i,j]) # energy

    # now add everything up
    for j in range(1,G.jl):
        for i in range(1,G.il):
            dw[i,j,:] = fs[i,j,:] -fs[i-1,j,:]


    # j direction
    for j in range(0,G.jl):
      for i in range(1,G.il):

          # normal vector
         dxx        = G.X[i,j,0]  -G.X[i-1,j,0]
         dyx        = G.X[i,j,1]  -G.X[i-1,j,1]
         # pressure average
         Pa        = P[i,j+1]  +P[i,j]
         # convective operator
         qsp       = porJ[i,j]*(dxx*w[i,j+1,2]  - dyx*w[i,j+1,1])/w[i,j+1,0]    
         qsm       = porJ[i,j]*(dxx*w[i,j,2]  - dyx*w[i,j,1])/w[i,j,0]
        # add up on faces
         fs[i,j,0] = qsp*w[i,j+1,0]  +qsm*w[i,j,0]
         fs[i,j,1] = qsp*w[i,j+1,1]  +qsm*w[i,j,1]  - dyx*Pa
         fs[i,j,2] = qsp*w[i,j+1,2]  +qsm*w[i,j,2]  + dxx*Pa
         fs[i,j,3] = qsp*(w[i,j+1,3]  +P[i,j+1]) +qsm*(w[i,j,3]  + P[i,j])
                   
    # now add everything up for j direction
    for j in range(1,G.jl):
        for i in range(1,G.il):
            dw[i,j,:] = dw[i,j,:] + fs[i,j,:] - fs[i,j-1,:]
    