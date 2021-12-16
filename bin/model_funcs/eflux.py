# python implementation of eflux.f

from Field import Field, max
from Grid import Grid
from Workspace import Workspace
import numpy as np

def eflux(model, ws, state, dw):
    # take a workspace ws and calculate convective fluxes

    pad = model.padding

    def edge(i, j, side):
        return ws.edge(i-pad, j-pad, side)
    def normal(i, j, side):
        return ws.edge_normal(i-pad, j-pad, side)

    w = state # state vector
    porJ = ws.get_field('porJ', model.className) # porosity
    p = ws.get_field('p', model.className)   # pressure
    x = ws.get_field('x')
    n = state.dim() # number of quantities being convected

    [nx, ny] = ws.field_size()
    il = nx+1
    jl = ny+1
    ie = nx+2
    je = ny+2
    ib = nx+3
    jb = ny+4
    
    # i direction
    fs = Field((ib+1, jb+1), n)
    for j in range(pad,ny+pad):
        for i in range(pad-1,nx+pad):

            # normal vector 
            [dxy, dyy] = normal(i, j, 'e')
            
            # pressure averaging
            pa = p[i+1-pad,j-pad] + p[i-pad,j-pad]

            # flux operator
            qsp       = (dyy*w[i+1,j,1]  -dxy*w[i+1,j,2]) / w[i+1,j,0]
            qsm       = (dyy*w[i,j,1]    -dxy*w[i,j,2])    / w[i,j,0]

            # add up on faces
            fs[i,j,0] = qsp*w[i+1,j,0]  + qsm*w[i,j,0] # density
            fs[i,j,1] = qsp*w[i+1,j,1]  + qsm*w[i,j,1]  + dyy*pa # x - momentum
            fs[i,j,2] = qsp*w[i+1,j,2]  + qsm*w[i,j,2]  - dxy*pa # y - momentum
            fs[i,j,3] = qsp*(w[i+1,j,3] + p[+1,j]) + qsm*(w[i,j,3] + p[i,j]) # energy

    # now add everything up
    # for j in range(pad,pad+ny):
    #     for i in range(pad,pad+nx):
    #         dw[i,j,:] = fs[i,j,:] -fs[i-1,j,:]
    dw[pad:nx+pad, pad:ny+pad, :] = fs[pad:nx+pad, pad:ny+pad, :] - fs[1:nx+1, pad:ny+pad, :]


    # j direction
    for j in range(pad,pad+ny):
      for i in range(pad-1,pad+nx):

        # normal vector
        [dxx, dyx] = normal(i, j, 'n')

        # pressure average
        pa        = p[i,j+1]  +p[i,j]

        # convective operator
        qsp       = porJ[i-pad,j-pad] * (dxx*w[i,j+1,2]  - dyx*w[i,j+1,1]) / w[i,j+1,0]    
        qsm       = porJ[i-pad,j-pad] * (dxx*w[i,j,2]    - dyx*w[i,j,1])   / w[i,j,0]

        # add up on faces
        fs[i,j,0] = qsp*w[i,j+1,0]  +qsm*w[i,j,0]
        fs[i,j,1] = qsp*w[i,j+1,1]  +qsm*w[i,j,1]  - dyx*pa
        fs[i,j,2] = qsp*w[i,j+1,2]  +qsm*w[i,j,2]  + dxx*pa
        fs[i,j,3] = qsp*(w[i,j+1,3]  +p[i,j+1]) +qsm*(w[i,j,3]  + p[i,j])
                   
    # now add everything up for j direction
    # for j in range(pad,ny+pad):
    #     for i in range(pad,nx+pad):
    #         dw[i,j,:] = dw[i,j,:] + fs[i,j,:] - fs[i,j-1,:]
    dw[pad:nx+pad, pad:ny+pad, :] += fs[pad:nx+pad, pad:ny+pad, :] - fs[pad:nx+pad, 1:ny+1, :]