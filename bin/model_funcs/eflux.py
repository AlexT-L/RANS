# python implementation of eflux.f

from bin.Field import Field, max
from bin.Grid import Grid
from bin.Workspace import Workspace
import numpy as np

def eflux(model, ws, state, dw):
    """
    calculate convective fluxes 
    
    Args:
        model (NavierStokes): physics model
        workspace (Workspace): the relevant fields
        state (Field): containing the density, x-momentum, y-momentum, and energy
        dw (Field): to store new residuals after completing fluxes
        
    """
    
    # take a workspace ws and calculate convective fluxes

    pad = model.padding
    w = state # state vector
    porJ = ws.get_field('porJ', model.className) # porosity
    p = ws.get_field('p', model.className)   # pressure
    n = Field.dim(state) # number of quantities being convected

    [nx, ny] = ws.field_size()
    [nxp, nyp] = [pad+nx+pad, pad+ny+pad]
    ip = pad
    jp = pad
    il = nx+1
    jl = ny+1
    ie = nx+2
    je = ny+2
    ib = nx+3
    jb = ny+3

    # flux array
    fs = Field((ib+1, jb+1, n))

    # i direction
    dx = ws.edge_normals(0)
    dyx = dx[:,:,0]
    dyy = dx[:,:,1]
    p_avg = p[1:ie, jp:je] + p[ip:ib, jp:je]

    # flux operator
    qsp = (dyy*w[ip:ib, jp:je, 1] - dyx*w[ip:ib, jp:je, 2]) / w[ip:ib, jp:je, 0]
    qsm = (dyy*w[ 1:ie, jp:je, 1] - dyx*w[ 1:ie, jp:je, 2]) / w[ 1:ie, jp:je, 0]

    # add up on faces
    fs[1:ie, jp:je, 0] = qsp * w[ip:ib, jp:je, 0]                   + qsm*w[1:ie, jp:je, 0] # density
    fs[1:ie, jp:je, 1] = qsp * w[ip:ib, jp:je, 1]                   + qsm*w[1:ie, jp:je, 1]  + dyy*p_avg # x - momentum
    fs[1:ie, jp:je, 2] = qsp * w[ip:ib, jp:je, 2]                   + qsm*w[1:ie, jp:je, 2]  - dyx*p_avg # y - momentum
    fs[1:ie, jp:je, 3] = qsp *(w[ip:ib, jp:je, 3] + p[ip:ib, jp:je]) + qsm*(w[1:ie, jp:je, 3] + p[1:ie, jp:je]) # energy

    # now add everything up
    dw[ip:ie, jp:je, :] = fs[ip:ie, jp:je, :] - fs[1:il, jp:je, :]


    # j direction
    dx = ws.edge_normals(1)
    dxx = dx[:,:,0]
    dxy = dx[:,:,1]
    p_avg = p[ip:ie, 1:je] + p[ip:ie, jp:jb]

    # flux operator
    qsp = porJ * (dxx*w[ip:ie, jp:jb, 1] - dxy*w[ip:ie, jp:jb, 2]) / w[ip:ie, jp:jb, 0]
    qsm = porJ * (dxx*w[ip:ie,  1:je, 1] - dxy*w[ip:ie,  1:je, 2]) / w[ip:ie,  1:je, 0]

    # add up on faces
    fs[ip:ie, 1:je, 0] = qsp*w[ip:ie, jp:jb, 0]                     + qsm*w[ip:ie, 1:je, 0] # density
    fs[ip:ie, 1:je, 1] = qsp*w[ip:ie, jp:jb, 1]                     + qsm*w[ip:ie, 1:je, 1]  + dxy*p_avg # x - momentum
    fs[ip:ie, 1:je, 2] = qsp*w[ip:ie, jp:jb, 2]                     + qsm*w[ip:ie, 1:je, 2]  - dxx*p_avg # y - momentum
    fs[ip:ie, 1:je, 3] = qsp*(w[ip:ie,jp:jb, 3] + p[ip:ie, jp:jb]) + qsm*(w[ip:ie, 1:je, 3] + p[ip:ie, 1:je]) # energy

    # now add everything up for j direction
    dw[ip:ie, jp:je, :] += fs[ip:ie, jp:je, :] - fs[ip:ie, 1:jl, :]

    stop = 0