# python implementation of eflux.f

from Field import Field, max, minimum
from Grid import Grid
from Workspace import Workspace
import numpy as np

def dfluxc(model, ws, state, dw, rfil):
    # take a workspace ws and calculate dissipative fluxes

    # model parameters
    pad = model.padding
    n = state.dim() # number of quantities being convected

    # necessary fields
    def mget(varName):
        return ws.get_field(varName, model.className)
    w = state # state vector
    fw   = mget('fw')
    p    = mget('p')   # pressure
    porJ = mget('porJ') # porosity
    radI = mget('radI')
    radJ = mget('radJ')

    # grid dimensions
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

    # physical paramteters
    vis0 = model.params['vis0']

    fis0      = rfil*abs(vis0)/8.
    sfil      = 1.  -rfil

    # c
    # c     dissipation in the i direction
    # c
    dis = fis0*minimum(radI[ip:ib, jp:je], radI[1:ie, jp:je])

    fs  = dis*(w[ip:ib, jp:je] - w[1:ie, jp:je])
    fs += dis*(p[ip:ib, jp:je] - p[1:ie, jp:je])

    fw[:] = sfil*fw - fs[ip:ie, jp:je] + fs[1:il, jp:je]

    if ny < 3:
        return

    # c
    # c     dissipation in the j direction
    # c
    dis = fis0*porJ*minimum(radI[ip:ie, jp:jb], radI[ip:ie, 1:je])

    fs  = dis*(w[ip:ie, jp:jb] - w[ip:ie, 1:je])
    fs += dis*(p[ip:ie, jp:jb] - p[ip:ie, 1:je])

    fw += fs[ip:ie, jp:je] + fs[ip:ie, 1:jl]
