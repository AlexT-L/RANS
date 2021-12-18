# python implementation of eflux.f

from bin.Field import Field, max, minimum
from bin.Grid import Grid
from bin.Workspace import Workspace
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

    fis0      = rfil*np.abs(vis0).item()/8
    sfil      = 1.  -rfil

    # working array
    fs = Field((nxp, nyp, n))

    # c
    # c     dissipation in the i direction
    # c
    dis = fis0*minimum(radI[ip:ib, jp:je], radI[1:ie, jp:je])
    # dis = Field(dis)

    print("\ndfluxc")
    print(dis.shape)
    print(type(minimum(radI[ip:ib, jp:je], radI[1:ie, jp:je]).__module__))
    print(minimum(radI[ip:ib, jp:je], radI[1:ie, jp:je]).vals)
    print(type(radI).__module__)
    print(type(dis).__module__)
    print(type(dis.vals).__module__)
    print(dis.shape())

    fs[1:ie, jp:je]     = dis*(w[ip:ib, jp:je] - w[1:ie, jp:je])
    fs[1:ie, jp:je, 3] += dis*(p[ip:ib, jp:je] - p[1:ie, jp:je])

    fw[ip:ie, jp:je] = sfil*fw[ip:ie, jp:je] - fs[ip:ie, jp:je] + fs[1:il, jp:je]

    if ny < 3:
        return

    # c
    # c     dissipation in the j direction
    # c
    dis = fis0*porJ*minimum(radI[ip:ie, jp:jb], radI[ip:ie, 1:je])

    fs[ip:ie, 1:je]     = dis*(w[ip:ie, jp:jb] - w[ip:ie, 1:je])
    fs[ip:ie, 1:je, 3] += dis*(p[ip:ie, jp:jb] - p[ip:ie, 1:je])

    fw[ip:ie, jp:je] += fs[ip:ie, jp:je] + fs[ip:ie, 1:jl]

    # add to dw
    dw += fw