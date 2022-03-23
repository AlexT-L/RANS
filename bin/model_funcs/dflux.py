# python implementation of dflux.f

from operator import pos
from bin.Field import Field, is_field, maximum, minimum, abs, pos_diff, isscalar, max
from bin.Grid import Grid
from bin.Workspace import Workspace
import numpy as np

def dflux(model, ws, state, dw, rfil):
    """
    calculate artificial dissipation fluxes on finest mesh using blended first and 
    third order fluxes
    
    Args:
        model (NavierStokes): physics model
        workspace (Workspace): contains the relevant Fields
        state (Field): density, x-momentum, y-momentum, and energy
        dw (Field): to store new residuals after completing fluxes 
        rfil (float): relaxation factor determining balance between viscous and artificial dissipation fluxes
        
    """

    # take a workspace ws and calculate dissipative fluxes

    # model parameters
    pad = model.padding
    n = Field.dim(state) # number of quantities being convected

    # necessary fields
    def mget(varName):
        return ws.get_field(varName, model.className)
    w = state # state vector
    fw   = mget('fw')
    p    = mget('p')   # pressure
    porI = mget('porI')
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
    vis2 = model.params['vis2']
    vis4 = model.params['vis4']

    fis2      = 2*rfil*vis2
    fis4      = rfil*vis4/16
    sfil      = 1  -rfil
    plim      = 0.001
    tol       = 0.25

    # working arrays
    dp = Field((nxp,nyp))
    dis2 = Field((nxp, nyp))
    dis4 = Field((nxp, nyp))
    d    = Field((nxp, nyp, n))
    e    = Field((nxp, nyp, n))
    TOL  = Field((nx+1, ny+1), tol)

    # c
    # c     replace the energy by the enthalpy
    # c
    w[:,:,3] = w[:,:,3] + p

    # c
    # c     dissipation in the i direction
    # c
    dp[1:ib, jp:je] = abs((p[ip:nxp, jp:je] - 2*p[1:ib, jp:je] + p[0:ie, jp:je]) / \
                          (p[ip:nxp, jp:je] + 2*p[1:ib, jp:je] + p[0:ie, jp:je] + plim))

    dp[0, jp:je] = dp[1, jp:je]
    dp[ib, jp:je] = dp[ie, jp:je]

    rad = minimum(radI[ip:ib, jp:je], radI[1:ie, jp:je])

    max1 = maximum(dp[ip+1:nxp, jp:je], dp[ip:ib, jp:je])
    max2 = maximum(dp[1:ie, jp:je], dp[0:il, jp:je])
    dpsafe = minimum(TOL[:,0:ny], maximum(max1, max2))
    dis2[1:ie, jp:je] = dpsafe*rad*fis2

    dis4[1:ie, jp:je] = rad*fis4
    dis4 = pos_diff(dis4, dis2)

    # forward differencing (1st order)
    d[0:ib, jp:je] = w[1:nxp, jp:je] - w[0:ib, jp:je]

    # central differencing (third order)
    e[1:ie, jp:je] = d[ip:ib, jp:je] - 2*d[1:ie, jp:je] + d[0:il, jp:je]

    gs = porI*(dis2[1:ie, jp:je]*d[1:ie, jp:je] - dis4[1:ie, jp:je]*e[1:ie, jp:je])

    fw[ip:ie, jp:je] = fw[ip:ie, jp:je]*sfil + gs[0:nx] - gs[1:il]

    # c
    # c     dissipation in the j direction
    # c
    if ny >= 3:
        dp[ip:ie, 1:jb] = abs((p[ip:ie, jp:nyp] - 2*p[ip:ie, 1:jb] + p[ip:ie, 0:je]) / \
                              (p[ip:ie, jp:nyp] + 2*p[ip:ie, 1:jb] + p[ip:ie, 0:je] + plim))

        
        rad = minimum(radJ[ip:ie, 2], radJ[ip:ie, 1])

        max1 = maximum(dp[ip:ie, 3], dp[ip:ie, 2])
        max2 = dp[ip:ie, 1]
        dis2[ip:ie, 1] = rad*fis2*minimum(TOL[0:nx,0], maximum(max1, max2))

        dis4[ip:ie, 1] = rad*fis4
        dis4[ip:ie, 1] = pos_diff(dis4[ip:ie, 1], dis2[ip:ie, 1])


        rad = minimum(radJ[ip:ie, je], radJ[ip:ie, jl])

        max1 = maximum(dp[ip:ie, je], dp[ip:ie, jl])
        max2 = dp[ip:ie, ny]
        dis2[ip:ie, jl] = rad*fis2*minimum(TOL[0:nx,0], maximum(max1, max2))

        dis4[ip:ie, jl] = rad*fis4
        dis4[ip:ie, jl] = pos_diff(dis4[ip:ie, jl], dis2[ip:ie, jl])

        rad = min(radJ[ip:ie, jp+1:je], radJ[ip:ie, jp:jl])

        max1 = maximum(dp[ip:ie, jp+2:jb], dp[ip:ie, jp+1:je])
        max2 = maximum(dp[ip:ie, jp:jl], dp[ip:ie, 1:ny])
        dis2[ip:ie, jp:jl] = rad*fis2*minimum(TOL[0:nx,0:ny-1], maximum(max1, max2))

        dis4[ip:ie, jp:jl] = rad*fis4
        dis4[ip:ie, jp:jl] = pos_diff(dis4[ip:ie, jp:jl], dis2[ip:ie, jp:jl])


        # forward difference (first order)
        d[ip:ie, 0:jb] = w[ip:ie, 1:nyp] - w[ip:ie, 0:jb]

        # central difference (third order)
        e[ip:ie, 1:je] = d[ip:ie, jp:jb] - 2*d[ip:ie, 1:je] + d[ip:ie, 0:jl]

        gs = porJ*(dis2[ip:ie, 1:je]*d[ip:ie, 1:je] - dis4[ip:ie, 1:je]*e[ip:ie, 1:je])

        fw[ip:ie, jp:je] += gs[:, 0:ny] - gs[:, 1:jl]

    # c
    # c     replace the enthalpy by the energy
    # c
    w[:,:,3] = w[:,:,3] - p

    # add to flux field
    dw += fw
