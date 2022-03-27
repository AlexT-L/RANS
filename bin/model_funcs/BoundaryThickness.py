"""This module calculates boundary layer thickness for viscosity.

    Libraries/Modules:
        numpy\n
        """
import numpy as np

def boundary_thickness(model, ws, state, ynot, dsti):
    """Calculates the boundary layer thickness.

    Notes: 
        Adapted from subroutine delt"""

    # necessary fields
    def mget(varName):
        return ws.get_field(varName, model.className)
    w = state
    x = ws.get_field('x')
    xc = mget('xc')

    # set geometry parameters
    PAD = model.padding
    [nx, ny] = ws.field_size()
    [nxp, nyp] = [PAD+nx+PAD, PAD+ny+PAD]
    [il, jl] = [nx+1, ny+1]
    [ie, je] = [nx+2, ny+2]
    [ib, jb] = [nx+3, ny+3]

    qs = np.ones(nxp)
    ut = np.ones(nxp)
    dn = np.ones(nxp)
    ssmax = np.ones(nxp)
    
    js = 0.75*(ny - 4)
    js = int(np.floor(js))

    for i in range(PAD,nx+PAD):
        qs[1]     = 0.0
        ut[1]     = 0.0
        j         = js
        xy        = 0.5*(x[i-1,j-1,0]  -x[i-1,j-2,0]+x[i-2,j-1,0]  -x[i-2,j-2,0])
        yy        = 0.5*(x[i-1,j-1,1]  -x[i-1,j-2,1]+x[i-2,j-1,1]  -x[i-2,j-2,1])
        qs[j]     = (yy*w[i,j,1]  -xy*w[i,j,2])/(w[i,j,0])
        si        = np.copysign(1,qs[js])

        for j in range(2,js+1):
            xy        = 0.5*(x[i-1,j-1,0]  -x[i-1,j-2,0]+x[i-2,j-1,0]  -x[i-2,j-2,0])
            yy        = 0.5*(x[i-1,j-1,1]  -x[i-1,j-2,1]+x[i-2,j-1,1]  -x[i-2,j-2,1])
            dsi       = 1.0/np.sqrt(xy**2  +yy**2)
            qs[j]     = si*(yy*w[i,j,1]  -xy*w[i,j,2])
            dn[j]     = 1.0/dsi
            ut[j]     = qs[j]*dsi

        dsti[i]   = 0.0
        ynot[i]   = 0.0
        ssmax[i]  = 0.0
        cdu       = 0.98
        jmax      = np.argmax(ut)
        fx        = 0.6*ut[jmax]
        lend      = 2
        lbig      = 2
        locke     = False

        for  j in range(3,js+1):
            if ( ut[j-1] < 0 and ut[j] >= 0):
                lbig = j

            if not locke:
                if ( ut[j-1] >= cdu*ut[j] and ut[j] > fx):
                    locke     = True
                    lend      = j

        uinf      = 1./ut[lend]
        for j in range(lbig,lend+1):
            dsti [i]  = dsti[i] + (ut[lend]*dn[j] - qs[j])
            ra        = w[i,lend,0]/w[i,j,0]
            ssmax[i]  = ssmax[i] + ra*ut[j]*uinf*(dn[j]-qs[j]*uinf)

        dsti[i]   = dsti[i]*uinf
        
        dsti[i]  = max(dsti[i],1.0e-6)
        ra        = w[i,lend,0]/w[i,2,0]
        ssmax[i]  = max(ssmax[i],ra*qs[1]*uinf)
        fc        = .95*ut[lend]

        jse       = lend
        for j in range(2,jse+1):
            lend      = j  -1
            if (ut[j] > fc):
                break

        xbi       = 0.5*(x[i-1,0,0]+x[i-2,0,0])
        ybi       = 0.5*(x[i-1,0,1]+x[i-2,0,1])
        ycorr     = np.sqrt((xc[i,lend,0] - xbi)**2+(xc[i,lend,1]-ybi)**2)
        ynot[i]   = 1.5*(ycorr  +dn[lend]*(fc  -ut[lend])/(ut[lend+1]  -ut[lend]))
    return
