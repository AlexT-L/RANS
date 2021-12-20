"""This module calculates boundary layer thickness for viscosity.

    Libraries/Modules:
        numpy\n
        """
import numpy as np

def boundary_thickness(model, ws, state):
    """Calculates the boundary layer thickness.

    Notes: 
        Adapted from subroutine delt"""

    # necessary fields
    def mget(varName):
        return ws.get_field(varName, model.className)
    w = state
    x = mget('x')
    xc = mget('xc')

    # defining local variables:
    [nx, ny] = ws.field_size()
    il = nx+1

    qs = np.ones(nx)
    ut = np.ones(nx)
    dn = np.ones(nx)
    dsti = np.ones(nx)
    ynot = np.ones(nx)
    ssmax = np.ones(nx)
    
    js = 0.75*(ny - 4)
    js = int(np.floor(js))

    for i in range(1,il):
        qs[0]     = 0.
        ut[0]     = 0.
        j         = js
        xy        = .5*(x[i,j,0]  -x[i,j-1,0]+x[i-1,j,0]  -x[i-1,j-1,0])
        yy        = .5*(x[i,j,1]  -x[i,j-1,1]+x[i-1,j,1]  -x[i-1,j-1,1])
        qs[j]     = (yy*w[i,j,1]  -xy*w[i,j,2])/(w[i,j,0])
        si        = np.copysign(1,qs[js])

        for j in range(1,js):
            xy        = .5*(x[i,j,0]  -x[i,j-1,1]+x[i-1,j,0]  -x[i-1,j-1,0])
            yy        = .5*(x[i,j,1]  -x[i,j-1,2]+x[i-1,j,1]  -x[i-1,j-1,1])
            dsi       = 1./np.sqrt(xy**2  +yy**2)
            qs[j]     = si*(yy*w[i,j,1]  -xy*w[i,j,2])
            dn[j]     = 1./dsi
            ut[j]     = qs[j]*dsi

        dsti[i]   = 0.
        ynot[i]   = 0.
        ssmax[i]  = 0.
        cdu       = .98
        jmax      = np.argmax(ut)
        fx        = .6*ut[jmax]
        lend      = 2
        lbig      = 2
        locke     = False

        for  j in range(2,js):
            if ( ut[j-1] < 0 and ut[j] >= 0):
                lbig = j

            if not locke:
                if ( ut[j-1] >= cdu*ut[j] and ut[j] > fx):
                    locke     = True
                    lend      = j

        uinf      = 1./ut[lend]
        for j in range(lbig-1,lend):
            dsti [i]  = dsti[i] + (ut[lend]*dn[j] - qs[j])
            ra        = w[i,lend,0]/w[i,j,0]
            ssmax[i]  = ssmax[i] + ra*ut[j]*uinf*(dn[j]-qs[j]*uinf)

        dsti[i]   = dsti[i]*uinf
        
        dsti[i]  = max(dsti[i],1.e-6)
        ra        = w[i,lend,0]/w[i,1,0]
        ssmax[i]  = max(ssmax[i],ra*qs[1]*uinf)
        fc        = .95*ut[lend]

        jse       = lend
        for j in range(1,jse):
            lend      = j  -1
            if (ut[j] > fc):
                break

        xbi       = .5*(x[i,0,0]+x[i-1,0,0])
        ybi       = .5*(x[i,0,1]+x[i-1,0,1])
        ycorr     = np.sqrt((xc[i,lend,0] - xbi)**2+(xc[i,lend,1]-ybi)**2)
        ynot[i]   = 1.5*(ycorr  +dn[lend]*(fc  -ut[lend])/(ut[lend+1]  -ut[lend]))
    return
