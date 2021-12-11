import numpy as np
from Grid import Grid

class BoundaryThickness():
    def boundary_thickness(params, dims):
        # from subroutine delt
        # # calculates the boundary layer thickness
        # uses dims, flo_var, mesh_var, solv_var, flo_param, solv_param

        # variabes from inputs/dims/params
        ny = dims['ny']
        il = dims['il']

        x = params['x']
        w = params['w']
        xc = params['xc']

        # defining local variables:
        dim_var = 10
        qs = np.ones(dim_var)
        ut = np.ones(dim_var)
        dn = np.ones(dim_var)
        dsti = np.ones(dim_var)
        ynot = np.ones(dim_var)
        ssmax = np.ones(dim_var)
        




        # js        = 2*jl/3  -2
        js        = .75*(ny  -4)
        js = int(np.round(js))

        # one big loop 
        for i in range(1,il):
            qs[0]     = 0.
            ut[0]     = 0.
            j         = js
            xy        = .5*(x[i,j,0]  -x[i,j-1,0]+x[i-1,j,0]  -x[i-1,j-1,0])
            yy        = .5*(x[i,j,1]  -x[i,j-1,1]+x[i-1,j,1]  -x[i-1,j-1,1])
            qs[j]     = (yy*w[i,j,1]  -xy*w[i,j,2])/(w[i,j,0])
            # replacing sign function in fortran with np.copysign, should be the same
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
            # replacing function idmax(js,ut,1) with argmax
            # purpose is to "find the index of element having max value"
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
dim_var = 10
params = {
  "ie": dim_var,
  "je": dim_var,
  "kvis": 2,
  "gamma": 1,
  "rm": 1,
  "re": 1,
  "ncyc": dim_var,
  "rev": np.ones((dim_var+1,dim_var+1)),
  "cmesh": 1,
  "ncyci1": -1,
  "itl": dim_var-2, 
  "itu": dim_var-2,
  "x": np.random.rand(dim_var,dim_var,3),
  "w": np.ones((dim_var,dim_var,3)),
  "p": np.ones((dim_var,dim_var)),
  "vol": np.ones((dim_var,dim_var)),
  "xtran": 0,
  # in Visc but was not needed in BL:
  "scal": 1,
  "chord": 1,
  "t0": 1,
  "rmu0": 1,
  "p" : np.ones((dim_var,dim_var)),
  "mode": 0,
  "kturb": 1,
  "xc": np.ones((dim_var,dim_var,3))*2,
  "ynot": np.ones(dim_var),
  "rlv": np.ones((dim_var,dim_var)),
  "dsti": np.ones(dim_var),
  "ib": 1
  # 
}
dims = {
    "il": dim_var - 1, 
    "jl": dim_var - 1,
    "ny": dim_var,
}
BoundaryThickness.boundary_thickness(params, dims)