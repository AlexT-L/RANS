import sys
sys.path.append('../RANS/bin')
import time

import numpy as np
# from Grid import Grid
from BaldwinLomax import turbulent_viscosity
from BoundaryThickness import boundary_thickness

# @profile
def compute_viscosity(params, dims):
    z=0
    '''
     from subroutine viscf.f
    '''
    '''
    computes viscosity coefficients  
    '''             
    ie = params['ie'] # Mesh dimension
    je = params['je'] # Mesh dimension
    kvis = params['kvis']
    gamma = params['gamma']
    rm = params['rm']
    re = params['re']
    ncyc = params['ncyc']

    il = dims['il']
    jl = dims['jl']
    itl = params['itl']
    itu = params['itu']
    w = params['w']
    xtran = params['xtran'] # needs to be from flo_param

    # other parameters needed in this 
    scal = params['scal']
    chord = params['chord']
    t0 = params['t0']
    rmu0 = params['rmu0']
    p = params['p']
    mode = params['mode']
    kturb = params['kturb']
    rev = params['rev']
    x = params['x']
    xc = params['xc']
    ynot = params['ynot']
    rlv = params['rlv']
    dsti = params['dsti']
    ib = params['ib']

    # initializing, defined later
    rev0 = np.ones((dim_var,dim_var))
    astr = np.ones((dim_var,dim_var))
    u = np.ones((dim_var,dim_var))
    v = np.ones((dim_var,dim_var))

    # useful constants
    # ckr       = (.062/(2.*pi)**4)
    # ckr       = .01915
    ckr       = .0256
    # cwk       = .225
    cwk       = 0.
    #  scf       = re/(np.sqrt(gamma)*rm)
    scf       = (scal*re/chord)/(np.sqrt(gamma)*rm)
    visc_const_C1 = 1.461e-06
    visc_const_C2 = 110.3

    # compute the molecular viscosity
    # i>1:il, i+1>2:il+1, i-1>0:il-1
    # j>1:jlm, j+1>2:jlm+1, j-1>0:jlm-1

    tt       = p[0:ie,0:je]/w[0:ie,0:je,0]*t0
    rlv[0:ie,0:je] = visc_const_C1*tt*np.sqrt(tt)/((tt+visc_const_C2)*rmu0)


    # for laminar flows we are done.
    # for turbulent flows we are also done on the coarser grids.

    if (kvis <= 1) or (mode!=0):
        return
    # if we are using the baldwin and lomax model call turbbl and return

    aturb     = 1.
    if (ncyc > 25):
        aturb = .5
    if (kturb == 1): # if kturb is one, else  
        rev0[0:ie,0:je] = rev[0:ie,0:je]
        # call turbbl
        # call turb2
        '''
        If running laminar flows, calculation is more simple.
        '''
        # runs the tubr visc calculations
        '''
        Call either the Baldwin Lomax Model or run the RNG algebraic model.
        '''
        turbulent_viscosity(params, dims)


        rev[0:ie,0:je] = aturb*rev[0:ie,0:je]  +(1.  -aturb)*rev0[0:ie,0:je]
    #  else start the rng algebraic model
    else:
    # i>1:il, i+1>2:il+1, i-1>0:il-1
    # j>1:jlm, j+1>2:jlm+1, j-1>0:jlm-1

        u[0:ie,0:je]   = w[0:ie,0:je,1]/w[0:ie,0:je,0]
        v[0:ie,0:je]   = w[0:ie,0:je,2]/w[0:ie,0:je,0]

        u[itl:itu,0]   = -u[itl:itu,1]
        v[itl:itu,0]   = -v[itl:itu,1]

        dx13      = xc[0:il,0:jl,0]   - xc[1:il+1,1:jl+1,0]
        dy13      = xc[0:il,0:jl,1]   - xc[1:il+1,1:jl+1,1]
        dx24      = xc[1:il+1,0:jl,0] - xc[0:il,1:jl+1,0]
        dy24      = xc[1:il+1,0:jl,1] - xc[0:il,1:jl+1,1]
        du13      = u[0:il,0:jl] - u[1:il+1,1:jl+1]
        dv13      = v[0:il,0:jl] - v[1:il+1,1:jl+1]
        du24      = u[1:il+1,0:jl] - u[0:il,1:jl+1]
        dv24      = v[1:il+1,0:jl] - v[0:il,1:jl+1]
        # ua        = .25*(u[i,j] + u[i+1,j+1] + u[i+1,j] + u[i,j+1])
        # va        = .25*(v[i,j] + v[i+1,j+1] + v[i+1,j] + v[i,j+1])
        dsij      = 1./(dx13*dy24 - dx24*dy13)
        dvdx      =  dsij * (dv13*dy24 - dv24*dy13)
        dudy      = -dsij * (du13*dx24 - du24*dx13)
        dudx      =  dsij * (du13*dy24 - du24*dy13)
        dvdy      = -dsij * (dv13*dx24 - dv24*dx13)
        astr[0:il,0:jl] = (dudy+dvdx)**2. +2.*(dudx**2  +dvdy**2  -((dudx+dvdy)**2)/3.)

        #   to add: 
        # call delt
    '''
    Calculates the boundary layer thickness.
    '''
    boundary_thickness(params, dims) 



    for j in range(1,jl):
        for i in range(1,il):
            xbi       = .5*(x[i-1,0,0]  +x[i,0,0])
            ybi       = .5*(x[i-1,0,1]  +x[i,0,1])
            astra     = .25*(astr[i-1,j-1]  +astr[i-1,j]+astr[i,j-1]    +astr[i,j])
            if (i>=itl) and (i<=itu+1):
                a3        = 1./(.225*abs(ynot[i]))
                ysci      = np.sqrt((xc[i,j,0]  -xbi)**2  +(xc[i,j,1]  -ybi)**2)
                ysc       = w[i,2,0]/(ysci*w[i,j,0])
                csc       = 1./(ysc+a3)**2
            else:
                csc       = (cwk*ynot[i])**2


        #     set some parameters

            rnul      = rlv[i,j]/w[i,j,0]
            rnut0     = rev[i,j]/w[i,j,0]
            # rnul3     = rnul**3
            a11       = ckr*(csc*csc*scf*scf)/rnul**2
            a2        = 75.
            a1        = a11*(astra)
            rnut0     = rnul+rnut0

            '''
            Solves for the eddy viscosity
            '''
            # DIM(X,Y) function in fortran: returns the difference X-Y if the result is positive; otherwise returns zero.
            # Replacing with: max(X-Y, 0) for python
            
            if (max(rnut0*a1-a2,0) == 0.):
                rev[i,j]  = 0.
                continue 
            else:
                rnut      = np.sqrt(a1)
            
            k      = 0
            fac    = a2 - 1.

            while k<201: 
                z=z+1
                den    = 1./(4.*rnut*rnut*rnut + fac)
                rnut1  = rnut - (rnut**4+rnut*fac  -rnut0*rnut0*a1)*den
                
                if (abs((rnut1  -rnut))<=1.e-3):
                    rev[i,j] = w[i,j,1]*max(rnut1-rnul,0)
                    break # exits the while loop, then continues
                else:
                    k      = k  +1
                    if (k>199):
                        print([' iteration not converged ',i,j])
                        print([' rnut = ',rnut,' rnut1 =',rnut1])

                        rev[i,j]  = w[i,j,0]*max(rnut1-rnul,0)
                        break # exits the while loop, then continues

                    rnut   = rnut1

    #     adjust the near wake
    ii        = ie
    for i in range(1,itl+1):
        ii        = ii  -1
        for j in range(1,jl):
            pex       = -(xc[i,1,0]  -xc[itl+1,1,0])/(20.*dsti[itl+1])
            rev[i,j]  = rev[i,j]  +(rev[itl+1,j]  -rev[i,j])*np.exp(pex)
            pex       = -(xc[ii,1,0]  -xc[itu,1,0])/(20.*dsti[itu])
            rev[ii,j] = rev[ii,j]  +(rev[itu,j]  -rev[ii,j])*np.exp(pex)

    ii        = ib  -i
            
    rev[1:il,je] = rev[1:il,jl]
    rev[1:il,0]  = rev[ii,1]

    for i in range(itl,itu):
            if (xc[i,1,0] <= xtran):
                for j in range(0,jl):
                    rev[i,j]  = 0
            rev[i,1]  = -rev[i,2]

    rev[0,0:je]  = rev[1,0:je]
    rev[ie,0:je] = rev[il,0:je]

    return
dim_var = 500
params = {
  "ie": dim_var,
  "je": dim_var,
  "kvis": 2,
  "gamma": 1,
  "rm": 1,
  "re": 1,
  "ncyc": dim_var,
  "rev": np.random.rand(dim_var+1,dim_var+1),
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
  "xc": np.random.rand(dim_var,dim_var,3),
  "ynot": np.ones(dim_var),
  "rlv": np.ones((dim_var,dim_var)),
  "dsti": np.ones(dim_var),
  "ib": 1
}
dims = {
    "il": dim_var - 1, 
    "jl": dim_var - 1,
    "ny": dim_var,
}

t0 = time.time()
compute_viscosity(params,dims)
t1 = time.time()

total = t1-t0
print(total)