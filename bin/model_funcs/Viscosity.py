"""This module computes viscosity coefficients

    Libraries/Modules:
    numpy\n
    BaldwinLomax\n
    BoundaryThickness\n
    """
import sys
sys.path.append('../RANS/bin')
import time

import numpy as np
from bin.model_funcs.BaldwinLomax import turbulent_viscosity
from bin.model_funcs.BoundaryThickness import boundary_thickness

def compute_viscosity(model, ws, state):
    """Computes viscosity coefficients. 
        First, computes the molecular viscosity. 
        Then continues for turbulent, Baldwin Lomax Model 
        or runs the RNG algebraic model. 
        Next, calculates the boundary layer thickness
        Solves for the eddy viscosity. 

    Attributes:
        rlv: laminar viscosity
        rev: eddy viscosity

    Notes:
        Adapted from subroutine viscf.f"""    

    # set state
    w = state
    
    # get coordinates
    x = ws.get_field('x')
    xc = ws.get_field('xc')

    # get pressure
    p = ws.get_field('p', model.className)
    
    # get eddy viscosity
    ev = ws.get_field('ev', model.className)
    lv = ws.get_field('lv', model.className)
    rlv = lv
    rev = ev

    # set geometry parameters
    pad = model.padding
    [nx, ny] = ws.field_size()
    [nxp, nyp] = [pad+nx+pad, pad+ny+pad]
    [il, jl] = [nx+1, ny+1]
    [ie, je] = [nx+2, ny+2]
    [ib, jb] = [nx+3, ny+3]

    dims = ws.get_dims()
    itl = dims['itl'] + pad
    itu = dims['itu'] + pad
    ile = int(ie/2) + pad
    chord = x[itl,pad,1] - x[ile,pad,1] # Not sure if correct FIX!!
    
    geom = ws.get_geometry()
    scal = geom['scal']
    
    # physics parameters
    kvis = model.params['kvis']
    gamma = model.params['gamma']
    rm = model.params['rm']
    re = model.params['re']
    xtran = model.params['xtran'] 
    t0 = model.params['t0']
    rmu0 = model.params['mu0']
    
    # boundary layer parameters
    ynot = np.ones(nx)
    dsti = np.ones(nx)

    
    # select model
    kturb = 1
    
    # mg_param
    mode = 1
    if ws.is_finest():
        mode = 0

    # initializing, defined later
    rev0 = np.ones((nx,ny))
    astr = np.ones((nx,ny))
    u = np.ones((nx,ny))
    v = np.ones((nx,ny))

    # useful constants
    ckr       = .0256
    cwk       = 0.
    scf       = (scal*re/chord)/(np.sqrt(gamma)*rm)
    visc_const_C1 = 1.461e-06
    visc_const_C2 = 110.3

    # compute the molecular viscosity
    tt       = p[0:ie,0:je]/w[0:ie,0:je,0]*t0
    rlv[0:ie,0:je] = visc_const_C1*tt*np.sqrt(tt)/((tt+visc_const_C2)*rmu0)


    # for laminar flows we are done.
    # for turbulent flows we are also done on the coarser grids.

    if (kvis <= 1) or (mode!=0):
        return

    aturb     = 1.
    #if (ncyc > 25): commented out bc don't know how to implement
    #    aturb = .5
    if (kturb == 1): 
        rev0[0:ie,0:je] = rev[0:ie,0:je]

        '''
        If running laminar flows, calculation is more simple.
        Call either the Baldwin Lomax Model or run the RNG algebraic model.
        '''
        turbulent_viscosity(params, dims)


        rev[0:ie,0:je] = aturb*rev[0:ie,0:je]  +(1.  -aturb)*rev0[0:ie,0:je]
    else:
        '''
        Otherwise, start the rng algebraic model.
        '''
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
        dsij      = 1./(dx13*dy24 - dx24*dy13)
        dvdx      =  dsij * (dv13*dy24 - dv24*dy13)
        dudy      = -dsij * (du13*dx24 - du24*dx13)
        dudx      =  dsij * (du13*dy24 - du24*dy13)
        dvdy      = -dsij * (dv13*dx24 - dv24*dx13)
        astr[0:il,0:jl] = (dudy+dvdx)**2. +2.*(dudx**2  +dvdy**2  -((dudx+dvdy)**2)/3.)

    '''
    Calculates the boundary layer thickness.
    '''
    boundary_thickness(model, ws, state, ynot, dsti) 



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

            rnul      = rlv[i,j]/w[i,j,0]
            rnut0     = rev[i,j]/w[i,j,0]
            a11       = ckr*(csc*csc*scf*scf)/rnul**2
            a2        = 75.
            a1        = a11*(astra)
            rnut0     = rnul+rnut0

            '''
            Solves for the eddy viscosity
            '''
            
            if (max(rnut0*a1-a2,0) == 0.):
                rev[i,j]  = 0.
                continue 
            else:
                rnut      = np.sqrt(a1)
            
            k      = 0
            fac    = a2 - 1.

            while k<201: 
                den    = 1./(4.*rnut*rnut*rnut + fac)
                rnut1  = rnut - (rnut**4+rnut*fac  -rnut0*rnut0*a1)*den
                
                if (abs((rnut1  -rnut))<=1.e-3):
                    rev[i,j] = w[i,j,1]*max(rnut1-rnul,0)
                    break 
                else:
                    k      = k  +1
                    if (k>199):
                        print([' iteration not converged ',i,j])
                        print([' rnut = ',rnut,' rnut1 =',rnut1])

                        rev[i,j]  = w[i,j,0]*max(rnut1-rnul,0)
                        break 

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
