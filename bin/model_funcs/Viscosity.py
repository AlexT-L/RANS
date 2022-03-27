"""This module computes viscosity coefficients

    Libraries/Modules:
    numpy\n
    BaldwinLomax\n
    BoundaryThickness\n
    """
import sys
sys.path.append('../../../RANS/bin')

import numpy as np
from bin.Field import pos_diff, max, min
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
    xc = ws.get_field('xc', model.className)

    # get pressure
    p = ws.get_field('p', model.className)
    
    # get eddy viscosity
    ev = ws.get_field('ev', model.className)
    lv = ws.get_field('lv', model.className)
    rlv = lv
    rev = ev

    # set geometry parameters
    PAD = model.padding
    [nx, ny] = ws.field_size()
    [nxp, nyp] = [PAD+nx+PAD, PAD+ny+PAD]
    [il, jl] = [nx+1, ny+1]
    [ie, je] = [nx+2, ny+2]
    [ib, jb] = [nx+3, ny+3]

    dims = ws.get_dims()
    itl = dims['itl'] + PAD
    itu = dims['itu'] + PAD
    ile = int(ie/2) + PAD
    
    geom = ws.get_geometry()
    scal = geom['scal']
    chord = geom['chord']
    
    # physics parameters
    kvis = model.params['kvis']
    gamma = model.params['gamma']
    rm = model.params['rm']
    re = model.params['re']
    xtran = model.params['xtran'] 
    t0 = model.params['t0']
    rmu0 = model.params['mu0']
    
    # boundary layer parameters
    ynot = np.ones(nxp)
    dsti = np.ones(nxp)

    
    # select model
    kturb = 0
    
    # mg_param
    mode = 1
    if ws.is_finest():
        mode = 0

    # initializing, defined later
    rev0 = np.ones((nxp,nyp))
    astr = np.ones((nxp,nyp))
    u = np.ones((nxp,nyp))
    v = np.ones((nxp,nyp))

    # useful constants
    ckr       = .0256
    cwk       = 0.
    scf       = (scal*re/chord)/(np.sqrt(gamma)*rm)
    visc_const_C1 = 1.461e-06
    visc_const_C2 = 110.3

    # compute the molecular viscosity
    tt       = p[1:ie+1,1:je+1]/w[1:ie+1,1:je+1,0]*t0
    rlv[1:ie+1,1:je+1] = visc_const_C1*tt*np.sqrt(tt)/((tt+visc_const_C2)*rmu0)

    # for laminar flows we are done.
    # for turbulent flows we are also done on the coarser grids.

    if (kvis <= 1) or (mode!=0):
        return

    aturb     = 1.
    #if (ncyc > 25): commented out bc don't know how to implement FIX!!!
    #    aturb = .5
    if (kturb == 1): 
        rev0[1:ie+1,1:je+1] = rev[1:ie+1,1:je+1]

        '''
        If running laminar flows, calculation is more simple.
        Call either the Baldwin Lomax Model or run the RNG algebraic model.
        '''
        turbulent_viscosity(model, ws, w) # CHECK THIS FILE FIX!!!


        rev[1:ie+1,1:je+1] = aturb*rev[1:ie+1,1:je+1]  +(1.0  -aturb)*rev0[1:ie+1,1:je+1]
    else:
        '''
        Otherwise, start the rng algebraic model.
        '''
        u[1:ie+1,1:je+1]   = w[1:ie+1,1:je+1,1]/w[1:ie+1,1:je+1,0]
        v[1:ie+1,1:je+1]   = w[1:ie+1,1:je+1,2]/w[1:ie+1,1:je+1,0]

        u[itl:itu,0]   = -u[itl:itu,1]
        v[itl:itu,0]   = -v[itl:itu,1]

        dx13      = xc[1:il+1,1:jl+1,0]   - xc[2:il+2,2:jl+2,0]
        dy13      = xc[1:il+1,1:jl+1,1]   - xc[2:il+2,2:jl+2,1]
        dx24      = xc[2:il+2,1:jl+1,0] - xc[1:il+1,2:jl+2,0]
        dy24      = xc[2:il+2,1:jl+1,1] - xc[1:il+1,2:jl+2,1]
        du13      = u[1:il+1,1:jl+1] - u[2:il+2,2:jl+2]
        dv13      = v[1:il+1,1:jl+1] - v[2:il+2,2:jl+2]
        du24      = u[2:il+2,1:jl+1] - u[1:il+1,2:jl+2]
        dv24      = v[2:il+2,1:jl+1] - v[1:il+1,2:jl+2]
        dsij      = 1.0/(dx13*dy24 - dx24*dy13)
        dvdx      =  dsij * (dv13*dy24 - dv24*dy13)
        dudy      = -dsij * (du13*dx24 - du24*dx13)
        dudx      =  dsij * (du13*dy24 - du24*dy13)
        dvdy      = -dsij * (dv13*dx24 - dv24*dx13)
        astr[1:il+1,1:jl+1] = (dudy+dvdx)**2 +2.0*(dudx**2  +dvdy**2  -((dudx+dvdy)**2)/3.0)

    '''
    Calculates the boundary layer thickness.
    '''
    boundary_thickness(model, ws, state, ynot, dsti)



    for j in range(PAD,ny+PAD):
        for i in range(PAD,nx+PAD):
            xbi       = 0.5*(x[i-PAD,0,0]  +x[i+1-PAD,0,0])
            ybi       = 0.5*(x[i-PAD,0,1]  +x[i+1-PAD,0,1])
            astra     = 0.25*(astr[i-1,j-1]  +astr[i-1,j]+astr[i,j-1]    +astr[i,j])
            if (i>=itl-1) and (i<=itu):
                a3        = 1.0/(0.225*abs(ynot[i]))
                ysci      = np.sqrt((xc[i,j,0]  -xbi)**2  +(xc[i,j,1]  -ybi)**2)
                ysc       = w[i,2,0]/(ysci*w[i,j,0])
                csc       = 1.0/(ysc+a3)**2
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
            
            if (pos_diff(rnut0*a1,a2) == 0.):
                rev[i,j]  = 0.
                continue 
            else:
                rnut      = np.sqrt(a1)
            
            fac    = a2 - 1.

            MAX_ITER=200
            for k in range(MAX_ITER):
                den    = 1.0/(4.0*rnut**3 + fac)
                rnut1  = rnut - (rnut**4+rnut*fac  -rnut0*rnut0*a1)*den
                
                if (abs((rnut1  -rnut))<=1.0e-3):
                    rev[i,j] = w[i,j,0]*pos_diff(rnut1,rnul)
                    break 
                else:
                    if (k==MAX_ITER-1):
                        print(" iteration not converged ("+str(i)+", "+str(j)+")")
                        print(" rnut = "+str(rnut)+" rnut1 = "+str(rnut1))
                        
                        rev[i,j]  = w[i,j,0]*pos_diff(rnut1,rnul)
                        break 

                    rnut   = rnut1

    #     adjust the near wake
    ii        = il+1
    for i in range(PAD,itl+1):
        ii        = ii  -1
        for j in range(PAD,ny+PAD):
            pex       = -(xc[i,2,0]  -xc[itl,2,0])/(20.0*dsti[itl])
            rev[i,j]  = rev[i,j]  +(rev[itl,j]  -rev[i,j])*np.exp(pex)
            pex       = -(xc[ii,2,0]  -xc[itu-1,2,0])/(20.0*dsti[itu-1])
            rev[ii,j] = rev[ii,j]  +(rev[itu-1,j]  -rev[ii,j])*np.exp(pex)

    for i in range(PAD,nx+PAD):
        ii           = ib  -i 
        rev[i,je] = rev[i,jl]
        rev[1,1]  = rev[ii,2]

    for i in range(itl,itu):
            if (xc[i,2,0] <= xtran):
                for j in range(1,jl+1):
                    rev[i,j]  = 0
            rev[i,1]  = -rev[i,2]

    rev[1 ,1:je+1] = rev[2 ,1:je+1]
    rev[ie,1:je+1] = rev[il,1:je+1]

    return
