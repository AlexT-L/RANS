import numpy as np
from numpy.linalg import inv

class GaussSeidel():        
    # Smooth method to perform Gauss-Seidel method
    def smooth(grid, params):
        '''
        I think there should be a neat way to combine main loop into 4 function
        calls, 1 for each left, right, top, and bottom face. But will require 
        some thinking about the correct indexing on things. 
        '''
        
        # Maybe these are actually params
        ie = params['ie'] # Mesh dimension
        je = params['je'] # Mesh dimension
        kvis = params['kvis']
        gamma = params['gamma']
        rm = params['rm']
        re = params['re']
        nstage = params['nstage']
        il = params['il']
        jl = params['jl']
        eps = params['eps']
        cfl = params['cfl']
        dtl = params['dtl']
        
        # Very incorrect ones I think
        dw = params['dw']
        rs = params['rs']
        x = params['x'] 
        w = params['w']
        p = params['p']
        
        # Used in inner functions so this has to get here somehow
        vol = params['vol']
        rlv = params['rlv']
        rev = params['rev']
        
        
        
        r00 = np.zeros((ie, je))
        u00 = np.zeros((ie, je))
        v00 = np.zeros((ie, je))
        p00 = np.zeros((ie, je))
        
        sgrmrei = 0
        if (kvis > 0):
            sgrmrei = np.sqrt(gamma)*rm/re
        
        niter = 2
        
        # save the velocities, density, and pressure
        
        if (nstage == 1):
            for i in range(ie):
                for j in range(je):
                    r00[i,j] = w[i,j,0]
                    u00[i,j] = w[i,j,1]/w[i,j,0]
                    v00[i,j] = w[i,j,2]/w[i,j,0]
                    p00[i,j] = p[i,j]
                    
        ### Setup and solve linear implicit system
        
        # Transform the residuals to the symmetrizing variables
        dw = res_cons_to_sym(params, dw)
        
        # set the residuals to zero at the boundaries
        dw = bcdw(params, dw)
        
        # set the modified residuals to zero
        for i in range(ie):
            for j in range(je):
                rs[i,j,:] = np.zeros(4)
        
        ### Solution by symmetric gauss-seidel iteration
        idir = -1
        
        for i in range(niter):
            idir = -idir
            jdir = idir
            
            rs = bcrs(params, rs)
            
            if (idir == 1):
                ia = 2
                ja = 2
                ibf = il
                jbf = jl
            else:
                ia = il
                ja = jl
                ibf = 2
                jbf = 2
            
            for ja in range(ja, jbf, jdir):
                for i in range(ia, ibf, idir):
                    xx = np.zeros((4,4))
                    
                    dtv2 = eps*cfl*dtl[i,j]
                    
                    di = -1
                    dj = 0
                    sxa = x[i-1,j-1,1] - x[i-1,j,1]
                    sya = x[i-1,j,0] - x[i-1,j-1,0]
                    (xx, rhs_left) = interface(params, i, j, di, dj, sxa, sya,
                                               dtv2, r00, u00, v00, p00,
                                               sgrmrei, rs, xx)
                    
                    di = 1
                    dj = 0
                    sxa = x[i,j,1] - x[i,j-1,1]
                    sya = x[i,j-1,0] - x[i,j,0]
                    (xx, rhs_right) = interface(params, i, j, di, dj, sxa, sya,
                                                dtv2, r00, u00, v00, p00,
                                               sgrmrei, rs, xx)
                    
                    di = 0
                    dj = -1
                    sxa = x[i,j-1,1] - x[i-1,j-1,1]
                    sya = x[i-1,j-1,0] - x[i,j-1,0]
                    (xx, rhs_bot) = interface(params, i, j, di, dj, sxa, sya,
                                              dtv2, r00, u00, v00, p00,
                                               sgrmrei, rs, xx)
                    
                    di = 0
                    dj = 1
                    sxa = x[i-1,j,1] - x[i,j,1]
                    sya = x[i,j,0] - x[i-1,j,0]
                    (xx, rhs_top) = interface(params, i, j, di, dj, sxa, sya,
                                              dtv2, r00, u00, v00, p00,
                                               sgrmrei, rs, xx)
                    
                    # Assemble complete right hand side
                    rhs = np.zeros(len(rhs_left))
                    for k in range(len(rhs_left)):
                        rhs[i] = dw[i,j,k] - rhs_left[k] - rhs_right[k] \
                                - rhs_bot[k] - rhs_top[k]
                    
                    for k in range(4):
                        xx[k,k] += 1
                        
                    
                    binv = inv(xx)
                    x = np.matmul(binv, rhs)
                    
                    # Save x in residual?
                    rs[i,j,:] = x
    
        # Not sure what to do with rs now
        return rs
        
def interface(params, i, j, di, dj, sxa, sya, dtv2, r00, u00, v00, p00, 
              sgrmrei, rs, xx):
    '''
    Does a bunch of interfaces. The important parameters are the i and j 
    indices along with the di and dj which tell us which direction we are 
    looking in. One of di or dj will be 0, the other +-1

    Returns
    -------
    xx : 4x4 float array
    
    rhs_left : float vector length 4

    '''
    vol  = params['vol']
    rlv  = params['rlv']
    rev  = params['rev']
    kvis = params['kvis']
    
    
    # Left variables
    rl = r00[i,j]
    pl = p00[i,j]
    ul = u00[i,j]
    vl = v00[i,j]
    
    # Left i interface
    a_vec = rs[i+di,j+dj,:]
    sar2 = sxa**2 + sya**2
    sar = np.sqrt(sar2)
    vnx = sxa/sar
    vny = sya/sar
    
    # Right variables
    rr = r00[i+di,j+dj]
    pr = p00[i+di,j+dj]
    ur = u00[i+di,j+dj]
    vr = v00[i+di,j+dj]
    
    # Viscous Coefficients
    rhoa = 0.5*(rl + rr)
    svol = 0.5*(vol[i,j] + vol[i+di,j+dj])
    scale = sgrmrei*sar2/(svol*rhoa)
    if (kvis == 0):
        scale = 0
    rmuel = 0.5*scale*(rlv[i,j] + rlv[i+di,j+dj])
    rmuet = 0.5*scale*(rev[i,j] + rev[i+di,j+dj])
    
    (xx, rhs_left) = rhs_face(params, vnx, vny, sar, dtv2, rmuel, 
                         rmuet, rl, pl, ul, vl, 
                         rr, pr, ur, vr, a_vec, xx)
    return (xx, rhs_left)

def rhs_face(params, vnx, vny, sar, dtv2, rmuel, rmuet, rl, pl, ul, vl, rr, pr,
             ur, vr, a_vec, xx):
    '''
    Calculate the face contributions to the right 
    hand side and preconiditioning matrix

    Returns
    -------
    None.

    '''
    # Unpack params
    gamma = params['gamma']
    kvis  = params['kvis']
    ncyc  = params['ncyc']
    prn   = params['prn']
    prt   = params['prt']
    diag  = params['diag']
    
    
    gm1 = gamma - 1
    dlim = 0.10
    if (kvis == 0):
        dlim = 0.25
    if (ncyc < 5):
        dlim = 0.5
        
    # Arithmetic averaging
    rhoa = 0.5*(rl + rr)
    ua = 0.5*(ul + ur)
    va = 0.5*(vl + vr)
    pa = 0.5*(pl + pr)
    qn = vnx*ua + vny*va
    cc = gamma*pa/rhoa
    c = np.sqrt(cc)
    
    # Viscous Coefficients
    epsv1 = 1.0
    epsv2 = 0.25
    amu   = epsv1*dtv2*(rmuel + rmuet)
    amupr = epsv2*dtv2*(rmuel/prn + rmuet/prt)
    
    # Eigenvalues
    e1 = abs(qn)
    e2 = abs(qn + c)
    e3 = abs(qn - c)
    a = dlim*c
    if (e1 < a):
        e1 = 0.5*(a + e1**2/a)
    if (e2 < a):
        e2 = 0.5*(a + e2**2/a)
    if (e3 < a):
        e3 = 0.5*(a + e3**2/a)
    a = 0.5*sar*dtv2
    aqsa = a*diag*abs(qn)
    
    # Form the right hand side with negative Jacobian
    q1 = a*(qn - e1)
    q2 = a*(qn + c - e2)
    q3 = a*(qn - c - e3)
    r1 = 0.5*(q2 + q3) - q1
    r2 = 0.5*(q2 - q3)
    
    vq1 = q1 - amu
    vr1 = r1 - amu/3.0
    
    a11 = r1 + q1 - gm1*amupr - aqsa
    a12 = vnx*r2
    a13 = vny*r2
    a14 = -amupr
    a22 = vnx**2*vr1 + vq1 - aqsa
    a23 = vnx*vny*vr1
    a33 = vny**2*vr1 + vq1 - aqsa
    a41 = -gm1*amupr
    a44 = q1 - amupr - aqsa
    
    rhs = np.zeros(4)
    a1 = a[0]
    a2 = a[1]
    a3 = a[2]
    a4 = a[3]
    
    rhs[0]   = a11*a1  +a12*a2  +a13*a3  +a14*a4
    rhs[1]   = a12*a1  +a22*a2  +a23*a3
    rhs[2]   = a13*a1  +a23*a2  +a33*a3
    rhs[3]   = a41*a1  +a44*a4
    
    # Form the preconditioning matrix
    xx[1,1] -= a11
    xx[1,2] -= a12
    xx[1,3] -= a13
    xx[1,4] -= a14
    
    xx[2,1] -= a12
    xx[2,2] -= a22
    xx[2,3] -= a23
    
    xx[3,1] -= a13
    xx[3,2] -= a23
    xx[3,3] -= a33
    
    xx[4,1] -= a14
    xx[4,4] -= a44
    
    return xx, rhs
    
def bcdw(params, dw):
    '''
    Set the residuals to zero at the boundaries and 
    match the residuals across the cut in the c-mesh

    Parameters
    ----------
    dw : 3D float array

    Returns
    -------
    dw : 3D float array
        boundaries that get set to 0.

    '''
    itl = params['itl']
    itu = params['itu']
    
    dw = border_zero(dw)
    xlen = np.shape(dw)[0]
    
    for i in range(xlen):
        if (i < itl or i > itu):
            ii  = xlen - i - 1
            dw[i,0,:] = dw[ii,1,:]
        if (i > itl or i < itl):
            dw[i,0,:] = dw[i,1,:]
    return dw

def bcrs(params, rs):
    '''
    set the modified residuals to zero at the boundaries
    and match the modified residuals across the cut in the c-mesh
    '''
    itl = params['itl']
    itu = params['itu']
    
    rs = border_zero(rs)
    xlen = np.shape(rs)[0]
    
    for i in range(xlen):
        if (i <= itl or i > itu):
            ii = xlen - i - 1
            rs[i,0,:] = rs[ii,2,:]
    
    return rs

def border_zero(dw):
    '''
    Sets the edges in the first 2 dimensions of dw to 0

    Parameters
    ----------
    dw : 3D array of floats

    Returns
    -------
    dw : 3D array of floats
        First and second indices of 0 and -1 have values 
        set to 0.

    '''
    
    # Length of indices
    xlen = np.shape(dw)[0]
    ylen = np.shape(dw)[1]
    zlen = np.shape(dw)[2]
    
    # Set upper and lower boundary to 0
    dw[:,0,:] = np.zeros((xlen, zlen,))
    dw[:,-1,:] = np.zeros((xlen, zlen))
    
    # Set left adn right boundary to 0
    dw[0,:,:] = np.zeros((ylen, zlen))
    dw[-1,:,:] = np.zeros((ylen, zlen))
    
    return dw

def res_cons_to_sym(params, dw):
    '''
    conservative residuals dw changed to nonconservative residuals

    Parameters
    ----------
    dw : 3D float array 

    Returns
    -------
    dw : 3D float array

    '''
    ilower = params['ilower']
    jlower = params['jlower']
    w = params['w']
    gm1 = params['gm1']
    gamma = params['gamma']
    p = params['p']
    
    for i in range(2, ilower):
        for j in range(2, jlower):
            # Unpack dw for ease
            r0 = dw[i,j,0]
            r1 = dw[i,j,1]
            r2 = dw[i,j,2]
            r3 = dw[i,j,3]
            
            ua = w[i,j,1] / w[i,j,0]
            va = w[i,j,2] / w[i,j,0]
            qq = 0.5*(ua**2 + va**2)
            cc = gamma*p[i,j] / w[i,j,0]
            c  = np.sqrt(cc)
            
            dw[i,j,0] = gm1*(r3 + qq*r0 - ua*r1 - va*r2) / cc
            dw[i,j,1] = (r1 - ua*r0) / c
            dw[i,j,2] = (r2 - va*r0) / c
            dw[i,j,3] = dw[i,j,0] - r0
    
    return dw
    
def res_sym_to_cons(params, rs, dw):
    '''
    nonconservative residuals dw changed to conservative residuals

    Parameters
    ----------
    params : Dict of parameters
    rs : 3D float array
    dw : 3D float array

    Returns
    -------
    dw : 3D float array

    '''
    
    ilower = params['ilower']
    jlower = params['jlower']
    w = params['w']
    p = params['p']
    gamma = params['gamma']
    
    
    for j in range(2, jlower):
        for i in range(2, ilower):
            # Unpack
            r0 = rs[i,j,0]
            r1 = rs[i,j,1]
            r2 = rs[i,j,2]
            r3 = rs[i,j,3]
            
            ua = w[i,j,1] / w[i,j,0]
            va = w[i,j,2] / w[i,j,0]
            ha = (w[i,j,3] + p[i,j]) / w[i,j,0]
            qq = 0.5*(ua**2 + va**2)
            cc = gamma*p[i,j] / w[i,j,0]
            c = np.sqrt(cc)
            
            dw[i,j,0] = r0 - r3
            dw[i,j,1] = ua*(r0-r3) + c*r1
            dw[i,j,2] = va*(r0-r3) + c*r2
            dw[i,j,3] = ha*r0 + c*(ua*r1 + va*r2) - qq*r3
    
    
    return dw