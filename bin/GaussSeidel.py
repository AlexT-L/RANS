import numpy as np

class GaussSeidel():        
    # Smooth method to perform Gauss-Seidel method
    def smooth(grid, params):
        xx = np.zeros((4,4))
        rhs = np.zeros(4)
        

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
    
def bcdw(dw, itl, itu):
    '''
    Set the residuals to zero at the boundaries and 
    match the residuals across the cut in the c-mesh

    Parameters
    ----------
    dw : 3D float array
        DESCRIPTION.
    itl : TYPE
        DESCRIPTION.
    itu : TYPE
        DESCRIPTION.

    Returns
    -------
    dw : 3D float array
        boundaries that get set to 0.

    '''
    dw = border_zero(dw)
    xlen = np.shape(dw)[0]
    
    for i in range(xlen):
        if (i < itl or i > itu):
            ii  = xlen - i - 1
            dw[i,0,:] = dw[ii,1,:]
        if (i > itl or i < itl):
            dw[i,0,:] = dw[i,1,:]
    return dw

def bcrs(rs, itl, itu):
    '''
    set the modified residuals to zero at the boundaries
    and match the modified residuals across the cut in the c-mesh

    Parameters
    ----------
    rs : TYPE
        DESCRIPTION.
    itl : TYPE
        DESCRIPTION.
    itu : TYPE
        DESCRIPTION.

    Returns
    -------
    rs : TYPE
        DESCRIPTION.

    '''
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