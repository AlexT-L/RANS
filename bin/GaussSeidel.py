import numpy as np

class GaussSeidel():        
    # Smooth method to perform Gauss-Seidel method
    def smooth(grid, params):
        xx = np.zeros((4,4))
        rhs = np.zeros(4)
        

def rhs_face():
    '''
    Calculate the face contributions to the right 
    hand side and preconiditioning matrix

    Returns
    -------
    None.

    '''
    xx = np.zeros((4,4))
    
    
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