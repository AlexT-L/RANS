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

def res_cons_to_sym(dw):
    '''
    conservative residuals dw changed to nonconservative residuals

    Parameters
    ----------
    dw : TYPE
        DESCRIPTION.

    Returns
    -------
    dw : TYPE
        DESCRIPTION.

    '''
    
    for i in range(ilower):
        for j in range(jlower):
            # Do somethings
    
    
    return dw
    
    
    
    