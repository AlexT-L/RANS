import numpy as np

from bin.Field import mismatch_mul

def nsflux(self,ws,w,dw,rfil):
    
    # calculate viscous fluxes given a workspace

    # get geometry dictionary
    geom = ws.get_geometry()
    
    # dims
    PAD = 2
    [nx, ny] = ws.field_size()
    [il, jl] = [nx+1, ny+1]
    [ie, je] = [nx+2, ny+2]
    [ib, jb] = [nx+3, ny+3]
    
    # geometric parameters
    geom = ws.get_geometry()
    scal = geom['scal']
    chord = geom['chord']

    # flow related variabless
    def get(varName):
        return ws.get_field(varName, self.className)
    p = get('p') # pressure
    lv = get('lv') # laminar viscocity
    ev = get('ev') # eddy viscocity
    vw = get('vw') # storage for viscous residuals

    # flow parameters
    gamma = self.params['gamma']
    mach = self.params['rm']
    Re = self.params['re']
    prn = self.params['prn']
    prt = self.params['prt']
    
    # mg_param
    mode = 1
    if ws.is_finest():
        mode = 0

    # mesh related vars
    x = ws.get_field('x') # mesh vertices
    xc = get('xc') # mesh centers
    
    # get flux and velocity working vectors
    q = ws.get_field('q', self.className)
    u = ws.get_field('u', self.className)
    fs = ws.get_field('fs', self.className)

    # set up parameters
    sfil      = 1.  -rfil
    gm1       = gamma - 1.
    scf       = rfil*np.sqrt(gamma)*mach/(scal*Re/chord)
    
    # c
    # c     calcluate the velocities and temperature
    # c
    u[1:ib,1:jb,0] = w[1:ib,1:jb,1]/w[1:ib,1:jb,0]
    u[1:ib,1:jb,1] = w[1:ib,1:jb,2]/w[1:ib,1:jb,0]
    u[1:ib,1:jb,2] = p[1:ib,1:jb]/(gm1*w[1:ib,1:jb,0])
    
# c
# c     evaluate the viscous terms
# c

    # c
    # c     evaluate the velocity and temperature gradients
    # c
    dxi = xc[2:ib,2:jb,0] - xc[1:ie,2:jb,0] + xc[2:ib,1:je,0] - xc[1:ie,1:je,0]
    dxj = xc[2:ib,2:jb,0] + xc[1:ie,2:jb,0] - xc[2:ib,1:je,0] - xc[1:ie,1:je,0]
    dyi = xc[2:ib,2:jb,1] - xc[1:ie,2:jb,1] + xc[2:ib,1:je,1] - xc[1:ie,1:je,1]
    dyj = xc[2:ib,2:jb,1] + xc[1:ie,2:jb,1] - xc[2:ib,1:je,1] - xc[1:ie,1:je,1]
    dsj = 1.0/(dxi*dyj - dyi*dxj)
    
    for l in range(3):
        dui = u[2:ib,2:jb,l] - u[1:ie,2:jb,l] + u[2:ib,1:je,l] - u[1:ie,1:je,l]
        duj = u[2:ib,2:jb,l] + u[1:ie,2:jb,l] - u[2:ib,1:je,l] - u[1:ie,1:je,l]
        q[2:ib,2:jb,l,0] = (dui*dyj  -duj*dyi)*dsj
        q[2:ib,2:jb,l,1] = (duj*dxi  -dui*dxj)*dsj
        
    # c
    # c     evaluate the viscous stress tensor and dissipation
    # c
    ua        = u[2:ib,2:jb,0]  +u[1:ie,2:jb,0]  +u[2:ib,1:je,0]  +u[1:ie,1:je,0]
    va        = u[2:ib,2:jb,1]  +u[1:ie,2:jb,1]  +u[2:ib,1:je,1]  +u[1:ie,1:je,1]
    eva       = ev[2:ib,2:jb]  +ev[1:ie,2:jb]  +ev[2:ib,1:je]  +ev[1:ie,1:je]
    lva       = lv[2:ib,2:jb]  +lv[1:ie,2:jb]  +lv[2:ib,1:je]  +lv[1:ie,1:je]
    rmu       = 0.25*(eva+lva)
    rlam      = 2.0*rmu/3.0
    rk        = gamma*(lva/prn  +eva/prt)
    dq        = rlam*(q[1:ie,1:je,0,0] + q[1:ie,1:je,1,1])
    
    q[2:ib,2:jb,0,0]  = rmu*(q[2:ib,2:jb,0,0]  +q[2:ib,2:jb,0,0])  -dq
    q[2:ib,2:jb,1,1]  = rmu*(q[2:ib,2:jb,1,1]  +q[2:ib,2:jb,1,1])  -dq
    q[2:ib,2:jb,0,1]  = rmu*(q[2:ib,2:jb,0,1]  +q[2:ib,2:jb,1,0])
    q[2:ib,2:jb,1,0]  = q[2:ib,2:jb,0,1]
    q[2:ib,2:jb,2,0]  = 0.15*(rk*q[2:ib,2:jb,2,0]  +ua*q[2:ib,2:jb,0,0]  +va*q[2:ib,2:jb,1,0])
    q[2:ib,2:jb,2,1]  = 0.15*(rk*q[2:ib,2:jb,2,1]  +ua*q[2:ib,2:jb,0,1]  +va*q[2:ib,2:jb,1,1])

# c
# c     viscous fluxes in the i direction
# c
    dx = x[0:il,1:jl,0] - x[0:il,0:jl-1,0]
    dy = x[0:il,1:jl,1] - x[0:il,0:jl-1,1]
    fs[1:ie,2:je,:] = mismatch_mul(dy, q[1:ie,2:je,:,0] + q[1:ie,1:jl,:,0]) - mismatch_mul(dx, q[1:ie,2:je,:,1] + q[1:ie,1:jl,:,1])

    vw[PAD:nx+PAD,PAD:ny+PAD,1:4] = sfil*vw[PAD:nx+PAD,PAD:ny+PAD,1:4] - scf*(fs[2:ie,2:je,:] - fs[1:il,2:je,:])

# c
# c     viscous fluxes in the j direction
# c
    dx = x[1:il,0:jl,0] - x[0:il-1,0:jl,0]
    dy = x[1:il,0:jl,1] - x[0:il-1,0:jl,1]
    fs[2:ie,1:je,:] = mismatch_mul(dx, q[2:ie,1:je,:,1] + q[2:ie,1:je,:,1]) - mismatch_mul(dy, q[2:ie,1:je,:,0] + q[2:ie,1:je,:,0])
    
    vw[PAD:nx+PAD,PAD:ny+PAD,1:4] -= scf*(fs[2:ie,2:je,:] - fs[2:ie,1:jl,:])


    # add to flux field
    dw += vw
    
    return