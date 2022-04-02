import numpy as np
from bin.Field import Field, copy, isfinite, maximum, minimum, abs, pos_diff, sqrt, square, min

def stability(self, model, workspace, state):
    """Calculates timestep limits to maintain stability
    
    Args:
        model (Model):  The physics model
        workspace (Workspace): The current Workspace
        state (Field): Field containing current state
    """

    # padding
    PAD = self.padding
    
    # set stability parameters
    slim      = 0.001
    rlim      = 0.001
    b         = 4.0
    dtmin     = 0.0
    imin      = 0
    jmin      = 0
    iprec     = 0 # turns on gauss-seidel preconditioner

    # get relevant geometry parameters
    dims = workspace.get_dims()
    [nx, ny] = workspace.field_size()
    [nxp, nyp] = [PAD+nx+PAD, PAD+ny+PAD]
    ip = PAD
    jp = PAD
    il = nx+1
    jl = ny+1
    ie = nx+2
    je = ny+2
    ib = nx+3
    jb = ny+3
    itl = dims['itl'] + PAD
    itu = dims['itu'] + PAD

    # flo_param
    gamma = model.params['gamma']
    rm = model.params['rm']
    re = model.params['re']
    prn = model.params['prn']
    prt = model.params['prt']
    adis = model.params['adis']
    cfl = model.cfl
    kvis = model.params['kvis']

    # retrieve working arrays from model
    def get(varName):
            return workspace.get_field(varName, model.className)
    ev = get('ev')
    lv = get('lv')
    radI = get("radI")
    radJ = get("radJ")
    rfli = get("rfli")
    rflj = get("rflj")
    rfl = get("rfl")
    dtl = get("dtl")
    w = state
    p = get("p")
    vol = get("vol")
        
    
    # local working arrays
    def get_local(varName):
            return workspace.get_field(varName, self.className)
    s = get_local("s")
    dtlc = get_local("dtlc")

    # dim helper function
    def dim(a, b):
        diff = a-b
        return maximum(diff, diff*0)

    # c
    # c     permissible time step
    # c
    RLIM = Field.create((nx,ny), rlim)

    cc = gamma*p[ip:ie,jp:je]/maximum(w[ip:ie,jp:je,0],RLIM)

    i_edges = workspace.edges(0)
    j_edges = workspace.edges(1)

    dj = (i_edges[0:nx,:] + i_edges[1:nx+1,:])/2
    di = (j_edges[:,0:ny] + j_edges[:,1:ny+1])/2

    dix = di[:,:,0]
    diy = di[:,:,1]
    djx = dj[:,:,0]
    djy = dj[:,:,1]

    qs                  = (djy*w[ip:ie,jp:je,1] - djx*w[ip:ie,jp:je,2])/w[ip:ie,jp:je,0]
    cs                  = cc*(square(djx) + square(djy))
    radI[ip:ie,jp:je]   = abs(qs) + sqrt(cs)
    qs                  = (dix*w[ip:ie,jp:je,2] - diy*w[ip:ie,jp:je,1])/w[ip:ie,jp:je,0]
    cs                  = cc*(square(dix) + square(diy))
    radJ[ip:ie,jp:je]   = abs(qs) + sqrt(cs)
    dtl[ip:ie,jp:je]    = 1/(radI[ip:ie,jp:je] + radJ[ip:ie,jp:je])
    dtlc[ip:ie,jp:je]   = radI[ip:ie,jp:je] + radJ[ip:ie,jp:je]

    # c
    # c     pressure or entropy switch
    # c
    s[:] = copy(p)

    # c
    # c     adaptive time step
    # c
    CFL = Field.create((nx, ny), pos_diff(cfl, 1))

    dpi = abs((s[ip+1:ib, jp:je] - 2*s[ip:ie, jp:je] + s[1:il, jp:je]) / \
              (s[ip+1:ib, jp:je] + 2*s[ip:ie, jp:je] + s[1:il, jp:je] + slim))
    dpj = abs((s[ip:ie, jp+1:jb] - 2*s[ip:ie, jp:je] + s[ip:ie, 1:jl]) / \
              (s[ip:ie, jp+1:jb] + 2*s[ip:ie, jp:je] + s[ip:ie, 1:jl] + slim))
    
    rfl[ip:ie, jp:je] = 1/(1 + minimum(CFL, b*(dpi+dpj)))

    # c
    # c     fixed time step
    # c
    if not self.local_timestepping:
        dtmin = min(dtl)

        rfl[ip:ie, jp:je] = dtmin/dtl[ip:ie, jp:je]

    # c
    # c     option to rescale the dissipative coefficients
    # c
    if iprec != 0:
        rfli[ip:ie, jp:je] = (radJ[ip:ie, jp:je]/radI[ip:ie, jp:je])**(0.25)
        rflj[ip:ie, jp:je] = (radI[ip:ie, jp:je]/radJ[ip:ie, jp:je])**(0.25)

    a                   = (radI[ip:ie, jp:je]/radJ[ip:ie, jp:je])**adis
    radI[ip:ie, jp:je]  = radI[ip:ie, jp:je]*(1 + 1.0/a)
    radJ[ip:ie, jp:je]  = radJ[ip:ie, jp:je]*(1 + a)

    if iprec == 0:
        rfli[ip:ie, jp:je] = radI[ip:ie, jp:je]/dtlc[ip:ie, jp:je]
        rflj[ip:ie, jp:je] = radJ[ip:ie, jp:je]/dtlc[ip:ie, jp:je]
    
    print("iprec: "+str(iprec))
    # c
    # c
    # c     reduce the artificial dissipation for viscous flows
    # c     and adjust the time step estimate
    # c
    if (kvis >= 0):

        v1        = np.sqrt(gamma).item()*rm/re
        v2        = 0.0
        if (kvis > 1): 
            v2 = v1

        k = gamma *(lv[ip:ie, jp:je] * v1/prn + ev[ip:ie, jp:je] * v2/prt)/w[ip:ie, jp:je,0]
        mu = (lv[ip:ie, jp:je]*v1 + ev[ip:ie, jp:je]*v2)/w[ip:ie, jp:je,0]

        dsi = square(djx) + square(djy)
        dsj = square(dix) + square(diy)
        vsi = (k*dsi + mu*sqrt(dsi*dsj)/6)/vol[ip:ie, jp:je]
        vsj = (k*dsj + mu*sqrt(dsi*dsj)/6)/vol[ip:ie, jp:je]
        dtv = dtlc[ip:ie, jp:je] + 4*(vsi+vsj)
        
        dtl[ip:ie, jp:je] = 1/dtv
        radI[ip:ie, jp:je] = pos_diff(radI[ip:ie, jp:je], vsi)
        radJ[ip:ie, jp:je] = pos_diff(radJ[ip:ie, jp:je], vsj)

    # c
    # c     set boundary values at i=1 and i=ie
    # c
    radI[1, jp:je]  = radI[2, jp:je]
    radI[ie, jp:je] = radI[il, jp:je]
    rfl[1, jp:je]   = rfl[2, jp:je]
    rfl[ie, jp:je]  = rfl[il, jp:je]
    rfli[1, jp:je]  = rfli[2, jp:je]
    rfli[ie, jp:je] = rfli[il, jp:je]
    rflj[1, jp:je]  = rflj[2, jp:je]
    rflj[ie, jp:je] = rflj[il, jp:je]
    dtl[1, jp:je]   = dtl[2, jp:je]
    dtl[ie, jp:je]  = dtl[il, jp:je]

    # c
    # c     set boundary values at j=1 and j=je
    # c
    radJ[1:ib,1]   = radJ[1:ib,2]
    radJ[1:ib,je]  = radJ[1:ib,jl]
    rfl[1:ib,1]    = rfl[1:ib,2]
    rfl[1:ib,je]   = rfl[1:ib,jl]
    rfli[1:ib,1]   = rfli[1:ib,2]
    rfli[1:ib,je]  = rfli[1:ib,jl]
    rflj[1:ib,1]   = rflj[1:ib,2]
    rflj[1:ib,je]  = rflj[1:ib,jl]
    dtl[1:ib,1]    = dtl[1:ib,2]
    dtl[1:ib,je]   = dtl[1:ib,jl]

    # c
    # c     set boundary values along the cut
    # c
    radJ[ie:itu-1:-1,1]  = radJ[1:itl, 2]
    radJ[1:itl, 1]   = radJ[ie:itu-1:-1,2]
    rfl[ie:itu-1:-1,1]   = rfl[1:itl, 2]
    rfl[1:itl, 1]    = rfl[ie:itu-1:-1,2]
    rfli[ie:itu-1:-1,1]  = rfli[1:itl, 2]
    rfli[1:itl, 1]   = rfli[ie:itu-1:-1,2]
    rflj[ie:itu-1:-1,1]  = rflj[1:itl, 2]
    rflj[1:itl, 1]   = rflj[ie:itu-1:-1,2]
    dtl[ie:itu-1:-1,1]   = dtl[1:itl, 2]
    dtl[1:itl, 1]    = dtl[ie:itu-1:-1,2]

    return dtlc