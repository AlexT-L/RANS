import numpy as np

def stability(self, model, workspace, state):
    """Calculates timestep limits to maintain stability
    
    Args:
        model (Model):  The physics model
        workspace (Workspace): The current Workspace
        state (Field): Field containing current state
    """


    # padding
    pad = self.padding
    
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
    grid = workspace.get_grid()
    [nx, ny] = workspace.field_size()
    il = nx+1
    jl = ny+1
    ie = nx+2
    je = ny+2
    ib = nx+3
    jb = nx+3
    itl = dims['itl']
    itu = dims['itu']

    # flo_param
    gamma = model.params['gamma']
    rm = model.params['rm']
    rho0 = model.params['rho0']
    p0 = model.params['p0']
    h0 = model.params['h0']
    c0 = model.params['c0']
    u0 = model.params['u0']
    v0 = model.params['v0']
    ca = model.params['ca']
    sa = model.params['sa']
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
        return max(diff, 0)

    # edge function decorator
    def edge(i, j, side):
        pad = model.padding
        return workspace.edge(i-pad, j-pad, side)

    # c
    # c     permissible time step
    # c
    for j in range(pad,ny+pad):
        for i in range(pad,nx+pad):
            cc        = gamma*p[i,j]/np.maximum(w[i,j,0],rlim)

            # get side lengths
            [xn, yn] = edge(i,j,'n')
            [xs, ys] = edge(i,j,'s')
            [xe, ye] = edge(i,j,'e')
            [xw, yw] = edge(i,j,'w')

            xx = np.mean([xn, xs])
            yx = np.mean([yn, ys])
            xy = np.mean([xe, xw])
            yy = np.mean([ye, yw])

            # xx        = 0.5*(x(i,j,1) -x(i-1,j,1) +x(i,j-1,1) -x(i-1,j-1,1))
            # yx        = 0.5*(x(i,j,2) -x(i-1,j,2) +x(i,j-1,2) -x(i-1,j-1,2))
            # xy        = 0.5*(x(i,j,1) -x(i,j-1,1) +x(i-1,j,1) -x(i-1,j-1,1))
            # yy        = 0.5*(x(i,j,2) -x(i,j-1,2) +x(i-1,j,2) -x(i-1,j-1,2))
            qs        = (yy*w[i,j,1] -xy*w[i,j,2])/w[i,j,0]
            cs        = cc*(xy**2  +yy**2)
            radI[i,j] = abs(qs)  +np.sqrt(cs)
            qs        = (xx*w[i,j,2]  -yx*w[i,j,1])/w[i,j,0]
            cs        = cc*(xx**2  +yx**2)
            radJ[i,j] = abs(qs)  +np.sqrt(cs)
            dtl[i,j]  = 1.0/(radI[i,j]  +radJ[i,j])
            dtlc[i,j] = radI[i,j]  +radJ[i,j]

    # c
    # c     pressure or entropy switch
    # c
    # c     if (kvis == 0) then
    for j in range(pad+ny+pad):
        for i in range(pad+nx+pad):
            s[i,j]    = p[i,j]
    # c     else
    # c        do j=0,jb
    # c        do i=0,ib
    # c           s[i,j]    = p[i,j]/w(i,j,1)**gamma
    # c        end do
    # c        end do
    # c     end if
    # c
    # c     adaptive time step
    # c
    for j in range(pad,ny+pad):
        for i in range(pad,nx+pad):
         dpi       = abs((s[i+1,j]  -2.0*s[i,j]  +s[i-1,j])/ \
                         (s[i+1,j]  +2.0*s[i,j]  +s[i-1,j]  +slim))
         dpj       = abs((s[i,j+1]  -2.0*s[i,j]  +s[i,j-1])/ \
                         (s[i,j+1]  +2.0*s[i,j]  +s[i,j-1]  +slim))
         rfl[i,j]  = 1.0/(1.0  + min(dim(cfl,1.0), b*(dpi+dpj)))

    #       if (vt == 0.) go to 11
    # c
    # c     fixed time step
    # c
    dtmin = dtl[imin,jmin]
    if not self.local_timestepping:

        for j in range(pad,ny+pad):
            for i in range(pad,nx+pad):
                if (dtl[i,j] <= dtmin):
                    dtmin     = dtl[i,j]
                    imin      = i
                    jmin      = j

        for j in range(pad,ny+pad):
            for i in range(pad,nx+pad):
                rfl[i,j]  = dtmin/dtl[i,j]

    # c
    # c     option to rescale the dissipative coefficients
    # c
    #    11 do j=2,jl
    #       do i=2,il
    for j in range(pad,ny+pad):
        for i in range(pad,nx+pad):
            if (iprec != 0):
    # c         rfli[i,j] = math.sqrt(radJ[i,j]/radI[i,j])
    # c         rflj[i,j] = math.sqrt(radI[i,j]/radJ[i,j])
                rfli[i,j] = (radJ[i,j]/radI[i,j])**(0.25)
                rflj[i,j] = (radI[i,j]/radJ[i,j])**(0.25)
    # c         rfli[i,j] = 1.
    # c         rflj[i,j] = 1.
    
            a         = (radI[i,j]/radJ[i,j])**adis
            radI[i,j] = radI[i,j]*(1.0  +1.0/a)
            radJ[i,j] = radJ[i,j]*(1.0  +a)
            if (iprec == 0):
                rfli[i,j] = radI[i,j]/dtlc[i,j]
                rflj[i,j] = radJ[i,j]/dtlc[i,j]

    # c
    # c
    # c     reduce the artificial dissipation for viscous flows
    # c     and adjust the time step estimate
    # c
    if (kvis >= 0):

        v1        = np.sqrt(gamma)*rm/re
        v2        = 0.0
        if (kvis > 1): 
            v2 = v1

        for j in range(pad,ny+pad):
            for i in range(pad,nx+pad):
                rk        = gamma *(lv[i,j]* v1/prn + ev[i,j] * v2/prt)/w[i,j,0]
                rmu       = (lv[i,j]*v1 + ev[i,j]*v2)/w[i,j,0]
            
                # get side lengths
                [xn, yn] = edge(i,j,'n')
                [xs, ys] = edge(i,j,'s')
                [xe, ye] = edge(i,j,'e')
                [xw, yw] = edge(i,j,'w')

                xx = np.mean([xn, xs])
                yx = np.mean([yn, ys])
                xy = np.mean([xe, xw])
                yy = np.mean([ye, yw])

                # xx        = 0.5*(x(i,j,1) -x(i-1,j,1) +x(i,j-1,1) -x(i-1,j-1,1))
                # yx        = 0.5*(x(i,j,2) -x(i-1,j,2) +x(i,j-1,2) -x(i-1,j-1,2))
                # xy        = 0.5*(x(i,j,1) -x(i,j-1,1) +x(i-1,j,1) -x(i-1,j-1,1))
                # yy        = 0.5*(x(i,j,2) -x(i,j-1,2) +x(i-1,j,2) -x(i-1,j-1,2))
                dsi       = xy**2  +yy**2
                dsj       = xx**2  +yx**2
                vsi       = (rk*dsi + rmu*np.sqrt(dsi*dsj)/6.)/vol[i,j]
                vsj       = (rk*dsj + rmu*np.sqrt(dsi*dsj)/6.)/vol[i,j]
                dtv       = dtlc[i,j]+4.*(vsi+vsj)
                dtl[i,j]  = 1.0/dtv
                radI[i,j] = dim(radI[i,j],vsi)
                radJ[i,j] = dim(radJ[i,j],vsj)

    # c
    # c     set boundary values at i=1 and i=ie
    # c
    for j in range(pad,ny+pad):
        radI[1,j]   = radI[2,j]
        radI[ie,j]  = radI[il,j]
        rfl[1,j]    = rfl[2,j]
        rfl[ie,j]   = rfl[il,j]
        rfli[1,j]   = rfli[2,j]
        rfli[ie,j]  = rfli[il,j]
        rflj[1,j]   = rflj[2,j]
        rflj[ie,j]  = rflj[il,j]
        dtl[1,j]    = dtl[2,j]
        dtl[ie,j]   = dtl[il,j]

    # c
    # c     set boundary values at j=1 and j=je
    # c
    for i in range(pad-1,nx+pad+1):
        radJ[i,1]   = radJ[i,2]
        radJ[i,je]  = radJ[i,jl]
        rfl[i,1]    = rfl[i,2]
        rfl[i,je]   = rfl[i,jl]
        rfli[i,1]   = rfli[i,2]
        rfli[i,je]  = rfli[i,jl]
        rflj[i,1]   = rflj[i,2]
        rflj[i,je]  = rflj[i,jl]
        dtl[i,1]    = dtl[i,2]
        dtl[i,je]   = dtl[i,jl]

    # c
    # c     set boundary values along the cut
    # c
    for i in range(pad-1, itl+pad):
        ii        = ib  -i
        radJ[ii,1]  = radJ[i,2]
        radJ[i,1]   = radJ[ii,2]
        rfl[ii,1]   = rfl[i,2]
        rfl[i,1]    = rfl[ii,2]
        rfli[ii,1]  = rfli[i,2]
        rfli[i,1]   = rfli[ii,2]
        rflj[ii,1]  = rflj[i,2]
        rflj[i,1]   = rflj[ii,2]
        dtl[ii,1]   = dtl[i,2]
        dtl[i,1]    = dtl[ii,2]


# edge vector of control volume in positive i or j direction
def edge(workspace, i, j, side):
    i1 = i; i1 = i; j2 = j; j2 = j
    
    if side == "n":
        i1 = i  ; j1 = j+1
        i2 = i+1; j2 = j+1
    if side == "s":
        i1 = i  ; j1 = j
        i2 = i+1; j2 = j
    if side == "e":
        i1 = i; j1 = j
        i2 = i; j2 = j+1
    if side == "w":
        i1 = i+1; j1 = j
        i2 = i+1; j2 = j+1

    X = workspace.get_field('x')
    x = X[:,:,0]
    y = X[:,:,1]

    x1 = x[i1, j1]
    y1 = y[i1, j1]
    x2 = x[i2, j2]
    y2 = y[i2, j2]

    dx = x2-x1
    dy = y2-y1
    
    return [dx, dy]
