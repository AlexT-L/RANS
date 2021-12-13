from numpy.core.fromnumeric import mean
# import bcfar_fort, bcwall_fort, halo_fort, math
from Field import Field

def init_state(self, model, workspace, state):

    p = workspace.get_field("p", model.className)

    # set initial values
    rho0 = model.params['rho0']
    u0 = model.params['u0']
    v0 = model.params['v0']
    h0 = model.params['h0']
    p0 = model.params['p0']

    [lenx, leny] = p.size()
    for i in range(lenx):
        for j in range(leny):
            # print("i,j")
            # print([i,j])
            state[i,j,0] = rho0
            state[i,j,1] = rho0 * u0
            state[i,j,2] = rho0 * v0
            state[i,j,3] = rho0 * h0 - p0
            p[i,j]       = p0

# set porosity
def set_porosity(self, workspace):
    # get relevant geometry parameters
    geom = workspace.get_geometry()
    grid = workspace.get_grid()
    [nx, ny] = workspace.field_size()
    il = nx+1
    jl = ny+1
    itl = grid.itl
    itu = grid.itu

    # get porosity
    pori = workspace.get_field("pori", self.className)
    porj = workspace.get_field("porj", self.className)

    # c
    # c     set the porosity to unity
    # c
    for j in range(jl):
        for i in range(il):
            pori[i,j] = 1.0
            porj[i,j] = 1.0

    # c
    # c     flag the wall and far field points at the j boundaries
    # c
    for i in range(itl,itu):
        porj[i,0]   = 0.0

# update rev and rlv
def update_physics(self, model, workspace, state):
    pass
    ##### TO DO #####
    
    # This should call Baldwin Lomax

# update stability
def update_stability(self, model, workspace, state):

    # padding
    p = self.padding
    
    # set stability parameters
    slim      = 0.001
    rlim      = 0.001
    b         = 4.0
    dtmin     = 0.0
    imin      = 0
    jmin      = 0
    iprec     = 0 # turns on gauss-seidel preconditioner

    # get relevant geometry parameters
    geom = workspace.get_geometry()
    grid = workspace.get_grid()
    [nx, ny] = workspace.field_size()
    il = nx+1
    jl = ny+1
    ie = nx+2
    je = ny+2
    ib = nx+3
    jb = nx+3
    itl = grid.itl
    itu = grid.itu

    # flo_param
    gamma = model.params['gamma']
    rm = model.params['rm']
    rho0 = model.params['rho0']
    p0 = model.params['p0']
    h0 = model.params['h0']
    c0 = model.params['co']
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
    rev = get("rev")
    rlv = get("rlv")
    radi = get("radi")
    radj = get("radj")
    rfli = get("rfli")
    rflj = get("rflj")
    rfl = get("rfl")
    dtl = get("dtl")
    w = state
    p = get("p")
    vol = get("vol")

    # local working arrays
    def get_local(varName):
            return workspace.get_field(varName, model.className)
    s = get_local("s")
    dtlc = get_local("dtlc")

    # dim helper function
    def dim(a, b):
        diff = a-b
        return max(diff, 0)

    # edge function decorator
    def edge(i, j, side):
        p = model.padding
        return workspace.edge(i-p, j-p, side)

    # c
    # c     permissible time step
    # c
    for j in range(p,ny+p):
        for i in range(p,nx+p):
            cc        = gamma*p[i,j]/max(w(i,j,1),rlim)

            # get side lengths
            [xn, yn] = edge(i,j,'n')
            [xs, ys] = edge(i,j,'s')
            [xe, ye] = edge(i,j,'e')
            [xw, yw] = edge(i,j,'w')

            xx = mean([xn, xs])
            yx = mean([yn, ys])
            xy = mean([xe, xw])
            yy = mean([ye, yw])

            # xx        = 0.5*(x(i,j,1) -x(i-1,j,1) +x(i,j-1,1) -x(i-1,j-1,1))
            # yx        = 0.5*(x(i,j,2) -x(i-1,j,2) +x(i,j-1,2) -x(i-1,j-1,2))
            # xy        = 0.5*(x(i,j,1) -x(i,j-1,1) +x(i-1,j,1) -x(i-1,j-1,1))
            # yy        = 0.5*(x(i,j,2) -x(i,j-1,2) +x(i-1,j,2) -x(i-1,j-1,2))
            qs        = (yy*w(i,j,2) -xy*w(i,j,3))/w(i,j,1)
            cs        = cc*(xy**2  +yy**2)
            radi[i,j] = abs(qs)  +math.sqrt(cs)
            qs        = (xx*w(i,j,3)  -yx*w(i,j,2))/w(i,j,1)
            cs        = cc*(xx**2  +yx**2)
            radj[i,j] = abs(qs)  +math.sqrt(cs)
            dtl[i,j]  = 1.0/(radi[i,j]  +radj[i,j])
            #dtlc[i,j] = radi[i,j]  +radj[i,j]

    # c
    # c     pressure or entropy switch
    # c
    # c     if (kvis == 0) then
    for j in range(p+ny+p):
        for i in range(p+nx+p):
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
    for j in range(p,ny+p):
        for i in range(p,nx+p):
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

        for j in range(p,ny+p):
            for i in range(p,nx+p):
                if (dtl[i,j] <= dtmin):
                    dtmin     = dtl[i,j]
                    imin      = i
                    jmin      = j

        for j in range(p,ny+p):
            for i in range(p,nx+p):
                rfl[i,j]  = dtmin/dtl[i,j]

    # c
    # c     option to rescale the dissipative coefficients
    # c
    #    11 do j=2,jl
    #       do i=2,il
    for j in range(p,ny+p):
        for i in range(p,nx+p):
            if (iprec != 0):
    # c         rfli[i,j] = math.sqrt(radj[i,j]/radi[i,j])
    # c         rflj[i,j] = math.sqrt(radi[i,j]/radj[i,j])
                rfli[i,j] = (radj[i,j]/radi[i,j])**(0.25)
                rflj[i,j] = (radi[i,j]/radj[i,j])**(0.25)
    # c         rfli[i,j] = 1.
    # c         rflj[i,j] = 1.
    
            a         = (radi[i,j]/radj[i,j])**adis
            radi[i,j] = radi[i,j]*(1.0  +1.0/a)
            radj[i,j] = radj[i,j]*(1.0  +a)
            if (iprec == 0):
                rfli[i,j] = radi[i,j]/dtlc[i,j]
                rflj[i,j] = radj[i,j]/dtlc[i,j]

    # c
    # c
    # c     reduce the artificial dissipation for viscous flows
    # c     and adjust the time step estimate
    # c
    if (kvis >= 0):

        v1        = math.sqrt(gamma)*rm/re
        v2        = 0.0
        if (kvis > 1): 
            v2 = v1

        for j in range(p,ny+p):
            for i in range(p,nx+p):
                rk        = gamma *(v1*rlv[i,j]/prn + v2*rev[i,j]/prt)/w(i,j,1)
                rmu       = (v1*rlv[i,j]+v2*rev[i,j])/w(i,j,1)
            
                # get side lengths
                [xn, yn] = edge(i,j,'n')
                [xs, ys] = edge(i,j,'s')
                [xe, ye] = edge(i,j,'e')
                [xw, yw] = edge(i,j,'w')

                xx = mean([xn, xs])
                yx = mean([yn, ys])
                xy = mean([xe, xw])
                yy = mean([ye, yw])

                # xx        = 0.5*(x(i,j,1) -x(i-1,j,1) +x(i,j-1,1) -x(i-1,j-1,1))
                # yx        = 0.5*(x(i,j,2) -x(i-1,j,2) +x(i,j-1,2) -x(i-1,j-1,2))
                # xy        = 0.5*(x(i,j,1) -x(i,j-1,1) +x(i-1,j,1) -x(i-1,j-1,1))
                # yy        = 0.5*(x(i,j,2) -x(i,j-1,2) +x(i-1,j,2) -x(i-1,j-1,2))
                dsi       = xy**2  +yy**2
                dsj       = xx**2  +yx**2
                vsi       = (rk*dsi + rmu*math.sqrt(dsi*dsj)/6.)/vol[i,j]
                vsj       = (rk*dsj + rmu*math.sqrt(dsi*dsj)/6.)/vol[i,j]
                dtv       = dtlc[i,j]+4.*(vsi+vsj)
                dtl[i,j]  = 1.0/dtv
                radi[i,j] = dim(radi[i,j],vsi)
                radj[i,j] = dim(radj[i,j],vsj)

    # c
    # c     set boundary values at i=1 and i=ie
    # c
    for j in range(p,ny+p):
        radi[1,j]   = radi[2,j]
        radi[ie,j]  = radi[il,j]
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
    for i in range(p-1,nx+p+1):
        radj[i,1]   = radj[i,2]
        radj[i,je]  = radj[i,jl]
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
    for i in range(p-1, itl+p):
        ii        = ib  -i
        radj[ii,1]  = radj[i,2]
        radj[i,1]   = radj[ii,2]
        rfl[ii,1]   = rfl[i,2]
        rfl[i,1]    = rfl[ii,2]
        rfli[ii,1]  = rfli[i,2]
        rfli[i,1]   = rfli[ii,2]
        rflj[ii,1]  = rflj[i,2]
        rflj[i,1]   = rflj[ii,2]
        dtl[ii,1]   = dtl[i,2]
        dtl[i,1]    = dtl[ii,2]


def bc_far(self, model, workspace, state):
    # define helper function for getting fields
    def get(varName):
        return workspace.get_field(varName, model.className)

    # get geometry dictionary
    geom = workspace.get_geometry()
    grid = workspace.get_grid()
    
    # dims
    [nx, ny] = workspace.field_size()
    il = nx+1
    jl = ny+1
    ie = nx+2
    je = ny+2
    ib = nx+3
    jb = nx+3
    itl = grid.itl
    itu = grid.itu
    
    # flo_var
    w = state.get_vals()
    rlv = get("rlv").get_vals()
    rev = get("rev").get_vals()
    p = get("p").get_vals()
    
    # mesh_var
    coords = workspace.get_field("x")
    x = coords.get_vals()
    coords = workspace.get_field("xc")
    xc = coords.get_vals()
    
    # out_var
    cp = workspace.get_field("cp", self.className)
    cp = cp.get_vals()
    cf = workspace.get_field("cf", self.className)
    cf = cf.get_vals()
    
    # flo_param
    gamma = model.params['gamma']
    rm = model.params['rm']
    rho0 = model.params['rho0']
    p0 = model.params['p0']
    h0 = model.params['h0']
    c0 = model.params['co']
    u0 = model.params['u0']
    v0 = model.params['v0']
    ca = model.params['ca']
    sa = model.params['sa']
    re = model.params['re']
    prn = model.params['prn']
    prt = model.params['prt']
    scal = geom['scal']
    chord = geom['chord']
    xm = geom['xm']
    ym = geom['ym']
    kvis = model.params['kvis']
    
    # solv_param
    bc = self.bc
    
    # mg_param
    mode = 1
    if workspace.is_finest():
        mode = 0
    
    # bcfar_fort.bcfar(il, jl, ie, je, itl+1, itu+1,
    #                     w, p, rlv, rev,
    #                     x, xc, 
    #                     cp, cf,
    #                     gamma,rm,rho0,p0,h0,c0,u0,v0,ca,sa,re,prn,prt,scal,chord,xm,ym,kvis,
    #                     bc,
    #                     mode)


def bc_wall(self, model, workspace, state):
    # define helper function for getting fields
    def get(varName):
        return workspace.get_field(varName, model.className)

    # get geometry dictionary
    geom = workspace.get_geometry()
    grid = workspace.get_grid()
    
    # dims
    [nx, ny] = workspace.field_size()
    il = nx+1
    jl = ny+1
    ie = nx+2
    je = ny+2
    ib = nx+3
    jb = nx+3
    itl = grid.itl
    itu = grid.itu
    
    # flo_var
    w = state.get_vals()
    p = get("p").get_vals()
    rev = get("rev").get_vals()
    
    # mesh_var
    coords = workspace.get_field("x")
    x = coords.get_vals()
    
    # flo_param
    rm = model.params['rm']
    sa = model.params['sa']
    kvis = model.params['kvis']
    
    # solv_param
    isym = geom.isym
    
    # bcwall_fort.bcwall(ny, il, ie, ib, itl+1, itu+1, 
    #                     w, p, rev,
    #                     x,
    #                     rm, sa, kvis,
    #                     isym)


def halo(self, model, workspace, state):
    # define helper function for getting fields
    def get(varName):
        return workspace.get_field(varName, model.className)

    # get geometry dictionary
    geom = workspace.get_geometry()
    grid = workspace.get_grid()
    
    # dims
    [nx, ny] = workspace.field_size()
    il = nx+1
    jl = ny+1
    ie = nx+2
    je = ny+2
    ib = nx+3
    jb = nx+3
    itl = grid.itl
    itu = grid.itu
    
    # flo_var
    w = state.get_vals()
    p = get("p").get_vals()

    # mesh_var
    coords = workspace.get_field("x")
    x = coords.get_vals()
    vol = get("vol").get_vals()
    
    # halo_fort.halo(il, jl, ie, je, ib, jb, itl+1, itu+1,
    #         w, p,
    #         x, vol)


def transfer_down(self, model, workspace1, workspace2):
    # get padding
    p = self.padding

    # get geometry dictionary
    geom1 = workspace1.get_geometry()
    geom2 = workspace2.get_geometry()

    # variables
    rev = workspace1.get_field("rev", model.className)
    rlv = workspace1.get_field("rlv", model.className)
    revc = workspace2.get_field("rev", model.className)
    rlvc = workspace2.get_field("rlv", model.className)

    # dims
    [nx, ny] = workspace1.field_size()

    # coarse mesh dims
    ratio = 2
    nxc = nx/ratio
    nyc = ny/ratio

    # parameters
    kvis = model.params['kvis']

    #     if (kvis.gt.0) then
    # c
    # c     transfer the molecular and turbulent viscosity to the coarse grid
    # c
    if kvis > 0:

        jj        = 1
        for j in range(p, ny+p, 2):
            jj        = jj  +1
            ii        = 1
            for i in range(p, nx+p, 2):
                ii        = ii  +1
                rlvc[ii, jj] = mean(rlv[i:i+2, j:j+2])
                revc[ii, jj] = mean(rlv[i:i+2, j:j+2])

        # c
        # c     set the boundary values at i=1 and i=ie
        # c
        jj        = 1
        for j in range(p, ny+p, ratio):
            jj        = jj  +1

            rlvc[1    ,jj] = mean(rlv[1   , j:j+ratio])
            rlvc[nxc+p,jj] = mean(rlv[nx+p, j:j+ratio])
            revc[1    ,jj] = mean(rev[1   , j:j+ratio])
            revc[nxc+p,jj] = mean(rev[nx+p, j:j+ratio])

        # c
        # c     set the boundary values at j=1 and j=je
        # c
        ii        = 1
        for i in range(p, nx+p, ratio):
            ii        = ii  +1

            rlvc[ii,1    ] = mean(rlv[i:i+ratio, 1   ])
            rlvc[ii,nyc+p] = mean(rlv[i:i+ratio, ny+p])
            rlvc[ii,1    ] = mean(rlv[i:i+ratio, 1   ])
            rlvc[ii,nyc+p] = mean(rlv[i:i+ratio, ny+p])
