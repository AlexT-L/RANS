from Field import norm, pos_diff

def wall(bcmodel, model, workspace, state):

    # retrieve variables
    w = state
    p = workspace.get_field('p', model.className)
    x = workspace.get_field('x')
    
    # set geometry parameters
    pad = model.padding
    [nx, ny] = workspace.field_size()
    [il, jl] = [nx+1, ny+1]
    [ie, je] = [nx+2, ny+2]
    [ib, jb] = [nx+3, ny+3]

    dims = workspace.get_dims()
    itl = dims['itl']
    itu = dims['itu']

    geom = workspace.get_geometry()
    isym = geom['isym']

    # get physical parameters
    mp = model.params
    u0 = mp['u0']
    v0 = mp['v0']
    rho0 = mp['rho0']
    p0 = mp['p0']
    h0 = mp['h0']
    kvis = mp['kvis']


    # c
    # c     set values below the cut in the c mesh
    # c
    w[0:itl, 1,:]    = w[ib:itu-1:-1, 2,:]
    p[0:itl, 1]      = p[ib:itu-1:-1, 2]
    w[ib:itu-1:-1, 1,:] = w[0:itl, 2,:]
    p[ib:itu-1:-1, 1]   = p[0:itl, 2]
    
    # c
    # c     set the viscous no-slip boundary condition
    # c
    if kvis > 0:
        w[itl:itu,1,:] = w[itl:itu,2,:]
        p[itl:itu,1]   = p[itl:itu,2]

    # c
    # c     set the euler boundary condition
    # c
    else:
        a = x[itl:itu+1,1,0] - x[itl-1:itu,1,0]
        b = x[itl:itu+1,1,1] - x[itl-1:itu,1,1]
        sx = 1/norm(a,b)
        xx = a*sx
        yx = b*sx
        
        # c
        # c     calculate the wall curvature
        # c
        xxx = (xx[itl+1:itu+1] - xx[itl-1:itu-1])
        yxx = (yx[itl+1:itu+1] - yx[itl-1:itu-1])

        # c
        # c     special treatment of a sharp leading edge
        # c
        if isym == 2:
            ile       = int((pad+nx+pad)/2) - 1
            a         = x(ile,2,0)  -x(ile,1,0)
            b         = x[ile,2,1]  -x[ile,1,1]
            sxn       = 1./norm(a,b)
            xxn       = a*sxn
            yxn       = b*sxn
            xxx[ile]    = .5*(xxn  -xx[ile-1])
            yxx[ile]    = .5*(yxn  -yx[ile-1])
            xxx[ile+1]  = .5*(xx[ile+2]  +xxn)
            yxx[ile+1]  = .5*(yx[ile+2]  +yxn)
        
        # c
        # c     extrapolation using normal pressure gradient at surface
        # c
        qt = xx*w[itl:itu,2,1] + yx*w[itl:itu,2,2]
        qn = yx*w[itl:itu,2,1] + xx*w[itl:itu,2,2]

        w[itl:itu,1,0] = w[itl:itu,2,0]
        w[itl:itu,1,1] = xx*qt - yx*qn
        w[itl:itu,1,2] = xx*qn - yx*qt
        w[itl:itu,1,3] = w[itl:itu,2,3]

        rho_a   = w[itl:itu,2,0] + w[itl:itu,1,0]
        rhoU_a  = w[itl:itu,2,1] + w[itl:itu,1,1]
        rhoV_a  = w[itl:itu,2,2] + w[itl:itu,1,2]

        px      = (p[itl+1:itu+1,2] - p[itl-1:itu-1,2])/2
        xy      = (x[itl:itu,2,0] - x[itl:itu,1,0] + x[itl-1:itu-1,2,0] - x[itl-1:itu-1,1,0])/2
        yy      = (x[itl:itu,2,1] - x[itl:itu,1,1] + x[itl-1:itu-1,2,1] - x[itl-1:itu-1,1,1])/2

        qs      = (yy*rhoU_a - xy*rhoV_a)/rho_a
        gxy     = xx*xy + yx*yy
        qxy     = qs*(xxx*rhoV_a - yxx*rhoU_a)/2

        if ny < 3:
            py = p[itl:itu,3] - p[itl:itu,2]

        p[itl:itu,1] = pos_diff(p[itl:itu,2], py)
        w[itl:itu,1,3] = w[itl:itu,2,3] + p[itl:itu,2] - p[itl:itu,1]
    
    # c
    # c     update eddy viscosity boundary conditions
    # c   
    if kvis and not workspace.is_finest():
        rev = workspace.get_field('rev', model.className)
        rev[itl:itu,1] = -rev[itl:itu,2]
