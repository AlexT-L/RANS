from bin.Field import Field, norm, pos_diff

def wall(bcmodel, model, workspace, state):
    """
    set values at the wall
    
    Args:
        bcmodel (NS_Arifoil): boundary condition object
        model (NavierStokes): physics model
        workspace (Workspace): the relevant fields
        state (Field): containing the density, x-momentum, y-momentum, and energy
        
    """
    # retrieve variables
    w = state
    p = workspace.get_field('p', model.className)
    x = workspace.get_field('x')
    
    # set geometry parameters
    pad = model.padding
    [nx, ny] = workspace.field_size()
    [nxp, nyp] = [pad+nx+pad, pad+ny+pad]
    [il, jl] = [nx+1, ny+1]
    [ie, je] = [nx+2, ny+2]
    [ib, jb] = [nx+3, ny+3]

    dims = workspace.get_dims()
    itl = dims['itl'] + pad
    itu = dims['itu'] + pad

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
        a = x[itl:itu+2,1,0] - x[itl-1:itu+1,1,0]
        b = x[itl:itu+2,1,1] - x[itl-1:itu+1,1,1]
        sx = 1/norm(a,b)
        xx = a*sx
        yx = b*sx
        
        # c
        # c     calculate the wall curvature
        # c
        xxx = (xx[2:itu-itl+2] - xx[0:itu-itl])/2
        yxx = (yx[2:itu-itl+2] - yx[0:itu-itl])/2

        # c
        # c     special treatment of a sharp leading edge
        # c
        if isym == 2:
            ile       = int(ie/2)
            a         = x[ile,2,0]  -x[ile,1,0]
            b         = x[ile,2,1]  -x[ile,1,1]
            sxn       = 1/norm(a,b)
            xxn       = a*sxn
            yxn       = b*sxn
            xxx[ile-itl-1]  = .5*(xxn  -xx[ile-itl-1])
            yxx[ile-itl-1]  = .5*(yxn  -yx[ile-itl-1])
            xxx[ile-itl]    = .5*(xx[ile-itl+2]  +xxn)
            yxx[ile-itl]    = .5*(yx[ile-itl+2]  +yxn)
        
        # c
        # c     extrapolation using normal pressure gradient at surface
        # c
        xx = xx[1:itu-itl+1]
        yx = yx[1:itu-itl+1]
        
        qt = xx*w[itl:itu,2,1] + yx*w[itl:itu,2,2]
        qn = yx*w[itl:itu,2,1] - xx*w[itl:itu,2,2]

        w[itl:itu,1,0] = w[itl:itu,2,0]
        w[itl:itu,1,1] = xx*qt - yx*qn
        w[itl:itu,1,2] = xx*qn + yx*qt
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
        py        = (px*gxy  +qxy)*sx[1:itu-itl+1]
        
        if ny < 3:
            py = p[itl:itu,3] - p[itl:itu,2]

        p[itl:itu,1] = pos_diff(p[itl:itu,2], py)
        w[itl:itu,1,3] = w[itl:itu,2,3] + p[itl:itu,2] - p[itl:itu,1]
    
    # c
    # c     update eddy viscosity boundary conditions
    # c   
    if kvis and not workspace.is_finest():
        ev = workspace.get_field('ev', model.className)
        ev[itl:itu,1] = - ev[itl:itu,2]
