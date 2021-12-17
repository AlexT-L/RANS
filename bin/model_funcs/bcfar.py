

def far_field(bcmodel, model, workspace, state):

    # retrieve variables
    w = state
    p = workspace.get_field('p', model.className)
    x = workspace.get_field('x')
    
    # set geometry parameters
    pad = model.padding
    [nx, ny] = workspace.field_size()
    [il, jl] = [nx+1, ny+1]
    [ie, je] = [nx+2, ny+2]

    # get physical parameters
    mp = model.params
    u0 = mp['u0']
    v0 = mp['v0']
    rho0 = mp['rho0']
    p0 = mp['p0']
    h0 = mp['h0']

    # c
    # c     extrapolate all quantities at the outflow boundary
    # c     and fix all quantities at the inflow boundary
    # c     if the flow is supersonic
    # c

    xx = x[1:nx+1, ny, 0] - x[0:nx, ny, 0]
    yx = x[1:nx+1, ny, 1] - x[0:nx, ny, 1]
    qn = (xx*v0 - yx*u0) > 0

    w[pad:pad+nx, je, 0] += qn*w[pad:pad+nx, jl, 0]
    w[pad:pad+nx, je, 0] += (1-qn)*rho0
    w[pad:pad+nx, je, 1] = qn*w[pad:pad+nx, jl, 1] + (1-qn)*(rho0*u0)
    w[pad:pad+nx, je, 2] = qn*w[pad:pad+nx, jl, 2] + (1-qn)*(rho0*v0)
    w[pad:pad+nx, je, 3] = qn*w[pad:pad+nx, jl, 3] + (1-qn)*(rho0*h0 - p0)


    print(p.shape())
    print(p[pad:pad+nx, jl].shape())
    print((1-qn).shape())

    p[pad:pad+nx, je]    = qn*p[pad:pad+nx, jl]    + (1-qn)*p0

    w[1, pad:pad+ny, :] = w[2, pad:pad+ny, :]
    p[1, pad:pad+ny]    = p[1, pad:pad+ny]

    w[ie, pad:pad+ny, :] = w[il, pad:pad+ny, :]
    p[ie, pad:pad+ny]    = p[il, pad:pad+ny]

    # c
    # c     enforce bc at corner cells
    # c
    w[1,je,:] = w[2,je,:]
    p[1,je]     = p[2,je]

    w[ie,je,:]  = w[il,je,:]
    p[ie,je]    = p[il,je]

    w[1,1,:]    = w[ie,2,:]
    p[1,1]      = p[ie,2]

    w[ie,1,:]   = w[1,2,:]
    p[ie,1]     = p[1,2]

