

def halo(bcmodel, model, workspace, state):

    # retrieve variables
    w = state
    p = workspace.get_field('p', model.className)
    x = workspace.get_field('x')
    
    # set geometry parameters
    pad = model.padding
    [nx, ny] = workspace.field_size()
    [il, jl] = [nx+1, ny+1]
    [ie, je] = [nx+2, ny+2]
    dims = workspace.get_dims()
    itl = dims['itl']
    itu = dims['itu']

    # get physical parameters
    mp = model.params
    u0 = mp['u0']
    v0 = mp['v0']
    rho0 = mp['rho0']
    p0 = mp['p0']
    h0 = mp['h0']


    # c     extended values at i=0
    # c
    w[0, 1:je+1, :] = 2*w[1, 1:je+1, :] - w[2, 1:je+1, :]
    p[0, 1:je+1]    = 2*p[1, 1:je+1]    - p[2, 1:je+1]
    
    # c
    # c     extended values at i=ib
    # c
    w[ie+1, 1:je+1, :] = 2*w[ie, 1:je+1, :] - w[il, 1:je+1, :]
    p[ie+1, 1:je+1]    = 2*p[ie, 1:je+1]    - p[il, 1:je+1]

    # c
    # c     extended values at j=0
    # c     allowing for the cut in the c-mesh
    # c
    w[0:pad+nx+pad, 0, :] = w[pad+nx+pad-1::-1, 3, :]
    p[0:pad+nx+pad, 0]    = p[pad+nx+pad-1::-1, 3]

    # c
    # c     extended values at j=0 over the surface of the wing
    # c
    w[itl:itu,0,:] = 2*w[itl:itu,1,:] - w[itl:itu,2,:]
    p[itl:itu,0] = 2*p[itl:itu,1] - p[itl:itu,2]

    # c
    # c     extended values at j=jb
    # c
    w[0:pad+nx+pad, je+1, :] = 2*w[0:pad+nx+pad, je, :] - w[0:pad+nx+pad, jl, :]
    p[0:pad+nx+pad, je+1]    = 2*p[0:pad+nx+pad, je]    - p[0:pad+nx+pad, jl]