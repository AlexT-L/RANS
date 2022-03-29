

def halo(bcmodel, model, workspace, state):
    """
    assign values in the ghost cells
    
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
    PAD = model.padding
    [nx, ny] = workspace.field_size()
    [il, jl] = [nx+1, ny+1]
    [ie, je] = [nx+2, ny+2]
    [ib, jb] = [nx+3, ny+3]
    [nxp,nyp]= [nx+4, ny+4]
    dims = workspace.get_dims()
    ITL = dims['itl']
    ITU = dims['itu']
    
    # shift itl and itu
    [itl, itu] = [PAD+ITL, PAD+ITU]

    # get physical parameters
    mp = model.params
    u0 = mp['u0']
    v0 = mp['v0']
    rho0 = mp['rho0']
    p0 = mp['p0']
    h0 = mp['h0']


    # c     extended values at i=0
    # c
    w[0, 1:je+1, :] = 2.0*w[1, 1:je+1, :] - w[2, 1:je+1, :]
    p[0, 1:je+1]    = 2.0*p[1, 1:je+1]    - p[2, 1:je+1]
    
    # c
    # c     extended values at i=ib
    # c
    w[ib, 1:je+1, :] = 2.0*w[ie, 1:je+1, :] - w[il, 1:je+1, :]
    p[ib, 1:je+1]    = 2.0*p[ie, 1:je+1]    - p[il, 1:je+1]

    # c
    # c     extended values at j=0
    # c     allowing for the cut in the c-mesh
    # c
    w[0:nxp, 0, :] = w[nxp-1::-1, 3, :]
    p[0:nxp, 0]    = p[nxp-1::-1, 3]

    # c
    # c     extended values at j=0 over the surface of the wing
    # c
    w[itl:itu,0,:] = 2.0*w[itl:itu,1,:] - w[itl:itu,2,:]
    p[itl:itu,0] = 2.0*p[itl:itu,1] - p[itl:itu,2]

    # c
    # c     extended values at j=jb
    # c
    w[0:nxp, jb, :] = 2.0*w[0:nxp, je, :] - w[0:nxp, jl, :]
    p[0:nxp, jb]    = 2.0*p[0:nxp, je]    - p[0:nxp, jl]