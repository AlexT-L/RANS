def halo_geom(self, model, workspace):
        """Sets the geometry values in the halo
        
        Args:
        
        model:
            The physics model
        workspace:
            The Workspace
        """
    # get dimensions
    dims = workspace.get_dims()
    itl = dims['itl']
    itu = dims['itu']
    pad = self.padding
    [nx, ny] = workspace.field_size()
    nnx = pad+nx+pad
    nny = pad+ny+pad
    ib = nnx-1
    jb = nny-1
    ie = nnx-2
    je = nny-2
    il = nnx-3
    jl = nny-3

    # get model fields
    def get(varName):
        return workspace.get_field(varName, model.className)
    x = workspace.get_field('x')
    xc = get('xc')
    vol = get('vol')

    # set values at edges
    for i in range(pad,nx+pad):
        vol[i,1]    = vol[i,2]
        vol[i,je]   = vol[i,jl]
        xc[i,1,0]   = x[i-1,0,0]   +x[i-2,0,0]   -xc[i,2,0]
        xc[i,1,1]   = x[i-1,0,1]   +x[i-2,0,1]   -xc[i,2,1]
        xc[i,je,0]  = x[i-1,ny,0]  +x[i-2,ny,0]  -xc[i,jl,0]
        xc[i,je,1]  = x[i-1,ny,1]  +x[i-2,ny,1]  -xc[i,jl,1]

    for j in range(pad,ny+pad):
        vol[1,j]    = vol[2,j]
        vol[ie,j]   = vol[il,j]
        xc[1,j,0]   = x[0,j-1,0]   +x[0,j-2,0]   -xc[2,j,0]
        xc[1,j,1]   = x[0,j-1,1]   +x[0,j-2,1]   -xc[2,j,1]
        xc[ie,j,0]  = x[nx,j-1,0]  +x[nx,j-2,0]  -xc[il,j,0]
        xc[ie,j,1]  = x[nx,j-1,1]  +x[nx,j-2,1]  -xc[il,j,1]

    # reflective symmetry boundary
    for i in range(pad,itl+pad):
        vol[ib-i,1]   = vol[i,2]
        xc[ib-i,1,0]  = xc[i,2,0]
        xc[ib-i,1,1]  = xc[i,2,1]
        vol[i,1]    = vol[ib-i,2]
        xc[i,1,0]   = xc[ib-i,2,0]
        xc[i,1,1]   = xc[ib-i,2,1]

    # corner values
    vol[1,1]    = vol[2,1]
    vol[1,je]   = vol[2,je]
    vol[ie,1]   = vol[il,1]
    vol[ie,je]  = vol[il,je]

    xc[1,1,0]   = xc[2,1,0]    +xc[1,2,0]    -xc[2,2,0]
    xc[1,1,1]   = xc[2,1,1]    +xc[1,2,1]    -xc[2,2,1]
    xc[ie,1,0]  = xc[il,1,0]   +xc[ie,2,0]   -xc[il,2,0]
    xc[ie,1,1]  = xc[il,1,1]   +xc[ie,2,1]   -xc[il,2,1]
    xc[1,je,0]  = xc[2,je,0]   +xc[1,jl,0]   -xc[2,jl,0]
    xc[1,je,1]  = xc[2,je,1]   +xc[1,jl,1]   -xc[2,jl,1]
    xc[ie,je,0] = xc[il,je,0]  +xc[ie,jl,0]  -xc[il,jl,0]
    xc[ie,je,1] = xc[il,je,1]  +xc[ie,jl,1]  -xc[il,jl,1]

    # outer halo volumes
    vol[0:pad+nx+pad, 1] = vol[pad+nx+pad-1::-1, 2]
    vol[0:pad+nx+pad, 0] = vol[pad+nx+pad-1::-1, 3]

    # along airfoil
    vol[itl:itu,1] = vol[itl:itu,2]
