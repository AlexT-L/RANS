# import bcfar_fort, bcwall_fort, halo_fort, math
from Field import mean
import math

# set porosity
def set_porosity(self, workspace):
    # get relevant geometry parameters
    dims = workspace.get_dims()
    grid = workspace.get_grid()
    [nx, ny] = workspace.field_size()
    il = nx+1
    jl = ny+1
    itl = dims['itl']
    itu = dims['itu']

    # get porosity
    pori = workspace.get_field("pori", self.className)
    porj = workspace.get_field("porj", self.className)

    # c
    # c     set the porosity to unity
    # c
    pori[:,:] = 1.0
    porj[:,:] = 1.0

    # c
    # c     flag the wall and far field points at the j boundaries
    # c
    porj[itl:itu,0]   = 0.0

# update rev and rlv
def update_physics(self, model, workspace, state):
    pass
    ##### TO DO #####
    
    # This should call Baldwin Lomax

def bc_far(self, model, workspace, state):
    # define helper function for getting fields
    def get(varName):
        return workspace.get_field(varName, model.className)

    # get geometry dictionary
    geom = workspace.get_geometry()
    dims = workspace.get_dims()
    
    # dims
    [nx, ny] = workspace.field_size()
    il = nx+1
    jl = ny+1
    ie = nx+2
    je = ny+2
    ib = nx+3
    jb = nx+3
    itl = dims['itl']
    itu = dims['itu']
    
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
    dims = workspace.get_dims()
    
    # dims
    [nx, ny] = workspace.field_size()
    il = nx+1
    jl = ny+1
    ie = nx+2
    je = ny+2
    ib = nx+3
    jb = nx+3
    itl = dims['itl']
    itu = dims['itu']
    
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
    dims = workspace.get_dims()
    grid = workspace.get_grid()
    
    # dims
    [nx, ny] = workspace.field_size()
    il = nx+1
    jl = ny+1
    ie = nx+2
    je = ny+2
    ib = nx+3
    jb = nx+3
    itl = dims['itl']
    itu = dims['itu']
    
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
    pad = self.padding

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
    nxc = int(nx/ratio)
    nyc = int(ny/ratio)

    # parameters
    kvis = model.params['kvis']

    #     if (kvis.gt.0) then
    # c
    # c     transfer the molecular and turbulent viscosity to the coarse grid
    # c
    if kvis > 0:

        jj        = 1
        for j in range(pad, ny+pad, 2):
            jj        = jj  +1
            ii        = 1
            for i in range(pad, nx+pad, 2):
                ii        = ii  +1
                rlvc[ii, jj] = mean(rlv[i:i+2, j:j+2])
                revc[ii, jj] = mean(rlv[i:i+2, j:j+2])

        # c
        # c     set the boundary values at i=1 and i=ie
        # c
        jj        = 1
        for j in range(pad, ny+pad, ratio):
            jj        = jj  +1

            rlvc[1    ,jj] = mean(rlv[1   , j:j+ratio])
            rlvc[nxc+pad,jj] = mean(rlv[nx+pad, j:j+ratio])
            revc[1    ,jj] = mean(rev[1   , j:j+ratio])
            revc[nxc+pad,jj] = mean(rev[nx+pad, j:j+ratio])

        # c
        # c     set the boundary values at j=1 and j=je
        # c
        ii        = 1
        for i in range(pad, nx+pad, ratio):
            ii        = ii  +1

            rlvc[ii,1    ] = mean(rlv[i:i+ratio, 1   ])
            rlvc[ii,nyc+pad] = mean(rlv[i:i+ratio, ny+pad])
            rlvc[ii,1    ] = mean(rlv[i:i+ratio, 1   ])
            rlvc[ii,nyc+pad] = mean(rlv[i:i+ratio, ny+pad])


def halo_geom(self, model, workspace):
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

