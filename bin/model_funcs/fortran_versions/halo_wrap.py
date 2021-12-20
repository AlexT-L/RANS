# fortran module
import halo_fort

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
    
    halo_fort.halo(il, jl, ie, je, ib, jb, itl+1, itu+1,
            w, p,
            x, vol)
