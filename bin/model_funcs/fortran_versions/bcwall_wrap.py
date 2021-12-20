# fortran module
import bcwall_fort

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
    
    bcwall_fort.bcwall(ny, il, ie, ib, itl+1, itu+1, 
                        w, p, rev,
                        x,
                        rm, sa, kvis,
                        isym)
