# fortran module
import sys
sys.path.append('../../../')
from bin.model_funcs.fortran_versions import bcwall_fort

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
    w = state
    p = get("p")
    rev = get("ev")
    
    # mesh_var
    x = workspace.get_field("x")
    
    # flo_param
    rm = model.params['rm']
    sa = model.params['sa']
    kvis = model.params['kvis']
    
    # solv_param
    isym = geom['isym']
    
    # mg_param
    mode = 1
    if workspace.is_finest():
        mode = 0

    bcwall_fort.bcwall(ny, ie, itl+1, itu+1, \
                        w, p, rev, \
                        x, \
                        kvis, isym, mode, \
                        [il,jl,ie,je])
