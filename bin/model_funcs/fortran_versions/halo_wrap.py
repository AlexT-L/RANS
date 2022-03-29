# fortran module
import sys
sys.path.append('../../../')
from bin.model_funcs.fortran_versions import halo_fort 

def halo(self, model, workspace, state):
    # define helper function for getting fields
    def get(varName):
        return workspace.get_field(varName, model.className)

    # get geometry dictionary
    dims = workspace.get_dims()
    
    # dims
    [nx, ny] = workspace.field_size()
    [il, jl] = [nx+1, ny+1]
    [ie, je] = [nx+2, ny+2]
    [ib, jb] = [nx+3, ny+3]
    itl = dims['itl']
    itu = dims['itu']
    
    # flo_var
    w = state
    p = get("p")

    # mesh_var
    vol = get("vol")
    
    halo_fort.halo(il, jl, ie, je, itl+1, itu+1, \
                   w, p, vol, [ib,jb])
