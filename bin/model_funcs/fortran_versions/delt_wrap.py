# updates eddy viscosity (ev/rev)

# append to path so we can access Field class
import sys
sys.path.append("../../../")

# class dependencies
import numpy as np
from bin.Field import Field, max, abs, isfinite

# fortran module
from bin.model_funcs.fortran_versions import delt_fort


def thickness(model,ws,w):
    
    # grid parameters
    [nx, ny] = ws.field_size()
    [il, jl] = [nx+1, ny+1]
    [ie, je] = [nx+2, ny+2]
    [ib, jb] = [nx+3, ny+3]
    [nxp,nyp]= [nx+4, ny+4]
    
    # create variable to store boundary layer thickness
    ynot = Field.create(nxp)
    dsti = Field.create(nxp)
    
    # mesh vars
    x = ws.get_field('x')
    xc = ws.get_field('xc', model.className)

    # call turb
    delt_fort.delt(ny,w,ynot,dsti,x,xc,[il,jl,ib,jb])
    
    return [ynot, dsti]

