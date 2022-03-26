# convective flux calculation wrapping fortran file eflux.f

import sys
sys.path.append("../../../")

# class dependencies
from bin.Workspace import Workspace
from bin.Grid import Grid
from bin.Field import Field
import numpy as np

# fortran module
from bin.model_funcs.fortran_versions import eflux_fort


def eflux(model,ws,w,dw):
    [nx, ny] = ws.field_size()
    [il, jl] = [nx+1, ny+1]
    [ib, jb] = [nx+3, ny+3]

    # calculate convective fluxes given a workspace
    def get(varName):
        return ws.get_field(varName, model.className)
    porJ = get('porJ') # porosity
    p = get('p') # pressure
    x = ws.get_field('x')
    
    # residuals returned in Field dw
    eflux_fort.eflux(w,dw,p,x,porJ,il,jl,ib,jb)

