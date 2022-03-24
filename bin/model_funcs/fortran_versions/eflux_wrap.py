# convective flux calculation wrapping fortran file eflux.f

import sys
sys.path.append("../../../")

# class dependencies
from bin.Workspace import Workspace
from bin.Grid import Grid
from bin.Field import Field

# fortran module
from bin.model_funcs.fortran_versions import eflux_fort


def eflux(self,ws,w,dw):
    [il, jl] = ws.field_size()

    # calculate convective fluxes given a workspace
    def get(varName):
        return ws.get_field(varName, self.className)
    porJ = get('porJ') # porosity
    p = get('p') # pressure
    x = ws.get_field('x')

    # residuals returned in Field dw
    eflux_fort.eflux(w,dw,p,x,porJ,il,jl)

