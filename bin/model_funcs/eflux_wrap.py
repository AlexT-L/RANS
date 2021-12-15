# convective flux calculation wrapping fortran file eflux.f

import sys
sys.path.append("../../../")

# class dependencies
from Workspace import Workspace
from Grid import Grid
from Field import Field

# fortran module
import eflux_fort


def eflux(self,ws,w,dw):
    [il, jl] = ws.field_size()

    # calculate convective fluxes given a workspace
    def get(varName):
        return ws.get_field(varName, self.className)
    porJ = get('porJ') # porosity
    p = get('p') # pressure
    x = ws.get_field('x')

    # residuals returned in Field dw
    eflux_fort.eflux(w.vals,dw.vals,p.vals,x.vals,porJ.vals,il,jl)

