# convective flux calculation wrapping fortran file eflux.f

import sys
sys.path.append("../")

# class dependencies
from Workspace import Workspace
from Grid import Grid
from Field import Field

# fortran module
import eflux_fort 

def eflux(self,ws,w,dw):

    G = ws.grid
    il = G.dims['il']
    jl = G.dims['jl']

    # calculate convective fluxes given a workspace
    def get(varName):
        return ws.get_field(varName, self.className)
    porJ = get('porJ') # porosity
    P = get('p') # pressure
    x = ws.x()

    # residuals returned in Field dw
    eflux_fort.eflux(w.vals,dw.vals,P.vals,x.vals,porJ.vals,il,jl)

