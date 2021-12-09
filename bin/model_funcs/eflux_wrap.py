# convective flux calculation wrapping fortran file eflux.f

# class dependencies
from Workspace import Workspace
from Grid import Grid
from Field import Field

# fortran module
import eflux_fort 

def eflux(ws,dw):

    # calculate convective fluxes given a workspace
    G = ws.grd
    w = ws.flds['w'] # state
    porJ = ws.flds['porJ'] # porosity
    P = ws.flds['P'] # pressure

    # residuals returned in Field dw
    eflux_fort.eflux(w.vals,dw.vals,P.vals,G.X.vals,porJ.vals)

