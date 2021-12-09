# viscous flux calculation wrapping fortran file nsflux.f

# append to path so we can access Field class
import sys
sys.path.append("../")

# class dependencies
from Workspace import Workspace
from Grid import Grid
from Field import Field

# fortran module
import nsflux_fort 

def nsflux(ws,dw):

    # calculate viscous fluxes given a workspace

    # grab grid related parameters
    G = ws.grd
    il = G.dims['il']
    jl = G.dims['jl']
    ie = G.dims['ie']
    je = G.dims['je']
    itl = G.dims['itl']
    itu = G.dims['itu']

    # flow related vars
    w = ws.flds['w'] # state
    P = ws.flds['P'] # pressure

    # mesh related vars

    # solver related vars

    # flow 

    # residuals returned in Field dw
    nsflux_fort.nsflux(il,jl,ie,je,itl,itu, \
                        w,p,rlv,rev, \
                        x,xc, \
                        vw, \
                        gamma,rm,scal,re,chord,prn,prt, \
                        mode, rfil)
