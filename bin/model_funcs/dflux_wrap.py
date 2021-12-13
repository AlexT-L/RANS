# upwind biasing by artificial dissipation using  
# blended first and third order fluxes (calls dflux.f)

# append to path so we can access Field class
import sys

from WorkspaceClass import WorkspaceClass
sys.path.append("../")

# class dependencies
from Workspace import Workspace
from Grid import Grid
from Field import Field
import numpy as np

# fortran module
import dflux_fort 

def dflux(model,ws,w,dw,rfil):

    # calculate artificial dissipation fluxes given a workspace

    # grab grid related parameters
    G = ws.grid
    ny = G.dims['ny']
    il = G.dims['il']
    jl = G.dims['jl']
    ie = G.dims['ie']
    je = G.dims['je']

    # flow related variabless
    def get(varName):
        return ws.get_field(varName, model.className)
    P = get('p') # pressure
    porI = get('porI') # porosity in i 
    porJ = get('porJ') # porosity in j

    # mesh related vars
    x = ws.x() # cell vertices
    xc = ws.xc # cell centers

    # solver related vars
    fw = get('fw') # storage for viscous residuals?
    radI = get('radI') # some kind of stability metric in i
    radJ = get('radJ') # some kind of stability metric in j

    # solver params
    vis2 = model.params.vis2
    vis4 = model.params.vis4

    # residuals returned in Field vw
    dflux_fort.dflux(ny,il,jl,ie,je, \
                    w.vals,P.vals, \
                    porI.vals,porJ.vals, \
                    fw.vals, radI.vals, radJ.vals, \
                    rfil,vis2,vis4)


    # put in residuals
    dw = dw + fw.vals
