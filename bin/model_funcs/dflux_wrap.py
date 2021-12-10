# upwind biasing by artificial dissipation using  
# blended first and third order fluxes (calls dflux.f)

# append to path so we can access Field class
import sys
sys.path.append("../")

# class dependencies
from Workspace import Workspace
from Grid import Grid
from Field import Field
import numpy as np

# fortran module
import dflux_fort 

def dflux(self,ws,dw):

    # calculate artificial dissipation fluxes given a workspace

    # grab grid related parameters
    G = ws.grid
    ny = G.dims['ny']
    il = G.dims['il']
    jl = G.dims['jl']
    ie = G.dims['ie']
    je = G.dims['je']
    itl = G.dims['itl']
    itu = G.dims['itu']

    # flow related variabless
    w = ws.getField['w'] # state
    P = ws.getField['P'] # pressure
    porI = ws.getField['porI'] # porosity in i 
    porJ = ws.getField['porJ'] # porosity in j

    # mesh related vars
    x = ws.getField['x'] # mesh vertices
    xc = ws.getField['xc'] # mesh centers

    # solver related vars
    fw = ws.getField['fw'] # storage for viscous residuals?
    radI = ws.getField['radI'] # some kind of stability metric in i
    radJ = ws.getField['radJ'] # some kind of stability metric in j

    # solver params
    rfil = ws.rfil
    vis2 = ws.vis2
    vis4 = ws.vis4

    # residuals returned in Field vw
    dflux_fort.dflux(ny,il,jl,ie,je, \
                    w.vals,P.vals, \
                    porI.vals,porJ.vals, \
                    fw.vals, radI.vals, radJ.vals, \
                    rfil,vis2,vis4)


    # put in residuals
    dw = dw + fw.vals
