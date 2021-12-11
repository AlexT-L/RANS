# artificial dissipation using  
# first order fluxes scaled to spectral radius (calls dfluxc.f)

# append to path so we can access Field class
import sys
sys.path.append("../")

# class dependencies
from Workspace import Workspace
from Grid import Grid
from Field import Field
import numpy as np

# fortran module
import dfluxc_fort 

def dfluxc(self,ws,dw):

    # calculate artificial dissipation fluxes given a workspace

    # grab grid related parameters
    G = ws.grid
    ny = G.dims['ny']
    il = G.dims['il']
    jl = G.dims['jl']

    # flow related variabless
    w = ws.getField['w'] # state
    P = ws.getField['P'] # pressure
    porI = ws.getField['porI'] # porosity in i 
    porJ = ws.getField['porJ'] # porosity in j

    # solver related vars
    fw = ws.getField['fw'] # storage for viscous residuals?
    radI = ws.getField['radI'] # some kind of stability metric in i
    radJ = ws.getField['radJ'] # some kind of stability metric in j

    # solver params
    rfil = self.rfil
    vis0 = self.vis0

    # residuals returned in Field vw
    dfluxc_fort.dfluxc(ny,il,jl, \
                        w.vals,P.vals, \
                        porJ.vals, \
                        fw.vals, radI.vals, radJ.vals, \
                        rfil,vis0)


    # put in residuals
    dw = dw + fw.vals
