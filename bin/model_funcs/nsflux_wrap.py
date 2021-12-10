# viscous flux calculation wrapping fortran file nsflux.f

# append to path so we can access Field class
import sys
sys.path.append("../")

# class dependencies
from Workspace import Workspace
from Grid import Grid
from Field import Field
import numpy as np

# fortran module
import nsflux_fort 

def nsflux(self,ws,dw):

    # calculate viscous fluxes given a workspace

    # grab grid related parameters
    G = ws.grid
    il = G.dims['il']
    jl = G.dims['jl']
    ie = G.dims['ie']
    je = G.dims['je']
    itl = G.dims['itl']
    itu = G.dims['itu']

    # flow related vars
    w = ws.getField['w'] # state
    P = ws.getField['P'] # pressure
    rlv = ws.getField['rlv'] # laminar viscocity
    rev = ws.getField['rev'] # eddy viscocity

    # mesh related vars
    x = ws.getField['x'] # mesh vertices
    xc = ws.getField['xc'] # mesh centers

    # solver related vars
    vw = ws.getField['vw']


    # residuals returned in Field dw
    nsflux_fort.nsflux(il,jl,ie,je,itl,itu, \
                        w.vals,P.vals,rlv.vals,rev.vals, \
                        x.vals,xc.vals, \
                        vw.vals, \
                        self.gamma,self.rm,self.scal,self.re,self.chord, \
                        self.prn,self.prt, self.mode, self.rfil)
