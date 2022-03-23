# viscous flux calculation wrapping fortran file nsflux.f

# append to path so we can access Field class
import sys
sys.path.append("../")

# class dependencies
from Workspace import Workspace
from Grid import Grid
from Field import Field
from Model import Model
from NavierStokes import NavierStokes

import numpy as np

# fortran module
import nsflux_fort 

def nsflux(self,ws,w,vw,rfil):

    # calculate viscous fluxes given a workspace

    # grab grid related parameters
    G = ws.grid
    il = G.dims['il']
    jl = G.dims['jl']
    ie = G.dims['ie']
    je = G.dims['je']
    itl = G.dims['itl']
    itu = G.dims['itu']

    # flow related variabless
    P = ws.get('P',self.className) # pressure
    lv = ws.get('lv',self.className) # laminar viscocity
    ev = ws.get('ev',self.className) # eddy viscocity
    #vw = ws.get('vw',self.className) # storage for viscous residuals
    #w = ws.get('w',self.className) # state

    # mesh related vars
    x = ws.get_field('x') # mesh vertices
    xc = ws.get_field('xc') # mesh centers

    # residuals returned in Field vw
    nsflux_fort.nsflux(il, jl, ie, je, \
                       w, P, lv, ev,  \
                       x, xc, \
                       vw,
                       self.flo_params['gamma'],self.flo_params['rm'],self.flo_params['scal'], \
                       self.flo_params['re'],self.flo_params['chord'], \
                       self.flo_params['prn'],self.flo_params['prt'], self.flo_params['mode'], \
                       rfil)
