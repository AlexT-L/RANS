# viscous flux calculation wrapping fortran file nsflux.f

# append to path so we can access Field class
import sys
sys.path.append("../../../")

# class dependencies
from bin.Field import Field

# fortran module
from bin.model_funcs.fortran_versions import nsflux_fort 

def nsflux(self,ws,w,dw,rfil):

    # calculate viscous fluxes given a workspace

    # grab grid related parameters
    dims = ws.get_dims()
    il = dims['il']
    jl = dims['jl']
    ie = dims['ie']
    je = dims['je']
    itl = dims['itl']
    itu = dims['itu']
    
    # geometric parameters
    geom = ws.get_geometry()
    scal = geom['scal']
    chord = geom['chord']

    # flow related variabless
    p = ws.get('p',self.className) # pressure
    lv = ws.get('lv',self.className) # laminar viscocity
    ev = ws.get('ev',self.className) # eddy viscocity
    vw = ws.get('vw',self.className) # storage for viscous residuals

    # flow parameters
    gamma = self.params['gamma']
    mach = self.params['rm']
    Re = self.params['re']
    chord = self.params['chord']
    
    # output params (not important)
    prn = 0
    prt = 0
    
    # mode parameter says if we are at the finest mesh
    mode = ws.is_finest()

    # mesh related vars
    x = ws.get_field('x') # mesh vertices
    xc = ws.get_field('xc') # mesh centers

    # residuals returned in Field vw
    nsflux_fort.nsflux(il, jl, ie, je, \
                       w, p, lv, ev,  \
                       x, xc, \
                       vw,
                       gamma, mach, scal, \
                       Re, chord, \
                       prn, prt, mode, \
                       rfil)

    # add viscous contribution
    dw += vw