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

    # get geometry dictionary
    geom = ws.get_geometry()
    
    # dims
    [nx, ny] = ws.field_size()
    [il, jl] = [nx+1, ny+1]
    [ie, je] = [nx+2, ny+2]
    [ib, jb] = [nx+3, ny+3]
    
    # geometric parameters
    geom = ws.get_geometry()
    scal = geom['scal']
    chord = geom['chord']

    # flow related variabless
    def get(varName):
        return ws.get_field(varName, self.className)
    p = get('p') # pressure
    lv = get('lv') # laminar viscocity
    ev = get('ev') # eddy viscocity
    vw = get('vw') # storage for viscous residuals

    # flow parameters
    gamma = self.params['gamma']
    mach = self.params['rm']
    Re = self.params['re']
    
    # output params (not important)
    prn = 0
    prt = 0
    
    # mg_param
    mode = 1
    if ws.is_finest():
        mode = 0

    # mesh related vars
    x = ws.get_field('x') # mesh vertices
    xc = ws.get_field('xc', self.className) # mesh centers

    # residuals returned in Field vw
    nsflux_fort.nsflux(ie, je, \
                       w, p, lv, ev,  \
                       x, xc, \
                       vw,
                       gamma, mach, scal, \
                       Re, chord, \
                       prn, prt, \
                       rfil, [il,jl,ib,jb])

    # add viscous contribution
    dw += vw