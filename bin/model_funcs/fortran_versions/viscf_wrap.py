# upwind biasing by artificial dissipation using  
# blended first and third order fluxes (calls dflux.f)

# append to path so we can access Field class

import sys
sys.path.append("../../../")

# class dependencies
from bin.Field import Field, max, abs
from bin.model_funcs.Viscosity import KTURB 

# fortran module
from bin.model_funcs.fortran_versions import viscf_fort, turb2_wrap


def viscosity(model,ws,w,ncyc=0):
    
    # grid parameters
    [nx, ny] = ws.field_size()
    [il, jl] = [nx+1, ny+1]
    [ie, je] = [nx+2, ny+2]
    [ib, jb] = [nx+3, ny+3]

    dims = ws.get_dims()
    itl = dims['itl']
    itu = dims['itu']
    
    geom = ws.get_geometry()
    scal = geom['scal']
    chord = geom['chord']

    # flow related variabless
    def get(varName):
        return ws.get_field(varName, model.className)
    p = get('p') # pressure

    # mesh related vars
    x = ws.get_field('x')
    xc = get('xc')

    # solver related vars
    rev = get('ev') # eddy viscosity
    rlv = get('lv') # laminar viscosity

    # flow params
    par = model.params
    gamma = par['gamma']
    mach = par['rm']
    Re = par['re']
    t0 = par['t0']
    rmu0 = par['mu0']
    xtran = par['xtran']
    kvis = par['kvis']
    
    # viscosity model toggle
    kturb = KTURB
    
    # mg_param
    mode = 1
    if ws.is_finest():
        mode = 0
        
    assert max(abs(w)) > 0
    
    # call turb2
    if kturb == 1:
        turb2_wrap.turb_BL(model,ws,w,ncyc)

    # call viscf
    viscf_fort.viscf(ny,ie,je,itl+1,itu+1, \
                    w,p,rlv,rev,x,xc, \
                    gamma,mach,Re,t0,rmu0,xtran,scal,chord, \
                    kvis,kturb, \
                    ncyc,mode,[il,jl,ib,jb])
