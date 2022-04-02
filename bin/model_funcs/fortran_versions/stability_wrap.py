# artificial dissipation using  
# first order fluxes scaled to spectral radius (calls dfluxc.f)

# append to path so we can access Field class
import sys
sys.path.append("../../../")

# fortran module
from bin.model_funcs.fortran_versions import step_fort 

def stability(self,model,ws,w):

    # set stability parameters
    slim      = 0.001
    rlim      = 0.001
    b         = 4.0
    dtmin     = 0.0
    imin      = 0
    jmin      = 0
    iprec     = 0 # turns on gauss-seidel preconditioner

    # grab grid related parameters
    dims = ws.get_dims()
    [nx, ny] = ws.field_size()
    [il, jl] = [nx+1, ny+1]
    [ie, je] = [nx+2, ny+2]
    [ib, jb] = [nx+3, ny+3]
    itl = dims['itl']

    # retrieve working arrays from model
    def get(varName):
            return ws.get_field(varName, model.className)
    ev = get('ev')
    lv = get('lv')
    radi = get("radI")
    radj = get("radJ")
    rfli = get("rfli")
    rflj = get("rflj")
    rfl = get("rfl")
    dtl = get("dtl")
    p = get("p")
    vol = get("vol")
    x = ws.get_field('x')
    
    # dtlc only stored locally
    dtlc = ws.get_field('dtlc', self.className)

    # getter method for model
    def get(varName):
        return ws.get_field(varName, model.className)

    # flow related variabless
    p = get('p') # pressure
    
    # local timestep
    vt = 0
    if not self.local_timestepping:
        vt = 1

    # flo_param
    gamma = model.params['gamma']
    rm = model.params['rm']
    re = model.params['re']
    prn = model.params['prn']
    prt = model.params['prt']
    adis = model.params['adis']
    cfl = model.cfl
    kvis = model.params['kvis']

    # solver params
    vis0 = model.params['vis0']

    # residuals returned in Field vw
    step_fort.step(ie,je,itl+1, \
                        w,p,lv,ev, \
                        x,vol, \
                        rfl,rfli,rflj,radi,radj,dtl,dtlc, \
                        gamma,rm,re,prn,prt,kvis, \
                        iprec, cfl, vt, adis, \
                        [il,jl,ib,jb])

    return dtlc