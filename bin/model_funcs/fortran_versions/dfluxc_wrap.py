# artificial dissipation using  
# first order fluxes scaled to spectral radius (calls dfluxc.f)

# append to path so we can access Field class
import sys
sys.path.append("../../../")

# class dependencies
from bin.Field import Field

# fortran module
from bin.model_funcs.fortran_versions import dfluxc_fort 

def dfluxc(model,ws,w,dw,rfil):

    # calculate artificial dissipation fluxes given a workspace

    # grab grid related parameters
    [nx, ny] = ws.field_size()
    [il, jl] = [nx+1, ny+1]
    [ie, je] = [nx+2, ny+2]
    [ib, jb] = [nx+3, ny+3]

    # getter method for model
    def get(varName):
        return ws.get_field(varName, model.className)

    # flow related variabless
    p = get('p') # pressure
    porI = get('porI') # porosity in i 
    porJ = get('porJ') # porosity in j

    # solver related vars
    fw = get('fw') # storage for viscous residuals?
    radI = get('radI') # some kind of stability metric in i
    radJ = get('radJ') # some kind of stability metric in j

    # solver params
    vis0 = model.params['vis0']

    # residuals returned in Field vw
    dfluxc_fort.dfluxc(ny, \
                        w,p, \
                        porJ, \
                        fw, radI, radJ, \
                        rfil,vis0,
                        [il,jl,ib,jb])


    # put in residuals
    dw += fw
