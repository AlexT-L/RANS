# artificial dissipation using  
# first order fluxes scaled to spectral radius (calls dfluxc.f)

# append to path so we can access Field class
import sys
sys.path.append("../../../")

# class dependencies
from Workspace import Workspace
from Grid import Grid
from Field import Field
import numpy as np

# fortran module
import dfluxc_fort 

def dfluxc(model,ws,w,dw,fw,rfil):

    # calculate artificial dissipation fluxes given a workspace

    # grab grid related parameters
    [nx, ny] = ws.field_size()
    [il, jl] = ws.grid_size()

    # getter method for model
    def get(varName):
        return ws.get_field(varName, model.className)

    # flow related variabless
    p = get('p') # pressure
    porI = get('porI') # porosity in i 
    porJ = get('porJ') # porosity in j

    # solver related vars
    radI = get('radI') # some kind of stability metric in i
    radJ = get('radJ') # some kind of stability metric in j

    # solver params
    vis0 = model.params.vis0

    # residuals returned in Field vw
    dfluxc_fort.dfluxc(ny,il,jl, \
                        w.vals,p.vals, \
                        porJ.vals, \
                        fw.vals, radI.vals, radJ.vals, \
                        rfil,vis0)


    # put in residuals
    dw.store_sum(dw, fw)
