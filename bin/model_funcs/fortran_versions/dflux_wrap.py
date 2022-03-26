# upwind biasing by artificial dissipation using  
# blended first and third order fluxes (calls dflux.f)

# append to path so we can access Field class

import sys
sys.path.append("../../../")

# class dependencies
from bin.Workspace import Workspace
from bin.Field import Field, max, abs
import numpy as np

# fortran module
from bin.model_funcs.fortran_versions import dflux_fort


def dflux(model,ws,w,dw,rfil):
    
    
    print(max(abs(dw)))
    [nx, ny] = ws.field_size()
    [il, jl] = [nx+1, ny+1]
    [ie, je] = [nx+2, ny+2]
    [ib, jb] = [nx+3, ny+3]

    # flow related variabless
    def get(varName):
        return ws.get_field(varName, model.className)
    p = get('p') # pressure
    porI = get('porI') # porosity in i 
    porJ = get('porJ') # porosity in j

    # mesh related vars
    x = ws.get_field('x') # cell vertices
    xc = ws.get_field('xc') # cell centers

    # solver related vars
    fw = get('fw') # storage for viscous residuals?
    radI = get('radI') # some kind of stability metric in i
    radJ = get('radJ') # some kind of stability metric in j

    # solver params
    vis2 = model.params['vis2']
    vis4 = model.params['vis4']

    # residuals returned in Field vw
    dflux_fort.dflux(ny,il,jl,ie,je, \
                    w,p, \
                    porI,porJ, \
                    fw, radI, radJ, \
                    rfil,vis2,vis4)


    # put in residuals
    dw = dw + fw
