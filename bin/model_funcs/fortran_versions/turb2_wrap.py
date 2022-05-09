# updates eddy viscosity (ev/rev)

# append to path so we can access Field class
import sys
sys.path.append("../../../")

# class dependencies
import numpy as np
from bin.Field import Field, max, abs, isfinite

# fortran module
from bin.model_funcs.fortran_versions import turb2_fort


def turb_BL(model,ws,w,ncyc=0):
        
    # grid parameters
    [nx, ny] = ws.field_size()
    [il, jl] = [nx+1, ny+1]
    [ie, je] = [nx+2, ny+2]
    [ib, jb] = [nx+3, ny+3]
    
    dims = ws.get_dims()
    itl = dims['itl']
    itu = dims['itu']
    
    # flow related variabless
    def get(varName):
        return ws.get_field(varName, model.className)
    p = get('p') # pressure
    ev = get('ev') # eddy viscosity
    
    # mesh vars
    vol = get('vol')
    x = ws.get_field('x')

    # flow params
    gamma = model.params['gamma']
    mach = model.params['rm']
    Re = model.params['re']
    xtran = model.params['xtran']
    
    # call turb
    turb2_fort.turb2(ie,je,itl+1,itu+1, w,p,ev, x,vol, \
                     gamma,mach,Re,xtran, ncyc, [il,jl,ib,jb])
