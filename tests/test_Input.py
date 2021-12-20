
"""
Tests the Input object to make sure all inputs as as expected

Libraries/Modules:
    pytest \n
    Input

Notes:
    Runs the following tests:\n
    1. Checks that the input airfoil geometry is a closed curve  \n
    2. Checks that the airfoil geomtry is of the right length \n
    3. Check that all the required dims values are being read into the dims dictionary \n
    4. Check that all the required solv_param values are being read into the solv_param dictionary \n
    5. Check that all the required flo_param values are being read into the flo_param dictionary \n
    6. Check that all the required geo_param values are being read into the geo_param dictionary \n

Author(s)

Vedin Dewan \n

"""
import pytest
from bin.Input import Input
import numpy as np

#create input object
input = Input("rae9-s1.data")
#airfoil geometry
xn    = input.in_var["xn"]
yn    = input.in_var["yn"]

def test_geometry_closed():
    """
    Asserts that geometry is a closed curve
    
    """
    
    # check if first and last value of xn
    # and yn array are same
    first_x =xn[0]
    last_x  =xn[-1]
    first_y =yn[0]
    last_y  =yn[-1]
    
    # Assert if closed curve
    assert first_x == last_x, first_y==last_y

def test_geom_length():
    """
    Asserts that the number of points on the geometry
    is as expected
    
    """
    #number of points on upper and lower surface
    nu=input.geo_param["nu"]
    nl=input.geo_param["nl"]
    #total number of points
    total=nu+nl-1
    #length of xn and yn arrays 
    l_x=len(xn)
    l_y=len(yn)

    assert l_x==total,l_y==total

def test_dims():
    """
    Asserts that the all the required values in the 
    dims dictionary are being read in
    
    """
    required_dims=["nx","ny"]
    dims=input.dims
    read_in_dims=[*dims]

    assert set(required_dims)==set(read_in_dims)

def test_solv_param():
    """
    Asserts that the all the required values in the 
    solv_param dictionary are being read in
    
    """
    required_solv_param=["fcyc","fprnt","fout","ftim","gprnt","hprnt","hmesh",
           "cflf","cflim","vis2","vis4","adis","qdis","bc","hmf",
            "cstp","cdis","smoopi","smoopj","ksmoop","vt",
            "iprec","epsf","epsc","diag",
            "cflc","fcoll","fadd","vis0","hmc","fbc","lcyc","mstage"]
    solv_param=input.solv_param
    read_in_solv_param=[*solv_param]

    assert set(required_solv_param)==set(read_in_solv_param)

def test_flo_param():
    """
    Asserts that the all the required values in the 
    flo_param dictionary are being read in
    
    """
    required_flo_param=["rm","al","fcl","clt","cd0","re","prn","prt","t0","xtran","kvis",
                   "alpha","ca","sa","gamma","rho0","p0","c0","ei0","u0","v0","h0","mu0"]
    flo_param=input.flo_param
    read_in_flo_param=[*flo_param]

    assert set(required_flo_param)==set(read_in_flo_param)

def test_geo_param():
    """
    Asserts that the all the required values in the 
    geo_param dictionary are being read in
    
    """
    required_geo_param=["boundx","boundy","bunch","xte","ylim1","ylim2","ax","ay","sy",
           "aplus","ncut","isym","nu","nl","trail","slopt","xsing","ysing","nn"]
    geo_param=input.geo_param
    read_in_geo_param=[*geo_param]

    assert set(required_geo_param)==set(read_in_geo_param)




    
    


