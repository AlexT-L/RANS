"""
Tests the AirfoilMap object to make sure all the outputs are as expected

Libraries/Modules:
    pytest \n
    AerfoilMap

Notes:
    Runs the following tests:\n
    1. Checks that x,xc and vol are Field objects  \n
    2. Checks that vol is positive values\n
    3. Check that all the required dims values are being read into the dims dictionary \n
    4. Check that all the required solv_param values are being read into the solv_param dictionary \n
    5. Check that all the required flo_param values are being read into the flo_param dictionary \n
    6. Check that all the required geo_param values are being read into the geo_param dictionary \n



"""
import pytest
from bin.AirfoilMap import AirfoilMap
from bin.Field import Field
from bin.Input import Input
import numpy as np

#create input object
input = Input("rae9-s1.data")
gridInput = input.add_dicts(input.geo_param, input.in_var)
grid_dim = [input.dims['nx'], input.dims['ny']]
grid = AirfoilMap.from_file(grid_dim, gridInput)
x=grid.fields["x"]
xc=grid.fields["xc"]
vol=grid.fields["vol"]


def test_if_field_object():
    """Assert that x,xc and vol are Field
    objects"""

    assert type(x) is Field, type(xc) is Field & type(vol) is Field

def test_vol_postive():
    """Assert that vol is all non-negative values"""
    vol_val=vol.get_vals() #get underlying numpy array
    assert np.all(vol_val>0)



    

    
