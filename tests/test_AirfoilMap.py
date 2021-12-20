"""
Tests the AirfoilMap object to make sure all the outputs are as expected

Libraries/Modules:
    pytest \n
    AerfoilMap

Notes:
    Runs the following tests:\n
    1. Checks that x,xc and vol are Field objects  \n
    2. Checks that vol is positive values\n

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



    

    
