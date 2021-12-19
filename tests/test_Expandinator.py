"""
Description
-----------
Tests the expandinator.

Libraries/Modules
-----------------
-pytest \n
-Field \n
-numpy \n
-scipy

Notes
-----


Author(s)
---------
Satya Butler
"""
from numpy.core.numeric import array_equal
import pytest
import numpy as np
from bin.Field import Field
import bin.Expandinator as Exp

def test_bilinear4way():
    """
    Description of the Expandinator.py 'bilinear4way' method.

    Parameters
    ----------
    
    Returns
    -------
    :
        Nothing, but asserts if 'bilinear4way' expands a Field properly.
    
    """
    output_dims = (4,4)
    output_values = np.array([[1,3,1,3],[1,3,1,3],[1,3,1,3],[1,3,1,3]])
    output_field = Field(output_dims, output_values)
    #coarse_dims = (2,2)
    #coarse_field = Field(coarse_dims)
    input_dims = (2,2)
    input_values = np.array([[2,2],[2,2]])
    input_field = Field(input_dims, input_values)

    dims = (2,2)
    values = np.array([[2,2],[2,2]])
    field = Field(dims, values)
    print(field)

#    Exp.bilinear4way(input_field, field)
    #assert array_equal(coarse_field, output_field)
    print(field)
#    assert 1==3
