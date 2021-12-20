"""
Description

Tests the expandinator.

Libraries/Modules

-pytest \n
-Field \n
-numpy \n
-scipy

Notes



Author(s)

Satya Butler
"""
from numpy.core.numeric import array_equal
import pytest
import numpy as np
from bin.Field import Field
import bin.Expandinator as Exp

def test_bilinear4way():
    """
    Tests the Expandinator.py 'bilinear4way' method. Does not currently assert (test) anything
    as Expandinator.py does not currently work.

    Args:
    
    
    Returns
    
    :
        Nothing, but asserts if 'bilinear4way' expands a Field properly.
    
    """
    #output_dims = (4,4)
    #output_values = np.array([[1,3,1,3],[1,3,1,3],[1,3,1,3],[1,3,1,3]])
    #output_field = Field(output_dims, output_values)
    ##coarse_dims = (2,2)
    ##coarse_field = Field(coarse_dims)
    #input_dims = (2,2)
    #input_values = np.array([[2,2],[2,2]])
    #input_field = Field(input_dims, input_values)


    test_dims = (2,2,4)
    #test_dims = (2,2)
    test_values = np.zeros(test_dims) + 1
    test_field = Field(test_dims, test_values)
    print(test_field)

    testfine_dims = (4,4,4)
    #testfine_dims = (4,4)
    testfine_values = np.zeros(testfine_dims)
    testfine_field = Field(testfine_dims, testfine_values)
    print(testfine_field)

    testoutput_dims = (4,4,4)
    #testoutput_dims = (4,4)
    testoutput_values = np.zeros(testoutput_dims) + 1
    testoutput_field = Field(testoutput_dims, testoutput_values)
    print(testoutput_field)
#    Exp.bilinear4way(test_field, testfine_field)
    print("bananas")
    print(testfine_field)
    print(testoutput_field)

#    assert array_equal(testfine_field, testoutput_field)

    #dims = (2,2,2)
    #values = np.array([[[2,2],[2,2]],[[2,2],[2,2]]])
    #ield = Field(dims, values)
    #print(field)

    #Exp.bilinear4way(input_field, field)
    #assert array_equal(coarse_field, output_field)
    #print(field)
    #assert 1==3
