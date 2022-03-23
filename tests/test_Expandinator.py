"""
Description

Tests the expandinator.

Libraries/Modules

-pytest \n
-Field \n
-numpy \n
-scipy

Notes

"""
import pytest
import numpy as np
from bin.Field import Field, array_equal
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
    #output_field = Field.create(output_dims, output_values)
    ##coarse_dims = (2,2)
    ##coarse_field = Field.create(coarse_dims)
    #input_dims = (2,2)
    #input_values = np.array([[2,2],[2,2]])
    #input_field = Field.create(input_dims, input_values)


    test_dims = (4,4,4)
    test_field = Field.create(test_dims, 1)

    testfine_dims = (8,8,4)
    testfine_field = Field.create(testfine_dims, 2)

    testoutput_field = Field.create(testfine_dims, 1)
    testoutput_field[:,0] = 0
    Exp.bilinear4way(test_field, testfine_field)
    # print(testoutput_field)
    # print(testfine_field)
    # assert array_equal(testfine_field, testoutput_field)

    #dims = (2,2,2)
    #values = np.array([[[2,2],[2,2]],[[2,2],[2,2]]])
    #ield = Field.create(dims, values)
    #print(field)

    #Exp.bilinear4way(input_field, field)
    #assert array_equal(coarse_field, output_field)
    #print(field)
    #assert 1==3
