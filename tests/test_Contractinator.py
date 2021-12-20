"""
Description
-----------
Tests the contractinator.

Libraries/Modules
-----------------
-pytest \n
-Field \n
-numpy \n

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
import bin.Contractinator as Ctr

def test_simple():
    """
    Tests the Contractinator.py 'simple' method.

    Parameters
    ----------
    
    Returns
    -------
    :
        Nothing, but asserts if 'simple' deletes items from the Field as it should.
    """
    # test 2D zero field
    input_dims = (4,4)
    input_zeros = np.zeros(input_dims)
    input_field = Field(input_dims, input_zeros)
    coarse_dims = (2,2)
    coarse_field = Field(coarse_dims)
    output_dims = (2,2)
    output_zeros = np.zeros(output_dims)
    output_field = Field(output_dims, output_zeros)

    Ctr.simple(input_field, coarse_field)
    assert array_equal(coarse_field, output_field)

    # test 3D zero field
    
    input_dims = (4,4,4)
    input_zeros = np.zeros(input_dims)
    input_field = Field(input_dims, input_zeros)
    coarse_dims = (2,2,4)
    coarse_field = Field(coarse_dims)
    output_dims = (2,2,4)
    output_zeros = np.zeros(output_dims)
    output_field = Field(output_dims, output_zeros)

    Ctr.simple(input_field, coarse_field)
    assert array_equal(coarse_field, output_field)

    # test 2D ones field
    input_dims = (4,4)
    input_values = np.array([[1,1,1,1],[1,1,1,1],[1,1,1,1],[1,1,1,1]])
    input_field = Field(input_dims, input_values)
    coarse_dims = (2,2)
    coarse_field = Field(coarse_dims)
    output_dims = (2,2)
    output_values = np.array([[1,1],[1,1]])
    output_field = Field(output_dims, output_values)

    Ctr.simple(input_field, coarse_field)
    assert array_equal(coarse_field, output_field)

    # test 2D one-two field
    input_dims = (4,4)
    input_values = np.array([[1,2,1,2],[1,2,1,2],[1,2,1,2],[1,2,1,2]])
    input_field = Field(input_dims, input_values)
    coarse_dims = (2,2)
    coarse_field = Field(coarse_dims)
    output_dims = (2,2)
    output_values = np.array([[1,1],[1,1]])
    output_field = Field(output_dims, output_values)

    Ctr.simple(input_field, coarse_field)
    assert array_equal(coarse_field, output_field)

    # test 2D field
    input_dims = (4,4)
    input_values = np.array([[1,1,1,1],[2,2,2,2],[3,3,3,3],[4,4,4,4]])
    input_field = Field(input_dims, input_values)
    coarse_dims = (2,2)
    coarse_field = Field(coarse_dims)
    output_dims = (2,2)
    output_values = np.array([[1,1],[3,3]])
    output_field = Field(output_dims, output_values)

    Ctr.simple(input_field, coarse_field)
    assert array_equal(coarse_field, output_field)


def test_sum4way():
    """
    Tests the Contractinator.py 'sum4way' method.

    Parameters
    ----------
    
    Returns
    -------
    :
        Nothing, but asserts if 'sum4way' properly sums items from the Field as it should.
    """
    # test 2D array
    input_dims = (4,4)
    input_values = np.array([[1,1,1,1],[2,2,2,2],[3,3,3,3],[4,4,4,4]])
    input_field = Field(input_dims, input_values)
    coarse_dims = (2,2)
    coarse_field = Field(coarse_dims)
    output_dims = (2,2)
    output_values = np.array([[6,6],[14,14]])
    output_field = Field(output_dims, output_values)

    Ctr.sum4way(input_field, coarse_field)
    assert array_equal(coarse_field, output_field)


def test_conservative4way():
    """
    Tests the Contractinator.py 'conservative4way' method. Note: does not test weighted averaging.

    Parameters
    ----------
    
    Returns
    -------
    :
        Nothing, but asserts if 'conservative4way' properly averages items from the Field as it should.
    """
    # test 2D array
    input_dims = (4,4)
    input_values = np.array([[1,3,1,3],[1,3,1,3],[1,3,1,3],[1,3,1,3]])
    input_field = Field(input_dims, input_values)
    coarse_dims = (2,2)
    coarse_field = Field(coarse_dims)
    output_dims = (2,2)
    output_values = np.array([[2,2],[2,2]])
    output_field = Field(output_dims, output_values)

    Ctr.conservative4way(input_field, coarse_field)
    assert array_equal(coarse_field, output_field)