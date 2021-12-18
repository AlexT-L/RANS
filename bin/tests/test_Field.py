"""
Description
-----------
Tests the Field object to see if it works as expected

Libraries/Modules
-----------------
-pytest \n
-Field

Notes
-----
Runs the following tests:\n
1. Checks that a 2D Field can be created  \n
2. Checks that a 3D field can be created \n
3. Tests that we can set a whole field \n
4. Tests that we can set an individual elements \n
5. Tests store_product function \n
6. Tests store_difference function \n
7. Tests store_product function \n
8. Tests store_quotient function \n

Author(s)
---------
Andy Rothstein \n

"""
import sys
sys.path.append("../")

import pytest
from Field import Field
import numpy as np

from bin.Field import isfinite

def rng(dim):
    limit = np.max(dim)
    return np.random.randint(0, limit, dim)

def test_constructor():
    """
    Asserts that we can create a 8x1 field
    
    """
    # Use dimensions for 8x8
    dims = (8,8)
    field = Field(dims)

    # Assert correct shape
    assert type(field) is Field


def test_isfinite():
    """
    Asserts that we can create a 8x1 field
    
    """
    # Use dimensions for 8x8
    dims = (8,8)
    field = Field(dims)
    
    # Assert correct shape
    assert isfinite(field)
    

def test_constructor_1d():
    """
    Asserts that we can create a 8x1 field
    
    """
    # Use dimensions for 8x8
    dims = (8)
    field = Field(dims, 0)
    
    # Assert correct shape
    assert np.array_equal(field.size(), (dims,1))
    

def test_constructor_2d():
    """
    Asserts that we can create a 8x8 field
    
    """
    # Use dimensions for 8x8
    dims = (8,8)
    field = Field(dims)
    
    # Assert correct shape
    assert np.array_equal(field.size(), dims)
    
def test_constructor_3d():
    '''
    Asserts we can 8x8 field with state dimension of 4
    '''
    # Use dimensions for 8x8
    stateDim = 4
    dims = (8,8,stateDim)
    field = Field(dims)
    
    # Assert correct shape
    assert np.array_equal(field.shape(), dims)
    
def test_constructor_wrap_2arg():
    '''
    Asserts we wrap a numpy array with state dimension of 4
    '''
    # Use dimensions for 8x8
    stateDim =  4
    dims = (8,8,stateDim)
    test = rng(dims)
    field = Field(dims, test)
    
    # Assert correct shape
    assert np.array_equal(field.vals, test)
    
def test_constructor_wrap_1arg():
    '''
    Asserts we wrap a numpy array with state dimension of 4
    '''
    # Use dimensions for 8x8
    stateDim =  4
    dims = (8,8,stateDim)
    test = rng(dims)
    field = Field(test)
    
    # Assert correct shape
    assert np.array_equal(field.vals, test)

def test_set_item():
    '''
    Asserts that we can change a single element in a field
    '''
    # Create array of 0's
    dims = (8,8)
    zeros = np.zeros(dims)
    
    # Create field initially all 0's
    field = Field(dims, zeros)
    
    # Set index to 1 in zeros
    indx = 3
    indy = 7
    zeros[indx, indy] = 1
    
    # Set index to 1 in field
    field[indx, indy] = 1
    
    assert np.array_equal(zeros, field.get_vals())
    
    
def test_add_func():
    '''
    Tests that the add function works out
    '''
    # Make 2 random arrays
    dims = (8,8)
    rand1 = rng(dims)
    rand2 = rng(dims)
    
    # Make fields
    field1 = Field(dims, rand1)
    field2 = Field(dims, rand2)
    
    # Add arrays
    newfield = field1 + field2
    baseline = rand1 + rand2
    
    assert np.array_equal(newfield.get_vals(), baseline)
    
    
def test_difference_func():
    '''
    Tests that the difference function works out
    '''
    # Make 2 random arrays
    dims = (8,8)
    rand1 = rng(dims)
    rand2 = rng(dims)
    
    # Make fields
    field1 = Field(dims, rand1)
    field2 = Field(dims, rand2)
    
    # Add arrays
    newfield = field1 - field2
    baseline = rand1 - rand2
    
    assert np.array_equal(newfield.get_vals(), baseline)

def test_product_func():
    '''
    Tests that the product function works out
    '''
    # Make 2 random arrays
    dims = (8,8)
    rand1 = rng(dims)
    rand2 = rng(dims)
    
    # Make fields
    field1 = Field(dims, rand1)
    field2 = Field(dims, rand2)
    
    # Add arrays
    newfield = field1 * field2
    baseline = rand1 * rand2
    
    assert np.array_equal(newfield.get_vals().astype(int), baseline.astype(int))

def test_quotient_func():
    '''
    Tests that the add function works out
    '''
    # Make 2 random arrays
    dims = (8,8)
    rand1 = rng(dims)
    rand2 = np.ones(dims)*2
    
    # Make fields
    field1 = Field(dims, rand1)
    field2 = Field(dims, rand2)
    
    # Add arrays
    newfield = field1 / field2
    baseline = rand1 / rand2
    
    assert np.array_equal(newfield.get_vals().astype(int), baseline.astype(int))
