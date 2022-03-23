"""
Description

Tests the Field object to see if it works as expected

Libraries/Modules

-pytest \n
-Field

Notes

Runs the following tests:\n
1. Checks that a 2D Field can be created  \n
2. Checks that a 3D field can be created \n
3. Tests that we can set a whole field \n
4. Tests that we can set an individual elements \n
5. Tests + operation \n
6. Tests - operation \n
7. Tests * operation \n
8. Tests / operation \n

Author(s)

Andy Rothstein \n

"""
from numpy.core.numeric import NaN
import pytest
from bin.Field import Field
import numpy as np

from bin.Field import isfinite, mean, copy, array_equal

def rng(dim):
    limit = np.max(dim)
    return np.random.randint(0, limit, dim)

def test_isfinite():
    """
    Asserts that we can create a 8x8 field
    
    """
    # Use dimensions for 8x8
    dims = (8,8)
    field = Field.create(dims)
    
    # Assert correct shape
    assert isfinite(field)    


def test_constructor_1d():
    """
    Asserts that we can create a 8x1 field
    
    """
    # Use dimensions for 8x8
    dims = (8)
    field = Field.create(dims, 0)
    
    # Assert correct shape
    assert array_equal(Field.size(field), (dims,1))
    

def test_constructor_2d():
    """
    Asserts that we can create a 8x8 field
    
    """
    # Use dimensions for 8x8
    dims = (8,8)
    field = Field.create(dims)
    
    # Assert correct shape
    assert array_equal(Field.size(field), dims)
    
def test_constructor_3d():
    '''
    Asserts we can 8x8 field with state dimension of 4
    '''
    # Use dimensions for 8x8
    stateDim = 4
    dims = (8,8,stateDim)
    field = Field.create(dims)
    
    # Assert correct shape
    assert array_equal(field.shape, dims)
    
def test_constructor_wrap_2arg():
    '''
    Asserts we wrap a numpy array with state dimension of 4
    '''
    # Use dimensions for 8x8
    stateDim =  4
    dims = (8,8,stateDim)
    test = rng(dims)
    field = Field.create(dims, test)
    
    # Assert correct shape
    assert array_equal(field, test)
    
def test_constructor_wrap_1arg():
    '''
    Asserts we wrap a numpy array with state dimension of 4
    '''
    # Use dimensions for 8x8
    stateDim =  4
    dims = (8,8,stateDim)
    test = rng(dims)
    field = Field.create(0, test)
    
    # Assert correct shape
    assert array_equal(field, test)

def test_set_item():
    '''
    Asserts that we can change a single element in a field
    '''
    # Create array of 0's
    dims = (8,8)
    zeros = np.zeros(dims)
    
    # Create field initially all 0's
    field = Field.create(dims, zeros)
    
    # Set index to 1 in zeros
    indx = 3
    indy = 7
    zeros[indx, indy] = 1
    
    # Set index to 1 in field
    field[indx, indy] = 1
    
    assert array_equal(zeros, field)
    
    
def test_add_func():
    '''
    Tests that the add function works out
    '''
    # Make 2 random arrays
    dims = (8,8)
    rand1 = rng(dims)
    rand2 = rng(dims)
    
    # Make fields
    field1 = Field.create(dims, rand1)
    field2 = Field.create(dims, rand2)
    
    # Add arrays
    newfield = field1 + field2
    baseline = rand1 + rand2
    
    assert array_equal(newfield, baseline)
    
    
def test_difference_func():
    '''
    Tests that the difference function works out
    '''
    # Make 2 random arrays
    dims = (8,8)
    rand1 = rng(dims)
    rand2 = rng(dims)
    
    # Make fields
    field1 = Field.create(dims, rand1)
    field2 = Field.create(dims, rand2)
    
    # Add arrays
    newfield = field1 - field2
    baseline = rand1 - rand2
    
    assert array_equal(newfield, baseline)

def test_product_func():
    '''
    Tests that the product function works out
    '''
    # Make 2 random arrays
    dims = (8,8)
    rand1 = rng(dims)
    rand2 = rng(dims)
    
    # Make fields
    field1 = Field.create(dims, rand1)
    field2 = Field.create(dims, rand2)
    
    # Add arrays
    newfield = field1 * field2
    baseline = rand1 * rand2
    
    assert array_equal(newfield.astype(int), baseline.astype(int))

def test_quotient_func():
    '''
    Tests that the add function works out
    '''
    # Make 2 random arrays
    dims = (8,8)
    rand1 = rng(dims)
    rand2 = np.ones(dims)*2
    
    # Make fields
    field1 = Field.create(dims, rand1)
    field2 = Field.create(dims, rand2)
    
    # Add arrays
    newfield = field1 / field2
    baseline = rand1 / rand2
    
    assert array_equal(newfield.astype(int), baseline.astype(int))

def test_dimensional_mean():
    '''
    Tests that the add function works out
    '''
    
    # Make 2 random arrays
    dims = (8,8)
    field = Field.create(dims)
    compare = Field.create(8,0)
    for i in range(np.max(dims)):
        field[i,:] = i
        compare[i] = i
    
    # find maximum y values along x:
    y_max = mean(field, (1))
    
    assert array_equal(y_max, compare)

def test_copy():
    '''
    Tests we change the values of an array using the copy() method
    '''
    # Make a numpy array and a field
    dims = (8,8)
    twosnp = 2*np.ones(dims)
    ones = Field.create(dims, 1)
    twos = Field.create(dims, 2)

    # copy values from ones into twos
    pointer = ones
    pointer[:] = copy(twos)
    
    assert array_equal(ones, twosnp)

def test_not_equals():
    '''
    Tests that != operator works
    '''
    # Make a numpy array and a field
    dims = (8,8)
    twosnp = 2*np.ones(dims)
    ones = Field.create(dims, 1)
    twos = Field.create(dims, 2)

    assert (ones != twos).all()
    assert (ones != twosnp).all()
    assert (ones != 0).all()
    assert (twos == twosnp).all()
