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
import pytest
from .. import Field
import numpy as np

def test_constructor_2d():
    """
    Asserts that we can create a 8x8 field
    
    """
    # Use dimensions for 8x8
    dims = (8,8)
    field = Field(dims)
    
    # Assert correct shape
    assert field.shape() == dims
    
def test_constructor_3d():
    '''
    Asserts we can 8x8 field with state dimension of 4
    '''
    # Use dimensions for 8x8
    dims = (8,8)
    stateDim =  4
    field = Field(dims, stateDim=stateDim)
    
    # Assert correct shape
    assert field.shape() == (dims[0], dims[1], stateDim)


def test_set_field():
    ''''
    Asserts that we can give Field an array and it will set all the values to that array
    '''
    # Get random array
    dims = (8,8)
    rand = np.random.rand(dims[0], dims[1])
    
    # Create field object
    field = Field(dims)
    
    # Set values
    field.set_val(rand)
    
    assert rand == field.get_vals()
    
def test_set_item():
    '''
    Asserts that we can change a single element in a field
    '''
    # Create array of 0's
    dims = (8,8)
    zeros = np.zeros(dims)
    
    # Create field initially all 0's
    field = Field(dims)
    field.set_val(zeros)
    
    # Set index to 1 in zeros
    indx = 3
    indy = 7
    zeros[indx, indy] = 1
    
    # Set index to 1 in field
    field[indx, indy] = 1
    
    assert zeros == field.get_vals()
    
    
def test_add_func():
    '''
    Tests that the add function works out
    '''
    # Make 2 random arrays
    dims = (8,8)
    rand1 = np.random.rand(dims[0], dims[1])
    rand2 = np.random.rand(dims[0], dims[1])
    
    # Make fields
    field1 = Field(dims)
    field1.set_val(rand1)
    field2 = Field(dims)
    field2.set_val(rand2)
    field_sum = Field(dims)
    
    # Add arrays
    field_sum.store_sum(field1, field2)
    
    assert field_sum.get_vals() == (rand1+rand2)
    
    
def test_difference_func():
    '''
    Tests that the difference function works out
    '''
    # Make 2 random arrays
    dims = (8,8)
    rand1 = np.random.rand(dims[0], dims[1])
    rand2 = np.random.rand(dims[0], dims[1])
    
    # Make fields
    field1 = Field(dims)
    field1.set_val(rand1)
    field2 = Field(dims)
    field2.set_val(rand2)
    field_sum = Field(dims)
    
    # Subtract arrays
    field_sum.store_difference(field1, field2)
    
    assert field_sum.get_vals() == (rand1-rand2)

def test_product_func():
    '''
    Tests that the product function works out
    '''
    # Make 2 random arrays
    dims = (8,8)
    rand1 = np.random.rand(dims[0], dims[1])
    rand2 = np.random.rand(dims[0], dims[1])
    
    # Make fields
    field1 = Field(dims)
    field1.set_val(rand1)
    field2 = Field(dims)
    field2.set_val(rand2)
    field_sum = Field(dims)
    
    # Multiply arrays
    field_sum.store_product(field1, field2)
    
    assert field_sum.get_vals() == (rand1*rand2)

def test_quotient_func():
    '''
    Tests that the add function works out
    '''
    # Make 2 random arrays
    dims = (8,8)
    rand1 = np.random.rand(dims[0], dims[1])
    rand2 = np.random.rand(dims[0], dims[1])
    
    # Make fields
    field1 = Field(dims)
    field1.set_val(rand1)
    field2 = Field(dims)
    field2.set_val(rand2)
    field_sum = Field(dims)
    
    # Divide arrays
    field_sum.store_quotient(field1, field2)
    
    assert field_sum.get_vals() == (rand1/rand2)