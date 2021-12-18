import sys
sys.path.append("../")

import numpy as np
import bin.Field as binField
from numpy.core.numeric import isscalar

# determine if a variable is a numpy array
def is_numpy(var):
    return type(var).__module__ is np.__name__

def is_field(var):
    isField = type(var).__module__ is Field.__name__
    isBinField = type(var).__module__ is binField.__name__
    return isField or isBinField

# Field classs math methods
def mean(array, axis=None):
    assert is_field(array)
    result = np.mean(array.vals, axis)

    if np.isscalar(result):
        return result

    return Field(result)

def abs(array):
    assert is_field(array)
    return Field(np.abs(array.vals))

def max(array, axis=None):
    assert is_field(array)
    
    result = np.max(array.vals, axis)

    if np.isscalar(result):
        return result

    return Field(result)

def min(array, axis=None):
    assert is_field(array)
    
    result = np.min(array.vals, axis)

    if np.isscalar(result):
        return result

    return Field(result)

def sum(array):
    assert is_field(array)
    
    result = np.sum(array.vals)

    if np.isscalar(result):
        return result

    return Field(result)

def sqrt(array):
    assert is_field(array)
    result = array.vals**(0.5)
    return Field(result)

def square(array):
    assert is_field(array)
    result = np.square(array.vals)
    return Field(result)

def pow(array, power):
    assert is_field(array)
    result = np.power(array.vals, power)
    return Field(result)

def norm(array1, array2):
    assert is_field(array1)
    assert is_field(array2)
    result = (array1.vals**2 + array2.vals**2)**(0.5)
    return Field(result)

def pos_diff(array1, array2):
    assert is_field(array1)
    assert is_field(array2)
    diff = array1.vals - array2.vals
    zeros = np.zeros(array1.shape())
    result = np.maximum(diff, zeros)
    return Field(result)

def isfinite(array):
    assert is_field(array)
    
    result = np.isfinite(array.vals)

    if np.isscalar(result):
        return result

    return Field(result)

def minimum(array1, array2):
    assert is_field(array1)
    assert is_field(array2)
    result = np.minimum(array1.vals, array2.vals)
    return Field(result)

def maximum(array1, array2):
    assert is_field(array1)
    assert is_field(array2)
    result = np.maximum(array1.vals, array2.vals)
    return Field(result)



class Field:

    def __init__(self, shape, vals=None):
        if is_numpy(shape):
            vals = shape
            shape = vals.shape
        
        if vals is None:
            assert not np.isscalar(shape)
            vals = np.zeros(shape, order = 'F') # set fortran ordering for f2py

        if np.isscalar(vals):
            vals = np.zeros(shape, order = 'F')
        else:
            assert np.array_equal(shape, vals.shape)

        self.fieldShape = shape
        self.vals = vals

    # Allow fields to be indexed like numpy arrays
    def __getitem__(self,indx):
        # x = indx
        # y = 0
        # z = 0

        # if type(indx) is not int:
        #     x = indx[0]
        #     y = indx[1]
        
        #     if len(indx) == 3:
        #         z = indx[2]
            
        # indx = (x, y, z)

        # if len(self.vals.shape) == 2:
        #     indx = (x, y)
        if np.isscalar(self.vals):
            return self.vals

        vals = self.vals[indx]

        if np.isscalar(vals):
            return vals

        return Field(vals)

    # Allow fields values to be set
    def __setitem__(self,indx,value):
        if is_field(value):
            value = value.vals

        # x = indx
        # y = 0
        # z = 0

        # if type(indx) is not int:
        #     x = indx[0]
        #     y = indx[1]
        
        #     if len(indx) == 3:
        #         z = indx[2]
            
        # indx = (x, y, z)

        # if len(self.vals.shape) == 2:
        #     indx = (x, y)

        if not np.isscalar(value):
            value = np.array(value, order = 'F')
        
        self.vals[indx] = value
    
    def set_val(self, new_vals):
        if np.shape(new_vals) != np.shape(self.vals):
            raise ValueError('Dimensions of field do not match expected dimensions')
        self.vals = np.array(new_vals, order = 'F')  # make new fortran ordered array  
    
    # return underlying implementation of field values
    def get_vals(self):
        return self.vals

    # size of field
    def size(self):
        shape = self.fieldShape

        if np.isscalar(shape):
            shape = (shape, 1)

        return (shape[0], shape[1])

    # size of field
    def shape(self):
        return self.vals.shape

    # dimension of variable
    # size of field
    def dim(self):
        dim = 1
        shape = self.fieldShape

        if not np.isscalar(shape):
            if len(shape) == 3:
                dim = shape[2]

        return dim


#############################
#       Math Methods        #
#############################

# for store methods, would be ideal if var1 and/or var2 could be individual values (not a field)

    # store the sum of var1 and var2 in self
    def store_sum(self, var1, var2):

        if is_field(var1):
            var1 = var1.vals
        if is_field(var2):
            var2 = var2.vals

        self.vals = var1 + var2
        return


    # store the difference (var1 - var2) in self
    def store_difference(self, var1, var2):
        
        if is_field(var1):
            var1 = var1.vals
        if is_field(var2):
            var2 = var2.vals

        self.vals = var1 * var2
        return


    # store the elementwise product of var1 and var2 in self
    def store_product(self, var1, var2, axis=None):

        if is_field(var1):
            var1 = var1.vals
        if is_field(var2):
            var2 = var2.vals 

        if axis != None: 
            dim = var1.shape[axis]
            if len(var2.shape) > len(var1.shape):
                temp = var1
                var1 = var2
                var2 = temp
        else:
            self.vals = var1 * var2

        if axis == 0:
            for i in range(1,dim):
                self.vals[i] = var1[i] * var2
        if axis == 1:
            for i in range(1,dim):
                self.vals[:,i] = var1[:,i] * var2
        if axis == 2:
            for i in range(1,dim):
                self.vals[:,:,i] = var1[:,:,i] * var2

        return



    # store the elementwise quotient (var1/var2) in self
    def store_quotient(self, var1, var2, axis=None):

        if is_field(var1):
            var1 = var1.vals
        if is_field(var2):
            var2 = var2.vals

        if axis != None: 
            dim = var1.shape[axis]
            if len(var2.shape) > len(var1.shape):
                temp = var1
                var1 = var2
                var2 = temp
        else:
            self.vals = var1 * var2
        
        if axis == 0:
            for i in range(1,dim):
                self.vals[i] += var1[i] / var2
        if axis == 1:
            for i in range(1,dim):
                self.vals[:,i] += var1[:,i] / var2
        if axis == 2:
            for i in range(1,dim):
                self.vals[:,:,i] += var1[:,:,i] / var2

        return
    

    # elementwise copy self into copy
    def copy_to(self, copy):
        assert is_field(copy)

        copy.vals = np.copy(self.vals)
        return


    # elementwise copy self into copy
    def copy_from(self, source):
        assert is_field(source)

        self.vals = np.copy(source.vals)
        return


    # elementwise multiply self by k (could be field or scalar)
    def scale(self, k):
        self.store_product(self, k)


    # string representation
    def __str__(self):
        return "Field: (\n" + str(self.vals) + " )"



    # matrix operations

    def T(self):
        [nx, ny] = self.size()
        dim = self.dim()
        trans = Field(self.vals.T)
        return trans

    
    def __len__(self):
        return len(self.vals)
    

    # comparison operatos

    def __lt__(self, other):
        if is_field(other):
            other = other.vals
        result = self.vals < other
        return Field(self.shape(), result)
    
    def __le__(self, other):
        if is_field(other):
            other = other.vals
        result =  self.vals <= other
        return Field(self.shape(), result)
    
    def __gt__(self, other):
        if is_field(other):
            other = other.vals
        result = self.vals > other
        return Field(self.shape(), result)
    
    def __ge__(self, other):
        if is_field(other):
            other = other.vals
        result = self.vals >= other
        return Field(self.shape(), result)
    
    def __eq__(self, other):
        if is_field(other):
            other = other.vals
        result = self.vals == other
        return Field(self.shape(), result)
    
    def __ne__(self, other):
        if is_field(other):
            other = other.vals
        result = self.vals != other
        return Field(self.shape(), result)

    def __bool__(self):
        return bool(self.vals.any())


    # math operators
    
    def __add__(self, other):
        if is_field(other):
            other = other.vals

        result = self.vals + other

        output = Field(self.shape(), result)

        return output
    
    def __sub__(self, other):
        if is_field(other):
            other = other.vals

        result = self.vals - other

        output = Field(self.shape(), result)

        return output
    
    def __mul__(self, other):
        if is_field(other):
            other = other.vals

        if not np.isscalar(other):
            mismatch = len(self.vals.shape) - len(other.shape)
            if (mismatch):
                return mismatch_mul(self, other)

        result = self.vals * other

        output = Field(self.shape(), result)
        assert is_field(output)

        return output
    
    def __pow__(self, other):
        if is_field(other):
            other = other.vals

        result = self.vals ** other

        output = Field(self.shape(), result)

        return output
    
    def __truediv__(self, other):
        if is_field(other):
            other = other.vals

        if not np.isscalar(other):
            mismatch = len(self.vals.shape) - len(other.shape)
            if (mismatch):
                return mismatch_truediv(self, other)

        result = self.vals / other

        output = Field(self.shape(), result)

        return output
    
    def __floordiv__(self, other):
        if is_field(other):
            other = other.vals

        result = self.vals // other

        output = Field(self.shape(), result)

        return output
    
    def __mod__(self, other):
        if is_field(other):
            other = other.vals

        result = self.vals % other

        output = Field(self.shape(), result)

        return output

    def __iadd__(self, other):
        if is_field(other):
            other = other.vals

        self.vals = self.vals + other
        return self

    
    def __isub__(self, other):
        if is_field(other):
            other = other.vals

        self.vals = self.vals + other
        return self

    
    def __imul__(self, other):
        if is_field(other):
            other = other.vals
        
        if not np.isscalar(other):
            mismatch = len(self.vals.shape) - len(other.shape)
            if (mismatch != 0):
                return self.__mismatch_imul(other)

        self.vals = self.vals * other
        return self

    
    def __ipow__(self, other):
        if is_field(other):
            other = other.vals

        self.vals = self.vals ^ other
        return self

    def __itruediv__(self, other):
        if is_field(other):
            other = other.vals
        
        if not np.isscalar(other):
            mismatch = len(self.vals.shape) - len(other.shape)
            if (mismatch != 0):
                return self.__mismatch_itruediv(other)

        self.vals = self.vals / other
        return self


    # reverse math operations
    def __radd__(self, other):
        if is_field(other):
            other = other.vals

        result = other + self.vals

        output = Field(self.shape(), result)

        return output
    
    def __rsub__(self, other):
        if is_field(other):
            other = other.vals

        result = other - self.vals

        output = Field(self.shape(), result)

        return output
    
    def __rmul__(self, other):
        if is_field(other):
            other = other.vals
        
        if not np.isscalar(other):
            mismatch = len(self.vals.shape) - len(other.shape)
            if (mismatch):
                return mismatch_mul(other, self)

        result =  other * self.vals

        output = Field(self.shape(), result)

        assert is_field(output)

        return output
    
    def __rpow__(self, other):
        if is_field(other):
            other = other.vals

        result = other ^ self.vals

        output = Field(self.shape(), result)

        return output

    def __rtruediv__(self, other):
        if is_field(other):
            other = other.vals

        if not np.isscalar(other):
            mismatch = len(self.vals.shape) - len(other.shape)
            if (mismatch):
                return mismatch_truediv(other, self)

        result = other / self.vals

        output = Field(self.shape(), result)

        return output

    def __rfloordiv__(self, other):
        if is_field(other):
            other = other.vals

        result = other // self.vals

        output = Field(self.shape(), result)

        return output

    def __neg__(self):
        result = -self.vals
        output = Field(self.shape(), result)
        return output

    def __pos__(self):
        return self

    def __mismatch_imul(self, other):
        if is_field(other):
            other = other.vals
        
        k = self.shape()[2]
        result = self.vals[:,:,0] * other
        for i in range(1,k):
            result += self.vals[:,:,k] * other

        self.vals = result
        return self

    def __mismatch_itruediv(self, other):
        if is_field(other):
            other = other.vals
        
        k = self.shape()[2]
        result = self.vals[:,:,0] / other
        for i in range(1,k):
            result += self.vals[:,:,k] / other

        self.vals = result
        return self
    
def mismatch_mul(self, other):
    if is_field(other):
        other = other.vals
    if is_field(self):
        self = self.vals
    
    if np.isscalar(other) or np.isscalar(self):
        result = self*other
        return Field(result)

    result = 0

    diff = len(self.shape) - len(other.shape)

    assert np.abs(diff) == 1

    if diff < 0:
        temp = self
        self = other
        other = temp

    k = self.shape[-1]
    result = Field(self.shape)

    if len(self.shape) == 3:
        for i in range(1,k):
            result.vals[:,:,i] = self[:,:,i] * other
    else:
        for i in range(1,k):
            result.vals[:,i] = self[:,i] * other

    assert not np.isscalar(result)
    assert is_field(result)

    return result

def mismatch_truediv(self, other):
    if is_field(other):
        other = other.vals
    if is_field(self):
        self = self.vals

    if np.isscalar(other) or np.isscalar(self):
        return Field(self*other)

    result = 0

    if len(other.shape) == 3:
        result = Field(other.shape)
        k = other.shape[2]
        for i in range(1,k):
            result.vals[:,:,i] = other[:,:,i] / self.vals
    else:
        k = self.shape[2]
        result = Field(self.shape)
        for i in range(1,k):
            result.vals[:,:,i] = self.vals[:,:,i] / other

    assert not np.isscalar(result)
    assert is_field(result)

    return result
