import numpy as np
import bin.Field as binField
from numpy.core.numeric import Infinity, NaN

class Field:
    """ Holds numeric data on a Grid. Meant to be used in a similar fashion to a numpy array. 
        Can be indexed and operators are overloaded for basic math operations.

    Constructor:
        Args:
            shape (tuple): n dimensional array of Field dimensions

        Returns:
            A new Field object

        Notes:
            Check top of Input.py file to see the contents of each of the five dictionanries 

    Attributes:
        vals (np.ndarray): numeric values of the Field

    """


    def __init__(self, shape, vals=None):
       
        if is_field(vals):
            vals = vals.vals

        if vals is not None:
            if np.isscalar(vals):
                vals = vals*np.ones(shape, order = 'F')

            if is_numpy(vals):
                if len(vals.shape) == 1 and vals.shape[0] == 1:
                    vals = vals*np.ones(shape, order = 'F')
                else:
                        vals = np.array(vals, order = 'F')
            else:
                vals = np.array(vals, order = 'F')
        else:
            vals = np.zeros(shape, order = 'F')
        assert vals is not None

        self.vals = vals

    # size of field
    def size(self):
        """2-d size of field
            
            Returns
            
            :
                The 2-d size of the field.
                This is important for fields living on a 2-d grid
            """
        shape = self.vals.shape

        if isscalar(shape):
            shape = (shape, 1)
        else:
            if len(shape) <= 1:
                shape = (shape[0], 1)

        return (shape[0], shape[1])

    # size of field
    def shape(self):
        """shape of field
            
            Returns
            
            :
                The shape of the underlying numpy array.
            """
        return self.vals.shape

    # dimension of variable
    # size of field
    def dim(self):
        """dimensions of variable
            
            Returns
            
            :
                The dimensions of the field vector living at each point in a 2-d grid
                This value is 1 for 1-d and 2-d arrays
            """
        dim = 1
        shape = self.vals.shape

        if not isscalar(shape):
            if len(shape) == 3:
                dim = shape[-1]

        return dim
    
    # return underlying implementation of field values
    def get_vals(self):
        """get the underlying numpy representation
            
            Returns
            
            :
                The underlying numpy ndarray that stores the values
            """
        return self.vals

    # return Field as given type
    def astype(self, dtype=None):
        """return a field with data stored as the given type
            
            Returns
            
            :
                A field with values stored as the given type
            """
        result = self.vals.astype(dtype)
        return Field(0, result)

    # transpose
    def T(self):
        """return a transposed Field 
            
            Returns
            
            :
                A field of the size of the transposed input
            """
        trans = Field(0, self.vals.T)
        return trans

    
    # Special methods

    # Allow fields to be indexed like numpy arrays
    def __getitem__(self,indx):
        # recursively specialize slicing
        if isinstance(indx, slice):
            return self.__class__(self[x] for x in xrange(*indx.indices(len(self))))

        # indx is now a correct int, within range(len(self))
        return self.vals[indx]


        if indx is int:
            if indx >= len(self):
                raise IndexError('End')
        
        if isscalar(self.vals):
            return self.vals

        vals = self.vals[indx]

        if np.isscalar(vals):
            return vals

        return Field(0, vals)

    # Allow fields values to be set
    def __setitem__(self,indx,value):
        if is_field(value):
            value = value.vals

        if not np.isscalar(value):
            value = np.array(value, order = 'F')
        
        self.vals[indx] = value
    
    def set_val(self, new_vals):
        """assign values to a Field of the same size
            
            Returns
            
            :
                A field with given values 
            """
        if np.shape(new_vals) != np.shape(self.vals):
            raise ValueError('Dimensions of field do not match expected dimensions')
        self.vals = np.array(new_vals, order = 'F')  # make new fortran ordered array  

    # string representation
    def __str__(self):
        return "Field: (\n" + str(self.vals) + " )"

    # length
    def __len__(self):
        return len(self.vals)

    def __iter__(self):
        return (self[i] for i in range(len(self)))
    

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
            result += self.vals[:,:,i] * other

        self.vals = result
        return self

    def __mismatch_itruediv(self, other):
        if is_field(other):
            other = other.vals
        
        k = self.shape()[2]
        result = self.vals[:,:,0] / other
        for i in range(1,k):
            result += self.vals[:,:,i] / other

        self.vals = result
        return self
    


##### Field Class Math and Logic Operators #####
##### These work similar to numpy functions #####

# determine if a variable is a numpy array
def is_numpy(var):
    return type(var).__module__ is np.__name__

def is_field(var):
    isField = type(var).__module__ is Field.__name__
    isBinField = type(var).__module__ is binField.__name__
    return isField or isBinField

def array_equal(array1, array2):
    if is_field(array1):
        array1 = array1.vals
    if is_field(array2):
        array2 = array2.vals
    return np.array_equal(array1, array2)

def copy(array):
    assert is_field(array)
    copy = np.copy(array.vals)
    return Field(0, copy)

# Field class math methods
def mean(array, axis=None):
    assert is_field(array)
    result = np.mean(array.vals, axis)

    if np.isscalar(result):
        return result

    return Field(0, result)

def abs(array):
    assert is_field(array)
    return Field(0, np.abs(array.vals))

def max(array, axis=None):
    assert is_field(array)
    
    result = np.max(array.vals, axis)

    if np.isscalar(result):
        return result

    return Field(0, result)

def min(array, axis=None):
    assert is_field(array)
    
    result = np.min(array.vals, axis)

    if np.isscalar(result):
        return result

    return Field(0, result)

def sum(array):
    assert is_field(array)
    
    result = np.sum(array.vals)

    if np.isscalar(result):
        return result

    return Field(0, result)

def sqrt(array):
    assert is_field(array)
    result = array.vals**(0.5)
    return Field(0, result)

def square(array):
    assert is_field(array)
    result = np.square(array.vals)
    return Field(0, result)

def pow(array, power):
    assert is_field(array)
    result = np.power(array.vals, power)
    return Field(0, result)

def norm(array1, array2):
    assert is_field(array1)
    assert is_field(array2)
    result = (array1.vals**2 + array2.vals**2)**(0.5)
    return Field(0, result)

def pos_diff(array1, array2):
    assert is_field(array1)
    assert is_field(array2)
    diff = array1.vals - array2.vals
    zeros = np.zeros(array1.shape())
    result = np.maximum(diff, zeros)
    return Field(0, result)

def isfinite(array):
    assert is_field(array)
    
    result = np.isfinite(array.vals)

    if len(result) == 1:
        return result

    return Field(0, result)

def isscalar(array):
    if is_field(array):
        array = array.vals

    if np.isscalar(array):
        return True
    else:
        if is_field(array):
            array = array.vals
        if is_numpy(array):
            shape = array.shape
            if shape is int:
                return True
            else:
                return len(shape) <= 1
        return False

def minimum(array1, array2):
    assert is_field(array1)
    assert is_field(array2)
    result = np.minimum(array1.vals, array2.vals)
    return Field(0, result)

def maximum(array1, array2):
    assert is_field(array1)
    assert is_field(array2)
    result = np.maximum(array1.vals, array2.vals)
    return Field(0, result)


def mismatch_mul(self, other):
    if is_field(other):
        other = other.vals
    if is_field(self):
        self = self.vals
    
    if isscalar(other) or isscalar(self):
        result = self*other
        if np.isscalar(result):
            return result
        return Field(0, result)

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
        for i in range(k):
            result.vals[:,:,i] = self[:,:,i] * other
    else:
        for i in range(k):
            result.vals[:,i] = self[:,i] * other

    assert not np.isscalar(result)
    assert is_field(result)

    return result

def mismatch_truediv(self, other):
    if is_field(other):
        other = other.vals
    if is_field(self):
        self = self.vals

    if isscalar(other) or isscalar(self):
        result = self*other
        if np.isscalar(result):
            return result
        return Field(0, result)

    result = 0

    if len(other.shape) == 3:
        result = Field(other.shape)
        k = other.shape[2]
        for i in range(k):
            result.vals[:,:,i] = other[:,:,i] / self.vals
    else:
        k = self.shape[2]
        result = Field(self.shape)
        for i in range(k):
            result.vals[:,:,i] = self.vals[:,:,i] / other

    assert not np.isscalar(result)
    assert is_field(result)

    return result

Infinity = Infinity
NaN = NaN