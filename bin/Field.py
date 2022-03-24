import numpy as np
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
        vals (np.array): numeric values of the Field

    """


    def create(shape, vals=None):
       
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

        return vals

    # size of field
    def size(self):
        """2-d size of field
            
            Returns
            
            :
                The 2-d size of the field.
                This is important for fields living on a 2-d grid
            """
        shape = self.shape

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
        return self.shape

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
        shape = self.shape

        if not isscalar(shape):
            if len(shape) == 3:
                dim = shape[-1]

        return dim
    
    # return underlying implementation of field values
    def get_vals(self):
        """get the underlying numpy representation
            
            Returns
            
            :
                The underlying numpy array that stores the values
            """
        return self

    # return Field as given type
    def astype(self, dtype=None):
        """return a field with data stored as the given type
            
            Returns
            
            :
                A field with values stored as the given type
            """
        result = self.astype(dtype)
        return np.array(0, result)

    # transpose
    def T(self):
        """return a transposed Field 
            
            Returns
            
            :
                A field of the size of the transposed input
            """
        trans = np.array(0, self.T)
        return trans

    
    # Special methods

    # Allow fields to be indexed like numpy arrays
    def __getitem__(self, indx):
        # array = super(Field, self).__getitem__(indx)
        
        # if isscalar(array):
        #     return array
        
        # return np.array(array)


        if indx is int:
            if indx >= len(self):
                raise IndexError('End')
        
        if isscalar(self):
            return self

        vals = self[indx]

        if np.isscalar(vals):
            return vals

        return np.array(0, vals)

    # Allow fields values to be set
    def __setitem__(self,indx,value):
        if is_field(value):
            value = value

        if not np.isscalar(value):
            value = np.array(value, order = 'F')
        
        self[indx] = value
    
    def set_val(self, new_vals):
        """assign values to a Field of the same size
            
            Returns
            
            :
                A field with given values 
            """
        if np.shape(new_vals) != np.shape(self):
            raise ValueError('Dimensions of field do not match expected dimensions')
        self = np.array(new_vals, order = 'F')  # make new fortran ordered array  

    # string representation
    def __str__(self):
        return "Field: (\n" + str(self) + " )"

    # length
    def __len__(self):
        return len(self)

    def __iter__(self):
        return (self[i] for i in range(len(self)))
    

    # comparison operatos
    def __lt__(self, other):
        if is_field(other):
            other = other
        result = self < other
        return np.array(self.shape(), result)
    
    def __le__(self, other):
        if is_field(other):
            other = other
        result =  self <= other
        return np.array(self.shape(), result)
    
    def __gt__(self, other):
        if is_field(other):
            other = other
        result = self > other
        return np.array(self.shape(), result)
    
    def __ge__(self, other):
        if is_field(other):
            other = other
        result = self >= other
        return np.array(self.shape(), result)
    
    def __eq__(self, other):
        if is_field(other):
            other = other
        result = self == other
        return np.array(self.shape(), result)
    
    def __ne__(self, other):
        if is_field(other):
            other = other
        result = self != other
        return np.array(self.shape(), result)

    def __bool__(self):
        return bool(self.all())


    # math operators
    def __add__(self, other):
        if is_field(other):
            other = other

        result = self + other

        output = np.array(self.shape(), result)

        return output
    
    def __sub__(self, other):
        if is_field(other):
            other = other

        result = self - other

        output = np.array(self.shape(), result)

        return output
    
    def __mul__(self, other):
        if is_field(other):
            other = other

        if not np.isscalar(other):
            mismatch = len(self.shape) - len(other.shape)
            if (mismatch):
                return mismatch_mul(self, other)

        result = self * other

        output = np.array(self.shape(), result)
        assert is_field(output)

        return output
    
    def __pow__(self, other):
        if is_field(other):
            other = other

        result = self ** other

        output = np.array(self.shape(), result)

        return output
    
    def __truediv__(self, other):
        if is_field(other):
            other = other

        if not np.isscalar(other):
            mismatch = len(self.shape) - len(other.shape)
            if (mismatch):
                return mismatch_truediv(self, other)

        result = self / other

        output = np.array(self.shape(), result)

        return output
    
    def __floordiv__(self, other):
        if is_field(other):
            other = other

        result = self // other

        output = np.array(self.shape(), result)

        return output
    
    def __mod__(self, other):
        if is_field(other):
            other = other

        result = self % other

        output = np.array(self.shape(), result)

        return output

    def __iadd__(self, other):
        if is_field(other):
            other = other

        self = self + other
        return self

    
    def __isub__(self, other):
        if is_field(other):
            other = other

        self = self + other
        return self

    
    def __imul__(self, other):
        if is_field(other):
            other = other
        
        if not np.isscalar(other):
            mismatch = len(self.shape) - len(other.shape)
            if (mismatch != 0):
                return self.__mismatch_imul(other)

        self = self * other
        return self

    
    def __ipow__(self, other):
        if is_field(other):
            other = other

        self = self ^ other
        return self

    def __itruediv__(self, other):
        if is_field(other):
            other = other
        
        if not np.isscalar(other):
            mismatch = len(self.shape) - len(other.shape)
            if (mismatch != 0):
                return self.__mismatch_itruediv(other)

        self = self / other
        return self


    # reverse math operations
    def __radd__(self, other):
        if is_field(other):
            other = other

        result = other + self

        output = np.array(self.shape(), result)

        return output
    
    def __rsub__(self, other):
        if is_field(other):
            other = other

        result = other - self

        output = np.array(self.shape(), result)

        return output
    
    def __rmul__(self, other):
        if is_field(other):
            other = other
        
        if not np.isscalar(other):
            mismatch = len(self.shape) - len(other.shape)
            if (mismatch):
                return mismatch_mul(other, self)

        result =  other * self

        output = np.array(self.shape(), result)

        assert is_field(output)

        return output
    
    def __rpow__(self, other):
        if is_field(other):
            other = other

        result = other ^ self

        output = np.array(self.shape(), result)

        return output

    def __rtruediv__(self, other):
        if is_field(other):
            other = other

        if not np.isscalar(other):
            mismatch = len(self.shape) - len(other.shape)
            if (mismatch):
                return mismatch_truediv(other, self)

        result = other / self

        output = np.array(self.shape(), result)

        return output

    def __rfloordiv__(self, other):
        if is_field(other):
            other = other

        result = other // self

        output = np.array(self.shape(), result)

        return output

    def __neg__(self):
        result = -self
        output = np.array(self.shape(), result)
        return output

    def __pos__(self):
        return self

    def __mismatch_imul(self, other):
        if is_field(other):
            other = other
        
        k = self.shape()[2]
        result = self[:,:,0] * other
        for i in range(1,k):
            result += self[:,:,i] * other

        self = result
        return self

    def __mismatch_itruediv(self, other):
        if is_field(other):
            other = other
        
        k = self.shape()[2]
        result = self[:,:,0] / other
        for i in range(1,k):
            result += self[:,:,i] / other

        self = result
        return self
    


##### Field Class Math and Logic Operators #####
##### These work similar to numpy functions #####

# determine if a variable is a numpy array
def is_numpy(var):
    return type(var).__module__ is np.__name__

def is_field(var):
    isField = type(var).__module__ is Field.__name__
    return isField

def array_equal(array1, array2):
    if is_field(array1):
        array1 = array1
    if is_field(array2):
        array2 = array2
    return np.array_equal(array1, array2)

def copy(array):
    return np.array(array)

# Field class math methods
def mean(array, axis=None):
    return np.mean(array, axis)

def abs(array):
    return np.abs(array)

def max(array, axis=None):
    return np.max(array, axis)

def min(array, axis=None):
    return np.min(array, axis)

def sum(array):
    return np.sum(array)

def sqrt(array):
    return array**(0.5)

def square(array):
    return np.square(array)

def pow(array, power):
    return np.power(array, power)

def norm(array1, array2):
    return (array1**2 + array2**2)**(0.5)

def pos_diff(array1, array2):
    diff = array1 - array2
    zeros = np.zeros(array1.shape)
    return np.maximum(diff, zeros)

def isfinite(array):
    return np.isfinite(array).all()

def isscalar(array):
    return np.isscalar(array)

def minimum(array1, array2): 
    return np.minimum(array1, array2)

def maximum(array1, array2):
    return np.maximum(array1, array2)

def mismatch_mul(self, other):
    if isscalar(other) or isscalar(self):
        return self*other

    result = 0

    diff = len(self.shape) - len(other.shape)

    assert np.abs(diff) == 1

    if diff < 0:
        temp = self
        self = other
        other = temp

    k = self.shape[-1]
    result = np.zeros(self.shape)

    if len(self.shape) == 3:
        for i in range(k):
            result[:,:,i] = self[:,:,i] * other
    else:
        for i in range(k):
            result[:,i] = self[:,i] * other

    return result

def mismatch_truediv(self, other):
    if isscalar(other) or isscalar(self):
        return self*other

    result = 0

    if len(other.shape) == 3:
        result = np.array(other.shape)
        k = other.shape[2]
        for i in range(k):
            result[:,:,i] = other[:,:,i] / self
    else:
        k = self.shape[2]
        result = np.array(self.shape)
        for i in range(k):
            result[:,:,i] = self[:,:,i] / other

    assert not np.isscalar(result)

    return result

Infinity = Infinity
NaN = NaN