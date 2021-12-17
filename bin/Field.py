import numpy as np

# determine if a variable is a numpy array
def is_numpy(var):
    return type(var).__module__ is np.__name__

def is_field(var):
    return type(var).__module__ is Field.__name__

# Field classs math methods
def mean(array):
    assert is_field(array)
    return Field(array.shape(), np.mean(array.vals))

def abs(array):
    is_field(array)
    return Field(array.shape(), np.abs(array.vals))

def max(array):
    is_field(array)
    return Field(array.shape(), np.max(array.vals))

def min(array):
    is_field(array)
    return Field(array.shape(), np.min(array.vals))

def sum(array):
    is_field(array)
    return Field(array.shape(), np.sum(array.vals))

def sqrt(array):
    is_field(array)
    result = array.vals**(0.5)
    return Field(array.shape(), result)

def norm(array1, array2):
    assert is_field(array1)
    assert is_field(array2)
    result = (array1.vals**2 + array2.vals**2)**(0.5)
    return Field(array1.size(), array1.dim(), result)

def pos_diff(array1, array2):
    assert is_field(array1)
    assert is_field(array2)
    diff = array1.vals - array2.vals
    result = np.max(diff, 0)
    return Field(array1.size(), array1.dim(), result)

def isfinite(array):
    is_field(array)
    return Field(array.shape(), np.isfinite(array.vals))




class Field:

    def __init__(self, shape, vals=None):
        if is_numpy(shape):
            vals = shape
            shape = vals.shape
        
        if type(shape) is int:
            shape = (shape, 1)

        assert(len(shape) > 0 and len(shape) <= 3)
        
        if len(shape) < 3:
            if len(shape) < 2:
                shape = (shape, 1)
            shape = (shape[0], shape[1], 1)

        if vals is None:
            vals = np.zeros(shape, order = 'F') # set fortran ordering for f2py

        self.fieldShape = shape
        self.vals = vals

    # Allow fields to be indexed like numpy arrays
    def __getitem__(self,indx):
        x = indx
        y = 0
        z = 0

        if type(indx) is not int:
            x = indx[0]
            y = indx[1]
        
            if len(indx) == 3:
                z = indx[2]
            
        indx = (x, y, z)

        if len(self.vals.shape) == 2:
            indx = (x, y)
        
        vals = self.vals[indx]

        if np.isscalar(vals):
            return vals

        fieldSlice = Field(vals)
        return fieldSlice

    # Allow fields values to be set
    def __setitem__(self,indx,value):
        if is_field(value):
            value = value.vals

        x = indx
        y = 0
        z = 0

        if type(indx) is not int:
            x = indx[0]
            y = indx[1]
        
            if len(indx) == 3:
                z = indx[2]
            
        indx = (x, y, z)

        if len(self.vals.shape) == 2:
            indx = (x, y)

        if not np.isscalar(value):
            np.array(value, order = 'F')
        
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
        [nx, ny, nz] = self.fieldShape
        return (nx, ny)

    # size of field
    def shape(self):
        return self.fieldShape

    # dimension of variable
    # size of field
    def dim(self):
        [nx, ny, nz] = self.fieldShape
        return nz


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

        self.vals = var1 - var2
        return


    # store the elementwise product of var1 and var2 in self
    def store_product(self, var1, var2):

        if is_field(var1):
            var1 = var1.vals
        if is_field(var2):
            var2 = var2.vals

        self.vals = var1 * var2
        return


    # store the elementwise quotient (var1/var2) in self
    def store_quotient(self, var1, var2):

        if is_field(var1):
            var1 = var1.vals
        if is_field(var2):
            var2 = var2.vals

        self.vals = var1 / var2
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

        result = self.vals * other

        output = Field(self.shape(), result)

        return output
    
    def __pow__(self, other):
        if is_field(other):
            other = other.vals

        result = self.vals ^ other

        output = Field(self.shape(), result)

        return output
    
    def __truediv__(self, other):
        if is_field(other):
            other = other.vals

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

        self.vals = self.vals * other
        return self

    
    def __ipow__(self, other):
        if is_field(other):
            other = other.vals

        self.vals = self.vals ^ other
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

        result =  other * self.vals

        output = Field(self.shape(), result)

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
    