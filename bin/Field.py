import numpy as np


# Field classs math methods
def mean(array):
    assert(isinstance(array, Field))
    return Field(array.fieldDim, array.varDim, np.mean(array.vals))

def abs(array):
    assert(isinstance(array, Field))
    return Field(array.fieldDim, array.varDim, np.abs(array.vals))

def max(array):
    assert(isinstance(array, Field))
    return Field(array.fieldDim, array.varDim, np.max(array.vals))

def min(array):
    assert(isinstance(array, Field))
    return Field(array.fieldDim, array.varDim, np.min(array.vals))

def sum(array):
    assert(isinstance(array, Field))
    return Field(array.fieldDim, array.varDim, np.sum(array.vals))

def sqrt(array):
    assert(isinstance(array, Field))
    result = array.vals**(0.5)
    return Field(array.fieldDim, array.varDim, result)

def norm(array1, array2):
    assert(isinstance(array1, Field))
    assert(isinstance(array2, Field))
    result = (array1.vals**2 + array2.vals**2)**(0.5)
    return Field(array1.fieldDim, array1.varDim, result)

def pos_diff(array1, array2):
    assert(isinstance(array1, Field))
    assert(isinstance(array2, Field))
    diff = array1.vals - array2.vals
    result = np.max(diff, 0)
    return Field(array1.fieldDim, array1.varDim, result)




class Field:

    def __init__(self, field_size, stateDim=1, vals=None):
        if type(stateDim) is not int:
            vals = stateDim
            stateDim = 1

        self.fieldDim = field_size
        self.varDim = stateDim

        if vals is None:
            [nx, ny] = field_size
            full_dim = (nx, ny, stateDim)
            vals = np.zeros(full_dim, order = 'F') # set fortran ordering for f2py

        self.vals = vals

    # Allow fields to be indexed like numpy arrays
    def __getitem__(self,indx):
        x = indx
        y = 0
        z = 0
        dim = 1

        if type(indx) is not int:
            x = indx[0]
            y = indx[1]
        
            if len(indx) == 3:
                z = indx[2]
            
        indx = (x, y, z)

        if len(self.vals.shape) == 2:
            indx = (x, y)
        
        vals = self.vals[indx]
        shape = vals.shape

        if len(shape) == 3:
            dim = shape[2]
            shape = shape[0:2]

        fieldSlice = Field(shape, dim, vals)
        return fieldSlice

    # Allow fields values to be set
    def __setitem__(self,indx,value):
        if isinstance(value, Field):
            value = value.vals

        x = indx
        y = 0
        z = 0
        dim = 1

        if type(indx) is not int:
            x = indx[0]
            y = indx[1]
        
            if len(indx) == 3:
                z = indx[2]
            
        indx = (x, y, z)

        if len(self.vals.shape) == 2:
            indx = (x, y)
        
        self.vals[indx] = np.array(value, order = 'F')
    
    def set_val(self, new_vals):
        if np.shape(new_vals) != np.shape(self.vals):
            raise ValueError('Dimensions of field do not match expected dimensions')
        self.vals = np.array(new_vals, order = 'F')  # make new fortran ordered array  
    
    # return underlying implementation of field values
    def get_vals(self):
        return self.vals

    # size of field
    def size(self):
        vals = self.vals
        return vals.shape[0:2]

    # size of field
    def shape(self):
        vals = self.vals
        return vals.shape

    # dimension of variable
    # size of field
    def dim(self):
        vals = self.vals
        return vals.shape[2]


#############################
#       Math Methods        #
#############################

# for store methods, would be ideal if var1 and/or var2 could be individual values (not a field)

    # store the sum of var1 and var2 in self
    def store_sum(self, var1, var2):

        if isinstance(var1, Field):
            var1 = var1.vals
        if isinstance(var2, Field):
            var2 = var2.vals

        self.vals = var1 + var2
        return


    # store the difference (var1 - var2) in self
    def store_difference(self, var1, var2):

        if isinstance(var1, Field):
            var1 = var1.vals
        if isinstance(var2, Field):
            var2 = var2.vals

        self.vals = var1 - var2
        return


    # store the elementwise product of var1 and var2 in self
    def store_product(self, var1, var2):

        if isinstance(var1, Field):
            var1 = var1.vals
        if isinstance(var2, Field):
            var2 = var2.vals

        self.vals = var1 * var2
        return


    # store the elementwise quotient (var1/var2) in self
    def store_quotient(self, var1, var2):

        if isinstance(var1, Field):
            var1 = var1.vals
        if isinstance(var2, Field):
            var2 = var2.vals

        self.vals = var1 / var2
        return
    

    # elementwise copy self into copy
    def copy_to(self, copy):
        assert(isinstance(copy, Field))

        copy.vals = np.copy(self.vals)
        return


    # elementwise copy self into copy
    def copy_from(self, source):
        assert(isinstance(source, Field))

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
        [nx, ny] = self.fieldDim
        dim = self.varDim
        trans = Field((ny, nx), dim)
        trans.vals = self.vals.T
        return trans

    
    def __len__(self):
        return len(self.vals)
    

    # comparison operatos

    def __lt__(self, other):
        if isinstance(other, Field):
            other = other.vals
        result = self.vals < other
        return Field(self.fieldDim, self.varDim, result)
    
    def __le__(self, other):
        if isinstance(other, Field):
            other = other.vals
        result =  self.vals <= other
        return Field(self.fieldDim, self.varDim, result)
    
    def __gt__(self, other):
        if isinstance(other, Field):
            other = other.vals
        result = self.vals > other
        return Field(self.fieldDim, self.varDim, result)
    
    def __ge__(self, other):
        if isinstance(other, Field):
            other = other.vals
        result = self.vals >= other
        return Field(self.fieldDim, self.varDim, result)
    
    def __eq__(self, other):
        if isinstance(other, Field):
            other = other.vals
        result = self.vals == other
        return Field(self.fieldDim, self.varDim, result)
    
    def __ne__(self, other):
        if isinstance(other, Field):
            other = other.vals
        result = self.vals != other
        return Field(self.fieldDim, self.varDim, result)

    def __bool__(self):
        return bool(self.vals)


    # math operators
    
    def __add__(self, other):
        if isinstance(other, Field):
            other = other.vals

        result = self.vals + other

        output = Field(self.fieldDim, self.varDim, result)

        return output
    
    def __sub__(self, other):
        if isinstance(other, Field):
            other = other.vals

        result = self.vals + other

        output = Field(self.fieldDim, self.varDim, result)

        return output
    
    def __mul__(self, other):
        if isinstance(other, Field):
            other = other.vals

        result = self.vals * other

        output = Field(self.fieldDim, self.varDim, result)

        return output
    
    def __pow__(self, other):
        if isinstance(other, Field):
            other = other.vals

        result = self.vals ^ other

        output = Field(self.fieldDim, self.varDim, result)

        return output
    
    def __truediv__(self, other):
        if isinstance(other, Field):
            other = other.vals

        result = self.vals / other

        output = Field(self.fieldDim, self.varDim, result)

        return output
    
    def __floordiv__(self, other):
        if isinstance(other, Field):
            other = other.vals

        result = self.vals // other

        output = Field(self.fieldDim, self.varDim, result)

        return output
    
    def __mod__(self, other):
        if isinstance(other, Field):
            other = other.vals

        result = self.vals % other

        output = Field(self.fieldDim, self.varDim, result)

        return output

    def __iadd__(self, other):
        if isinstance(other, Field):
            other = other.vals

        self.vals = self.vals + other
        return self

    
    def __isub__(self, other):
        if isinstance(other, Field):
            other = other.vals

        self.vals = self.vals + other
        return self

    
    def __imul__(self, other):
        if isinstance(other, Field):
            other = other.vals

        self.vals = self.vals * other
        return self

    
    def __ipow__(self, other):
        if isinstance(other, Field):
            other = other.vals

        self.vals = self.vals ^ other
        return self


    # reverse math operations
    def __radd__(self, other):
        if isinstance(other, Field):
            other = other.vals

        result = other + self.vals

        output = Field(self.fieldDim, self.varDim, result)

        return output
    
    def __rsub__(self, other):
        if isinstance(other, Field):
            other = other.vals

        result = other - self.vals

        output = Field(self.fieldDim, self.varDim, result)

        return output
    
    def __rmul__(self, other):
        if isinstance(other, Field):
            other = other.vals

        result =  other * self.vals

        output = Field(self.fieldDim, self.varDim, result)

        return output
    
    def __rpow__(self, other):
        if isinstance(other, Field):
            other = other.vals

        result = other ^ self.vals

        output = Field(self.fieldDim, self.varDim, result)

        return output

    def __neg__(self):
        result = -self.vals
        output = Field(self.fieldDim, self.varDim, result)
        return output

    def __pos__(self):
        return self
    