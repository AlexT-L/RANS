from os import stat
import numpy as np

class Field:

    def __init__(self, field_size, stateDim=1):
        self.fieldDim = field_size
        self.varDim = stateDim
        nx = field_size[0]
        ny = field_size[1]
        full_dim = (nx, ny, stateDim)
        self.vals = np.zeros(full_dim, order = 'F') # set fortran ordering for f2py

    # Allow fields to be indexed like numpy arrays
    def __getitem__(self,indx):
        x = indx[0]
        y = indx[1]
        z = 0
        if len(indx) == 3:
            z = indx[2]
        indx = (x, y, z)
        fieldSlice = Field(self.fieldDim, self.varDim)
        fieldSlice.vals = self.vals[indx]
        return fieldSlice

    # Allow fields values to be set
    def __setitem__(self,indx,value):
        x = indx[0]
        y = indx[1]
        z = 0
        if len(indx) == 3:
            z = indx[2]
        indx = (x, y, z)
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
        return "Field: (\n" + np.array_str(self.vals) + " )"

    # comparison operatos

    def __lt__(self, other):
        if isinstance(other, Field):
            other = other.vals
        return self.vals < other
    
    def __le__(self, other):
        if isinstance(other, Field):
            other = other.vals
        return self.vals <= other
    
    def __gt__(self, other):
        if isinstance(other, Field):
            other = other.vals
        return self.vals > other
    
    def __ge__(self, other):
        if isinstance(other, Field):
            other = other.vals
        return self.vals >= other
    
    def __eq__(self, other):
        if isinstance(other, Field):
            other = other.vals
        return self.vals == other
    
    def __ne__(self, other):
        if isinstance(other, Field):
            other = other.vals
        return self.vals != other


    # math operators
    
    def __add__(self, other):
        if isinstance(other, Field):
            other = other.vals

        result = self.vals + other

        output = Field(self.fieldDim, self.varDim)
        output.vals = result

        return output
    
    def __sub__(self, other):
        if isinstance(other, Field):
            other = other.vals

        result = self.vals + other

        output = Field(self.fieldDim, self.varDim)
        output.vals = result

        return output
    
    def __mul__(self, other):
        if isinstance(other, Field):
            other = other.vals

        result = self.vals * other

        output = Field(self.fieldDim, self.varDim)
        output.vals = result

        return output
    
    def __pow__(self, other):
        if isinstance(other, Field):
            other = other.vals

        result = self.vals ^ other

        output = Field(self.fieldDim, self.varDim)
        output.vals = result

        return output
    
    def __truediv__(self, other):
        if isinstance(other, Field):
            other = other.vals

        result = self.vals / other

        output = Field(self.fieldDim, self.varDim)
        output.vals = result

        return output
    
    def __floordiv__(self, other):
        if isinstance(other, Field):
            other = other.vals

        result = self.vals // other

        output = Field(self.fieldDim, self.varDim)
        output.vals = result

        return output
    
    def __mod__(self, other):
        if isinstance(other, Field):
            other = other.vals

        result = self.vals % other

        output = Field(self.fieldDim, self.varDim)
        output.vals = result

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

    