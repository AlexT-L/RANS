import numpy as np

class Field:

    def __init__(self, field_size, stateDim=1):
        if stateDim != 1:
            field_size.append(stateDim)
        self.dims = field_size
        self.vals = np.NaN*np.ones(field_size, order = 'F') # set fortran ordering for f2py

    # Allow fields to be indexed like numpy arrays
    def __getitem__(self,indx): 
        return self.vals[indx]
    
    def set_val(self, new_vals):
        if np.shape(new_vals) != np.shape(self.vals):
            raise ValueError('Dimensions of field do not match expected dimensions')
        self.vals = np.array(new_vals,order = 'F')  # make new fortran ordered array  
    
    # return underlying implementation of field values
    def get_vals(self):
        return self.vals

    # size of field
    def size(self):
        return self.dims


#############################
#       Math Methods        #
#############################

# for store methods, would be ideal if var1 and/or var2 could be individual values (not a field)

    # store the sum of var1 and var2 in self
    def store_sum(self, var1, var2):
        pass

    # store the difference (var1 - var2) in self
    def store_difference(self, var1, var2):
        pass

    # store the elementwise product of var1 and var2 in self
    def store_product(self, var1, var2):
        pass

    # store the elementwise quotient (var1/var2) in self
    def store_quotient(self, var1, var2):
        pass

    # elementwise copy self into copy
    def copy_to(self, copy):
        pass

    # elementwise multiply self by k (could be field or scalar)
    def scale(self, k):
        self.storeProduct(self, k)




