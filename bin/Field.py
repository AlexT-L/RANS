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
    
    def get_vals(self):
        return self.vals

    # size of field
    def size(self):
        return self.dims


#############################
#       Math Methods        #
#############################

    def storeSum(self, var1, var2):
        pass

    def storeDifference(self, var1, var2):
        pass

    def storeProduct(self, var1, var2):
        pass

    def copyTo(self, copy):
        pass

    def scale(self, k):
        pass




