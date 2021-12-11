import numpy as np

class Field:

    def __init__(self, field_size, stateDim=1):
        if stateDim != 1:
            field_size.append(stateDim)
        self.dims = field_size
        self.vals = np.ones(field_size, order = 'F') # set fortran ordering for f2py

    # Allow fields to be indexed like numpy arrays
    def __getitem__(self,indx): 
        return self.vals[indx]
    
    def set_init_vals(self, vals):
        self.vals = vals
        
    def set_val(self, vals):
        if np.shape(vals) != np.shape(self.vals):
            raise ValueError('Dimensions of field do not match expected dimensions')
        self.vals = vals        
    
    def get_for_array(self):
        return np.array(self.vals, order = 'F')

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




