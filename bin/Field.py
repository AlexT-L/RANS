import numpy as np

class Field:

    def __init__(self, grid_size, stateDim=1):
        grid_size.append(stateDim)
        self.dims = grid_size
        self.vals = np.ones(grid_size, order = 'F') # set fortran ordering for f2py

    # Allow fields to be indexed like numpy arrays
    def __getitem__(self,indx): 
        return self.vals[indx]
    
    def set_init_vals(self, vals):
        self.vals = vals


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




