import numpy as np

class Field:

    def __init__(self, grid, stateDim=1):
        self.dims = [grid.nx, grid.ny]
        self.vals = np.ones(self.dims)

    # Allow fields to be indexed like numpy arrays
    def __getitem__(self,indx): 
        return self.vals[indx]
    
    def set_init_vals(self, vals):
        self.vals = vals

    





