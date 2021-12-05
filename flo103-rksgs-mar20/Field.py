import numpy as np

class Field:

    # constructor
    def __init__(self,dims,pad,init_vals):

        # assign dims and vals
        self.dims = dims
        self.pad = pad

        # main array contained in a Field
        self.vals = np.array(init_vals, copy = True)

    # Allow fields to be indexed like numpy arrays
    def __getitem__(self,indx): 
        return self.vals[indx]

class VectorField(Field):
    pass

class ScalarField(Field):
    pass







