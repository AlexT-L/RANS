from os import stat
import numpy as np

class Field:

    def __init__(self, field_size, stateDim=1):
        self.dims = field_size
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
        return self.vals[indx]

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
        [nx, ny, nz] = self.shape()

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if isinstance(var1, Field):
                        k1 = np.min((k, var1.dim()-1))

                        if isinstance(var2, Field):
                            k2 = np.min((k, var2.dim()-1))
                            self[i,j,k] = var1[i,j,k1] + var2[i,j,k2]
                        else:
                            self[i,j,k] = var1[i,j,k1] + var2
                    else:
                        if isinstance(var2, Field):
                            k2 = np.min((k, var2.dim()-1))
                            self[i,j,k] = var1 + var2[i,j,k2]
                        else:
                            self[i,j,k] = var1 + var2

    # store the difference (var1 - var2) in self
    def store_difference(self, var1, var2):
        [nx, ny, nz] = self.shape()

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if isinstance(var1, Field):
                        k1 = np.min((k, var1.dim()-1))

                        if isinstance(var2, Field):
                            k2 = np.min((k, var2.dim()-1))
                            self[i,j,k] = var1[i,j,k1] - var2[i,j,k2]
                        else:
                            self[i,j,k] = var1[i,j,k1] - var2
                    else:
                        if isinstance(var2, Field):
                            k2 = np.min((k, var2.dim()-1))
                            self[i,j,k] = var1 - var2[i,j,k2]
                        else:
                            self[i,j,k] = var1 - var2

    # store the elementwise product of var1 and var2 in self
    def store_product(self, var1, var2):
        [nx, ny, nz] = self.shape()

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if isinstance(var1, Field):
                        k1 = np.min((k, var1.dim()-1))

                        if isinstance(var2, Field):
                            k2 = np.min((k, var2.dim()-1))
                            self[i,j,k] = var1[i,j,k1] * var2[i,j,k2]
                        else:
                            self[i,j,k] = var1[i,j,k1] * var2
                    else:
                        if isinstance(var2, Field):
                            k2 = np.min((k, var2.dim()-1))
                            self[i,j,k] = var1 * var2[i,j,k2]
                        else:
                            self[i,j,k] = var1 * var2

    # store the elementwise quotient (var1/var2) in self
    def store_quotient(self, var1, var2):
        [nx, ny, nz] = self.shape()

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if isinstance(var1, Field):
                        k1 = np.min((k, var1.dim()-1))

                        if isinstance(var2, Field):
                            k2 = np.min((k, var2.dim()-1))
                            self[i,j,k] = var1[i,j,k1] / var2[i,j,k2]
                        else:
                            self[i,j,k] = var1[i,j,k1] / var2
                    else:
                        if isinstance(var2, Field):
                            k2 = np.min((k, var2.dim()-1))
                            self[i,j,k] = var1 / var2[i,j,k2]
                        else:
                            self[i,j,k] = var1 / var2

    # elementwise copy self into copy
    def copy_to(self, copy):
        [nx, ny, nz] = self.shape()

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    copy[i,j,k] = self[i,j,k]

    # elementwise copy self into copy
    def copy_from(self, source):
        [nx, ny, nz] = self.shape()

        TWO_D = False
        if isinstance(source, np.ndarray):
            if len(source.shape) == 2:
                TWO_D = True

        for i in range(nx):
            for j in range(ny):
                if TWO_D:
                    self[i,j] = source[i,j]
                else:
                    for k in range(nz):
                        self[i,j,k] = source[i,j,k]


    # elementwise multiply self by k (could be field or scalar)
    def scale(self, k):
        self.store_product(self, k)




