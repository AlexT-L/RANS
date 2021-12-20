import numpy as np
from bin.Field import Field
from bin.Grid import Grid
from bin.Model import Model
from bin.Workspace import Workspace
from bin.Field import copy, isfinite

class CellCenterWS(Workspace):

    # Return another instance of CellCenterWS
    def MakeNew(self, grid, isFinest=True):
        return CellCenterWS(grid, isFinest)
    
    # dimensions of field (# of control volumes)
    def field_size(self):
        [xv, yv] = self.grid_size()
        return (xv-1, yv-1)

    # volume of control volume
    def volume(self, i, j):
        vol = self.get_field('vol')
        return vol[i,j]

    # return x field
    def x(self):
        return self.get_field('x')

    # return xc field
    def xc(self):
        return self.get_field('xc')

    def edges(self, dim):
        varName = "dx" + str(dim)
        if not self.exists(varName):
            self.__calc_edges(dim)

        return self.get_field(varName)

    def edge_normals(self, dim):
        varName = "nx" + str(dim)
        if not self.exists(varName):
            self.__calc_normals(dim)
        
        return self.get_field(varName)
            
        
    def __calc_edges(self, dim):
        varName = "dx" + str(dim)
        X = self.get_field('x')
        x = copy(X)
        dx = 0

        assert(isfinite(x))
        
        [nx, ny] = self.field_size()

        # i edges
        if dim == 0:
            dx = x[0:nx+1, 1:ny+1, :] - x[0:nx+1, 0:ny, :]

        # j edges
        if dim == 1:
            dx = x[1:nx+1, 0:ny+1, :] - x[0:nx, 0:ny+1, :]

        assert(max(dx) != 0)
        assert(isfinite(dx))

        self.add_field(dx, varName)

    def __calc_normals(self, dim):
        dx = self.edges(dim)
        nx = Field(dx.shape())
        nx[:,:,0] = dx[:,:,1]
        nx[:,:,1] = -dx[:,:,0]

        varName = "nx" + str(dim)
        self.add_field(dx, varName)