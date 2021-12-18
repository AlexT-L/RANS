import numpy as np
from Field import Field
from Grid import Grid
from Model import Model
from Workspace import Workspace
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

# For edges: side = "n", "s", "e", "w"

    def edges(self, dim):
        varName = "dx" + str(dim)
        if not self.exists(varName):
            self.__calc_edges(dim)

        return self.get_field(varName)

    # edge vector of control volume in positive i or j direction
    def edge(self, i, j, side):
        g = self.grid
        i1 = i; i1 = i; j2 = j; j2 = j
        
        if side == "n":
            i1 = i  ; j1 = j+1
            i2 = i+1; j2 = j+1
        if side == "s":
            i1 = i  ; j1 = j
            i2 = i+1; j2 = j
        if side == "e":
            i1 = i; j1 = j
            i2 = i; j2 = j+1
        if side == "w":
            i1 = i+1; j1 = j
            i2 = i+1; j2 = j+1

        X = self.get_field('x')
        x = X[:,:,0]
        y = X[:,:,1]

        x1 = x[i1, j1]
        y1 = y[i1, j1]
        x2 = x[i2, j2]
        y2 = y[i2, j2]

        dx = x2-x1
        dy = y2-y1
        
        return [dx, dy]

    # normal vector of control volume edge in positive i or j direction
    def edge_normal(self, i, j, dim):
        [dx, dy] = self.edge(i, j, dim)

        return [dy, -dx]

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