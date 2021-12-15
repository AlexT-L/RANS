import numpy as np
from Field import Field
from Grid import Grid
from Model import Model
from Workspace import Workspace

class CellCenterWS(Workspace):

    # Return another instance of CellCenterWS
    def MakeNew(self, grid, isFinest=True):
        return CellCenterWS(grid, isFinest)
    
    # dimensions of field (# of control volumes)
    def field_size(self):
        [xv, yv] = self.grid_size()
        return [xv-1, yv-1]

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

        grid = self.grid
        x = self.x
        y = self.y

        x1 = x[i1, j1]
        y1 = y[i1, j1]
        x2 = x[i2, j2]
        y2 = y[i2, j2]

        dx = x2-x1
        dy = y2-y1
        
        return [dx, dy]

    # normal vector of control volume edge in positive i or j direction
    def edgeNormal(self, i, j, side):
        [dx, dy] = self.edge(i, j, side)

        return [dy, -dx]
        