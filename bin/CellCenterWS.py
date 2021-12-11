import numpy as np
from bin.Field import Field
from bin.Grid import Grid
from bin.Model import Model
from bin.Workspace import Workspace

class CellCenterWS(Workspace):

    # Return another instance of CellCenterWS
    def MakeNew(self, grid, isFinest=True):
        return CellCenterWS(grid, isFinest)
    
# Methods for returning control volume edges

    def get_size(self):
        [xv, yv] = self.grid.get_size()
        return [xv-1, yv-1]

    # control volume center coordinates
    def get_x(self, i, j):
        return self.grid.get_xc(i, j)

    def get_y(self, i, j):
        return self.grid.get_yc(i, j)

    # volume of control volume
    def volume(self, i, j):
        return self.grid.get_volume(i, j)

# For edges: side = "N", "S", "E", "W"

    # edge vector of control volume in positive i or j direction
    def edge(self, i, j, side):
        g = self.grid
        i1 = i; i1 = i; j2 = j; j2 = j
        
        if side == "N":
            i1 = i-1; j1 = j
            i2 = i  ; j2 = j
        if side == "S":
            i1 = i-1; j1 = j-1
            i2 = i  ; j2 = j-1
        if side == "E":
            i1 = i-1; j1 = j-1
            i2 = i-1; j2 = j
        if side == "W":
            i1 = i; j1 = j-1
            i2 = i; j2 = j

        x1 = g.get_x(i1, j1)
        y1 = g.get_y(i1, j1)
        x2 = g.get_x(i2, j2)
        y2 = g.get_y(i2, j2)

        dx = x2-x1
        dy = y2-y1
        
        return [dx, dy]

    # normal vector of control volume edge in positive i or j direction
    def edgeNormal(self, i, j, side):
        [dx, dy] = self.edge(i, j, side)

        return [dy, -dx]
        