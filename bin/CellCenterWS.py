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

    # control volume center coordinates
    def center(self, i, j):

        return [xc, yc]

    # volume of control volume
    def volume(self, i, j):
        
        return [vol]

    # edge vector of control volume
    # side = "N", "S", "E", "W"
    def edge(self, i, j, side):
        
        return [dx, dy]

    # normal vector of control volume edge
    def edgeNormal(self, i, j, side):
        [dx, dy] = edge(i, j, side)

        return [nx, ny]
        