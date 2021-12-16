from typing import Dict
import numpy as np
from numpy.core.numeric import isclose
from numpy.lib.function_base import iterable
from Field import Field
from Grid import Grid
from Input import Input

#class methods
from utils.airfoil_map import init_from_file
from utils.airfoil_map import init_from_grid
        
class AirfoilMap(Grid):
    
    def __init__(self, num_divisions):
        [nx, ny] = num_divisions
        self.divisions = (int(nx), int(ny))

    @classmethod
    def from_file(thisClass, num_divisions, input):
        grid = AirfoilMap(num_divisions)
        init_from_file(grid, num_divisions, input)
        return grid
    
    @classmethod
    def from_grid(thisClass, grid):
        assert(isinstance(grid, AirfoilMap))

        [nx, ny] = grid.divisions
        newGrid = AirfoilMap((int(nx/2), int(ny/2)))
        init_from_grid(newGrid, grid)
        return newGrid


    def get_dims(self):
        return self.dims 
   
    def get_geometry(self):
        return self.geo

    def get_size(self):
        [nx, ny] = self.divisions
        return (int(nx+1), int(ny+1))




# input=Input("rae1-s1.data")
# print("INPUT")
# grid = AirfoilMap(input)
# print("SQRT")




