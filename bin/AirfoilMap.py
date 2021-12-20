"""This module creates am AirfoilMap object containing x,xc and vol Field objects

    Libraries/Modules:
        numpy\n
        typing\n
        Field\n
        Grid\n
        Input\n
        airfoil_map\n
    
        """
from typing import Dict
import numpy as np
from numpy.core.numeric import isclose
from numpy.lib.function_base import iterable
from bin.Field import Field
from bin.Grid import Grid
from bin.Input import Input

#class methods
from bin.utils.airfoil_map import init_from_file
from bin.utils.airfoil_map import init_from_grid
        
class AirfoilMap(Grid):
    """Creates Airfoil map object containing x,xc and vol as Field objects.

    Constructor (not intended to be implemented directly):
        Args:
            num_divisions (list):Number of cells in the x and y directions.

        Returns:
            A new AirfoilMap object. 

        Notes:
            Directly use from_file() method to perform confromal mapping
            and use from_grid() method to convert grid to coarser version 

    Attributes:
        dim_p (list): List of paramters to get from input file to dims dict.
        solv_p (list): List of paramters to get from input file to solv_param dict.
        flo_p (list): List of paramters to get from input file to flo_param dict.
        geo_p (list): List of paramters to get from input file to geo_param dict.
        in_v (list): List of paramters to get from input file to in_var dict."""
    
    def __init__(self, num_divisions):
        
        [nx, ny] = num_divisions
        self.divisions = (int(nx), int(ny))

    @classmethod
    def from_file(thisClass, num_divisions, input):
        """Initializes new AirfoilMap from datafile input.
        
        Args:
            num_divisions (list):Number of cells in the x and y directions.
            input:Dictionary containing data-file values
        Returns:
            grid (obj): new AirfoilMap object
        """
        grid = AirfoilMap(num_divisions)
        init_from_file(grid, num_divisions, input)
        return grid
    
    @classmethod
    def from_grid(thisClass, grid):
        """Initializes new AirfoilMap from existing object. The new grid will be half the size.
        
        Args:
            grid (obj): AirfoilMap object
            
        Returns:
            newGrid (obj): new AirfoilMap object
        """

        assert type(grid) is AirfoilMap

        [nx, ny] = grid.divisions
        newGrid = AirfoilMap((int(nx/2), int(ny/2)))
        init_from_grid(newGrid, grid)
        return newGrid

    def get_dims(self):
        """ Gets dimensions of grid"""
        return self.dims 
   
    def get_geometry(self):
        """Gets geometry """
        return self.geo

    def get_size(self):
        """Gets size"""
        [nx, ny] = self.divisions
        return (int(nx+1), int(ny+1))




# input=Input("rae1-s1.data")
# print("INPUT")
# grid = AirfoilMap(input)
# print("SQRT")




