"""This module is an inherited class of workspace. 
    Libraries/Modules:
        Workspace\n
        numpy\n
        Field\n
"""
import numpy as np
from bin.Field import Field
from bin.Workspace import Workspace
from bin.Field import copy, isfinite, norm

class CellCenterWS(Workspace):
    """Implements a Workspace using cell-centered discretization of the grid

    Constructor:
        Args:
            grid (Grid): a Grid object specifying the geometry

    Attributes: 
        self.grid: Inputted grid
        self.flds: Dictionary of fields residing on grid        
    """

    # Return another instance of CellCenterWS
    def make_new(self, grid):
        """Creates a new workspace corresponding to a grid of half the size
        
        Args:
            grid (Grid):  The grid for the current workspace
            isFinest (bool): Whether or not this is
        """
        return CellCenterWS(grid, False)
    
    # dimensions of field (# of control volumes)
    def field_size(self):
        """Returns the 2-dimenstional size of the field __> (n, 1) for a 1-d Field
        
        """
        [xv, yv] = self.grid_size()
        return (xv-1, yv-1)

    def edges(self, dim):
        """Returns a Field containing the edge vectors

            Args:
                dim (0 or 1): Which edges will be returned (0 for i, 1 for j edges)
        
        """
        varName = "dx" + str(dim)
        if not self.exists(varName):
            self.__calc_edges(dim)

        return self.get_field(varName)

    def edge_normals(self, dim):
        """Returns a Field containing the unit normal vectors to the edges along a given dimension

            Args:
                dim (0 or 1): Which edges normals will be returned for (0 for i, 1 for j edges)
        
        """
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
        n = Field(dx.shape())
        nx = dx[:,:,1]
        ny = -dx[:,:,0]

        dn = norm(nx, ny)
        n[:,:,0] = nx/dn
        n[:,:,1] = ny/dn

        varName = "nx" + str(dim)
        self.add_field(n, varName)