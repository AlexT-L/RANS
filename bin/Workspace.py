"""This module is the ABC for workspace. 
    Libraries/Modules:
        abc\n
        numpy\n
        Field
"""
from abc import ABC, abstractmethod
import numpy as np
from bin.Field import Field, isfinite

class Workspace(ABC):
    """ Abstract base class for workspaces. 
    Constructor initializes fields array with Grid fields

    Attributes: 
        self.grid: Inputted grid
        self.flds: Dictionary of fields for grid        

    Notes: Would not be directly implemented (ABC)
    """
    def __init__(self, grid, isFinest=True):
        """
        
        """
        self.grid = grid
        self.flds = { 'Grid': {} }
        gridFields = self.flds['Grid']
        for fieldName in grid.fields:
            field = grid.fields[fieldName]
            assert(isfinite(field))
            gridFields[fieldName] = field

        self.isFinest = bool(isFinest)

    @abstractmethod
    def make_new(self, grid):
        """
        Returns another instance of a workspace
        """
        return Workspace(grid, False)

    def get_grid(self):
        """
        Returns grid object
        """
        return self.grid

    def get_dims(self):
        """
        Returns grid-level-specific geometry info
        """
        return self.grid.get_dims()

    def get_geometry(self):
        """
        Returns geometry info
        """
        return self.grid.get_geometry()

    def add_field(self, new_field, fieldName, className='Grid'):
        """
        Add field method for grid. Checks if already an instance.
        Raises an error if already in list, adds new field if not.  
        """
        if fieldName in list(self.flds[className].keys()):
            raise ValueError('Field already exists: '+fieldName)
        else: 
            classWorkspace = self.flds[className]
            classWorkspace[fieldName] = new_field
    
    # get field method
    def get_field(self, fieldName, className='Grid'):
        """
        Returns the field. Checks if field aleady exists.
        If already exists, raises error. 
        If does not exist, creats a field. 
        """
        if fieldName not in list(self.flds[className].keys()):
            raise ValueError('Field does not exist: ' + fieldName) 
        else: 
            classWorkspace = self.flds[className]
            field = classWorkspace[fieldName]
        return field
    
    def has_dict(self, className):
        """
        Check if a class dictionary exists.
        """
        return className in self.flds

    def exists(self, fieldName, className='Grid'):
        """
        Checks if a field exists in a class's dictionary.
        Exists/returns true if the it does have a dictionary,
        and if there is a fieldName in self,flds. 
        """
        if not self.has_dict(className):
            return False
        if not fieldName in self.flds[className]:
            return False
        return True

    def init_vars(self, className, vars):
        """
        Initialize class's stored fields,
        must give a dictionary of variables with the following structure:
        - keys are the variable name as a string (e.g. "w")
        - values are an array of [field_dimensions, state_dim]
            - field_dim is usually the grid size in [nx, ny] format
            - state_dim is 1 for scalars   
        Creates a class dictionary. 
        Creats fields and stores in dictionary.  
        """
        
        classDict = dict()
        self.flds[className] = classDict

        for varName in vars:
            [shape] = vars[varName]
            if np.isscalar(shape):
                classDict[varName] = Field(shape, 0)
            else:
                classDict[varName] = Field(shape)


    def is_finest(self):
        """
        Checks if finest level of mesh. 
        """
        return self.isFinest

    # Methods for getting geometric info


    @abstractmethod
    def field_size(self):
        """
        Returns dimensions of field (# of control volumes)
        """
        pass

    def grid_size(self):
        """
        Returns dimensions of grid (# of vertices)
        """
        return self.grid.get_size()

    @abstractmethod
    def edges(self, i, j, side):
        """Returns a Field containing the edge vectors

            Args:
                dim (0 or 1): Which edges will be returned (0 for i, 1 for j edges)
        
        """
        pass

    @abstractmethod
    def edge_normals(self, i, j, side):
        """Returns a Field containing the unit normal vectors to the edges along a given dimension

            Args:
                dim (0 or 1): Which edges normals will be returned for (0 for i, 1 for j edges)
        
        """
        pass
