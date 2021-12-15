from abc import ABC, abstractmethod
import numpy as np
from Field import Field
from Grid import Grid
from Model import Model

class Workspace(ABC):
    
    # Constructor
    def __init__(self, grid, isFinest=True):
        self.grid = grid

        # initialize fields array with Grid fields
        self.flds = { 'Grid': {} }
        for fieldName in grid.fields:
            self.flds[fieldName] = grid.fields[fieldName]
        # dim_vals = np.zeros(self.grid.dims)
        # for i in range(len(self.mdl.reqFields)): # loop over required fields
        #     self.flds['Grid'][self.mdl.reqFields[i]] = Field(dim_vals)

        self.isFinest = bool(isFinest)

    # Return another instance of a workspace
    @abstractmethod
    def MakeNew(self, grid, isFinest=True):
        return Workspace(grid, isFinest)

    # return grid object
    def get_grid(self):
        return self.grid

    # return geometry info
    def get_geometry(self):
        return self.grid.get_geometry()

    # add field method
    def add_field(self, new_field, fieldName, className='Grid'):

        # check if we already have it
        if fieldName in list(self.flds[className].keys()):
            raise ValueError('Field already exists: '+fieldName)
        else: 
            # Add new field to workspace 
            classWorkspace = self.flds[className]
            classWorkspace[fieldName] = new_field
    
    # get field method
    def get_field(self, fieldName, className='Grid'):

        # check that field exists
        if fieldName not in list(self.flds[className].keys()):
            raise ValueError('Field does not exist: ' + fieldName) 
        else: 
            # Return field
            classWorkspace = self.flds[className]
            field = classWorkspace[fieldName]
        return field
    
    def set_field(self, fieldName, fieldVal, className='Grid'):
        # Sets value of a field if field already exists
        if fieldName not in list(self.flds[className].keys()):
            raise ValueError('Field does not exist: ' + fieldName) 
        else: 
            # Return field
            classWorkspace = self.flds[className]
            classWorkspace[fieldName].set_val(fieldVal)

    # check if a class dictionary exists
    def has_dict(self, className):
        return className in self.flds

    # check if a field exists in a class's dictionary
    def exists(self, fieldName, className='Grid'):

        # check that class dictionary exists
        if not self.has_dict(className):
            return False

        # check that field exists
        if not fieldName in self.flds[className]:
            return False

        return True

    # initialize class's stored fields
    # must give a dictionary of variables with the following structure:
    #   - keys are the variable name as a string (e.g. "w")
    #   - values are an array of [field_dimensions, state_dim]
    #       - field_dim is usually the grid size in [nx, ny] format
    #       - state_dim is 1 for scalars
    def init_vars(self, className, vars):
        # create class dictionary
        classDict = dict()
        self.flds[className] = classDict

        # create fields and store in dictionary
        for varName in vars:
            [size, dim] = vars[varName]
            newField = Field(size, dim)
            classDict[varName] = newField


    def isFinest(self):
        return self.isFinest

    # Methods for getting geometric info

    # dimensions of field (# of control volumes)
    @abstractmethod
    def field_size(self):
        pass

    # dimensions of grid (# of vertices)
    def grid_size(self):
        return self.grid.get_size()

    @abstractmethod
    def volume(self, i, j):
        pass

    @abstractmethod
    def edge(self, i, j, side):
        pass

    @abstractmethod
    def edgeNormal(self, i, j, side):
        pass
