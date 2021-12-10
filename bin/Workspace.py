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
        dim_vals = np.zeros(self.grid.dims)
        for i in range(len(self.mdl.reqFields)): # loop over required fields
            self.flds['Grid'][self.mdl.reqFields[i]] = Field(dim_vals)

        self.isFinest = bool(isFinest)

    # Return another instance of a workspace
    @abstractmethod
    def MakeNew(self, grid, isFinest=True):
        return Workspace(grid, isFinest)

    # return grid object
    def grid(self):
        return self.grd

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

    def exists(self, fieldName, className='Grid'):

        # check that class dictionary exists
        if not className in self.flds:
            return False

        # check that field exists
        if not fieldName in self.flds[className]:
            return False

        return True

    def isFinest(self):
        return self.isFinest

    

# Methods for returning control volume edges
    @abstractmethod
    def center(self, i, j):
        pass

    @abstractmethod
    def volume(self, i, j):
        pass

    @abstractmethod
    def edge(self, i, j, side):
        pass

    @abstractmethod
    def edgeNormal(self, i, j, side):
        pass
