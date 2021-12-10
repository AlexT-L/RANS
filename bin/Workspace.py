from abc import ABC, abstractmethod
import numpy as np
from Field import Field
from Grid import Grid
from Model import Model

class Workspace(ABC):
    
    # Constructor
    def __init__(self, grid, isFinest=True):
        self.grd = grid

        # initialize fields array with Grid fields
        self.flds = { 'Grid':{} }
        init_vals = np.zeros(self.grd.dims)
        for i in range(0,len(self.mdl.reqFields)): # loop over required fields
            self.flds.Grid[self.mdl.reqFields[i]] = Field(init_vals)

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
        if fieldName in self.flds[className]:
            print("Field already exists")
            return 1
            
        else: 
            # Add new field to workspace 
            classWorkspace = self.flds[className]
            classWorkspace[fieldName] = new_field
        
        return 0
    
    # get field method
    def get_field(self, fieldName, className='Grid'):

        # check that field exists
        if not self.exists(fieldName, className):
            print("Field doesn't exist")
            return -1
 
        else: 
            # Return field
            classWorkspace = self.flds[className]
            field = classWorkspace[fieldName]
        
        return 0

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
