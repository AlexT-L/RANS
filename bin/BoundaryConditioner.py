from abc import ABC, abstractmethod

class BoundaryConditioner(ABC):
    
    @abstractmethod
    def __init__(self, input):
        pass

# Methods for applying boundary conditions
    @abstractmethod
    def bc_far(workspace, state, fields):
        pass

    @abstractmethod
    def bc_wall(workspace, state, fields):
        pass

    @abstractmethod
    def halo(workspace, state, fields):
        pass

    @abstractmethod
    def bc_all(workspace, state, fields):
        pass

    # Get porosity
    @abstractmethod
    def get_pori(self, i, j):
        pass

    @abstractmethod
    def get_porj(self, i, j):
        pass

    # Transfer information between workspaces
    @abstractmethod
    def transfer_down(self, workspace1, workspace2, fields1, fields2):
        pass