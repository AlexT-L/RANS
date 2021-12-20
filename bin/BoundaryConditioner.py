from abc import ABC, abstractmethod

class BoundaryConditioner(ABC):
    
    @abstractmethod
    def __init__(self, input):
        pass

    @abstractmethod
    def update_stability(self, model, workspace, state):
        """
        updates stability parameters for time step calculations
        
        Args:
            model: instance of class inheriting from Model
            workspace: instance of Workspace class (or child)
            state (Field): current state of the system (density, momentum, energy)
    
        """
        pass

    @abstractmethod
    def update_physics(self, model, workspace, state):
        """updates physical parameters for calculation of boundary conditions
        
        Args:
            model: instance of class inheriting from Model
            workspace: instance of Workspace class (or child)
            state (Field): current state of the system (density, momentum, energy)
    
        """
        pass

    @abstractmethod
    def bc_wall(self, model, workspace, state):
        """
        apply boundary condition along the wall
        
         Args:
            model: instance of class inheriting from Model
            workspace: instance of Workspace class (or child)
            state (Field): current state of the system (density, momentum, energy)
        
        
        """
        pass

    # Methods for applying boundary conditions
    @abstractmethod
    def bc_far(self, model, workspace, state):
        """
        apply boundary condition in the far field
        
        Args:
            model: instance of class inheriting from Model
            workspace: instance of Workspace class (or child)
            state (Field): current state of the system (density, momentum, energy)
    
        """
        pass

    @abstractmethod
    def halo(self, model, workspace, state):
        """
        set the values in the halo
        
         Args:
            model: instance of class inheriting from Model
            workspace: instance of Workspace class (or child)
            state (Field): current state of the system (density, momentum, energy)
        
        
        """
        pass

    @abstractmethod
    def bc_all(self, model, workspace, state):
        """
        do wall boundaries, far field and set halo values at once
        
         Args:
            model: instance of class inheriting from Model
            workspace: instance of Workspace class (or child)
            state (Field): current state of the system (density, momentum, energy)
        
        """
        pass

    # Transfer information between workspaces
    @abstractmethod
    def transfer_down(self, model, workspace1, workspace2):
        pass