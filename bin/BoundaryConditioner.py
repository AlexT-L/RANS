from abc import abstractmethod
from ws_class import WorkspaceClass

class BoundaryConditioner(WorkspaceClass):
    
    @abstractmethod
    def __init__(self, input):
        pass

    @abstractmethod
    def update_stability(self, model, workspace, state):
        pass

    @abstractmethod
    def update_physics(self, model, workspace, state):
        pass

    # Methods for applying boundary conditions
    @abstractmethod
    def bc_far(self, model, workspace, state):
        pass

    @abstractmethod
    def bc_wall(self, model, workspace, state):
        pass

    @abstractmethod
    def halo(self, model, workspace, state):
        pass

    @abstractmethod
    def bc_all(self, model, workspace, state):
        pass

    # Transfer information between workspaces
    @abstractmethod
    def transfer_down(self, model, workspace1, workspace2):
        pass