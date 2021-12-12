from abc import abstractmethod
from bin.WorkspaceClass import WorkspaceClass

class BoundaryConditioner(WorkspaceClass):
    
    @abstractmethod
    def __init__(self, input):
        pass

    @abstractmethod
    def init_state(self, model, workspace):
        pass

    @abstractmethod
    def update_stability(self, model, workspace):
        pass

    @abstractmethod
    def update_physics(self, model, workspace):
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

    # Get porosity
    def get_pori(self, workspace):
        self.__check_vars(workspace)
        return workspace.get_field("pori", self.className)

    def get_porj(self, workspace):
        self.__check_vars(workspace)
        return workspace.get_field("pori", self.className)

    # Transfer information between workspaces
    @abstractmethod
    def transfer_down(self, model, workspace1, workspace2):
        pass