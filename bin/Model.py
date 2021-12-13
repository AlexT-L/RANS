from abc import abstractmethod
from Field import Field
from WorkspaceClass import WorkspaceClass

class Model(WorkspaceClass):
    
    @abstractmethod
    def __init__(self, bcmodel, input):
        pass

    @abstractmethod
    def init_state(self, workspace, state):
        pass
    
    @abstractmethod
    def get_flux(self, workspace, state, output, update_factor=1):
        pass

    @abstractmethod
    def get_safe_timestep(self, workspace, state, dt):
        pass

    @abstractmethod
    def update_physics(self, workspace, state):
        pass

    @abstractmethod
    def update_stability(self, workspace, state):
        pass

    @abstractmethod
    def transfer_down(self, workspace1, workspace2):
        pass


    # return state dimesions
    def dim(self):
        return self.dim