
from abc import abstractmethod
from WorkspaceClass import WorkspaceClass

class Integrator(WorkspaceClass):
    
    @abstractmethod
    def __init__(self, model, input):
        pass
    
    @abstractmethod
    def step(self, workspace, state, forcing):
        pass