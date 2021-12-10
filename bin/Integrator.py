from abc import ABC, abstractmethod
from Input import Input

class Integrator(ABC):
    
    @abstractmethod
    def __init__(self, model, input):
        pass
    
    @abstractmethod
    def step(self, workspace, state, forcing):
        pass
    