from abc import ABC, abstractmethod

from bin.Field import Infinity

class Model(ABC):
    """
    Model Class

    Description
    
    A generic physics model for performing a multigrid cycle. 

    Libraries/Modules
    
    Field

    Notes
    
    NA

    Author(s)
    
    Satya Butler, Nick Conlin, Vedin Dewan, Andy Rothstein, Alex Taylor-Lash, and Brian Wynne. \n

    """
    @abstractmethod
    def __init__(self, bcmodel, input):
        self.dimensions = 0

    @abstractmethod
    def init_state(self, workspace, state):
        pass
    
    @abstractmethod
    def get_flux(self, workspace, state, output, update_factor=1):
        pass

    def update_cfl_limit(self, cfl_lim=Infinity):
        self.cfl_lim = cfl_lim

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
        return self.dimensions