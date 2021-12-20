
from abc import ABC, abstractmethod

class Integrator(ABC):
    """ Abstract base class, never directly instantiated
        NS_Airfoil is a child class of this ABC

        Constructor:
            Args:
                model (Model): physics model
                input: necessary input parameters
    """
    
    @abstractmethod
    def __init__(self, model, input):
        pass
    
    @abstractmethod
    def step(self, workspace, state, forcing):
        """Returns the local timestep such that stability is maintained.
        
        Args:
            workspace:  The Workspace object
            state:      A Field containing the current state
            forcing:    Field of values on the right hand side of the equation that "force" the ODE
        """
        pass