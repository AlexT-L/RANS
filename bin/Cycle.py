from abc import ABC
import numpy as np

class Cycle(ABC):
    """
    Description
    
    Class describing multigrid cycle.
    Contains information on cycle shape and depth.

    Attributes
    
    pattern:
        sequence of directions to perform a multigrid cycle

    levels:
        depth of multigrid cycle

    Libraries/Modules
    
    numpy
    abc

    Notes
    
    Could be expanded to include standard options like "V" and "W" cycle

    Author(s)
    
    Satya Butler, Nick Conlin, Vedin Dewan, Andy Rothstein, Alex Taylor-Lash, and Brian Wynne. \n

    """
    
    # Constructor
    def __init__(self, input):
        self.pattern = [-1,-1,1,-1,1,1]
        self.levels = 3

    # array of directions to follow
    def path(self):
        return np.copy(self.pattern)

    # depth of cycle
    def depth(self):
        return int(self.levels)