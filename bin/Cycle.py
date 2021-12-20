from abc import ABC
import numpy as np

class Cycle(ABC):
    """ Contains information about the shape and depth of multigrid cycle

    Constructor:
        Args:
            input (dictionary) : placeholder for potential input

        Returns:
            A new Cycle object

        Notes:
            Could be expanded to include default options such as "V" and "W" cycle

    Attributes:
        pattern (np.ndarray): array with sequence of directions for a cycle
        levels (int):  depth of cycle

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