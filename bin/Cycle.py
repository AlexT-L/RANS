from abc import ABC
import numpy as np

class Cycle(ABC):
    
    # Constructor
    def __init__(self, input):
        self.grids = []
        self.pattern = [-1,-1,1,-1,1,1]
        self.levels = 3

    def path(self):
        return np.copy(self.pattern)

    def depth(self):
        return int(self.levels)