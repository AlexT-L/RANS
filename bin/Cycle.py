from abc import ABC

class Cycle(ABC):
    
    # Constructor
    def __init__(self, input):
        self.grids = []
        self.pattern = [-1,1]
        self.levels = 2