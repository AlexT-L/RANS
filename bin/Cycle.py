from abc import ABC

class Cycle(ABC):
    
    # Constructor
    def __init__(self, input):
        self.grids = []
        self.pattern = []
        self.levels = 1