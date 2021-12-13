from abc import ABC

class Cycle(ABC):
    
    # Constructor
    def __init__(self):
        self.grids = []
        self.pattern = None