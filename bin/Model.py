from abc import ABC, abstractmethod

class Model(ABC):
    
    def __init__(self, input):
        # Save some vars
        pass
    
    @abstractmethod
    def do_thing(self):
        pass