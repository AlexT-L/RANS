from abc import ABC, abstractmethod
from Input import Input

class Integrator(ABC):
    
    def __init__(self, input):

        pass
    
    @abstractmethod
    def do_thing(self):
        pass


    