from abc import ABC, abstractmethod
from Input import Input

class Model(ABC):
    
    def __init__(self, input):

        # get the fields required by the model
        self.reqFields = input.mdl['fields'] 

        pass
    
    @abstractmethod
    def do_thing(self):
        pass

class NavierStokes(Model):
    pass

    def convect(self):
        return 0
            
    def viscous(self):
        return 0

    def do_thing(self):
        return 0
    