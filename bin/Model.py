from abc import ABC, abstractmethod

class Model(ABC):
    
    def __init__(self, input):

        # get the fields required by the model
        self.reqFields = input.mdl.Fields 

        pass
    
    @abstractmethod
    def do_thing(self):
        pass