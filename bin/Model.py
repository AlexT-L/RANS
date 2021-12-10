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

    # return state dimesions
    def dim(self):
        return self.dim