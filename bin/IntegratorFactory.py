import RKInt

class IntegratorFactory:
    
    def __init__(self, integ_name, model, integ_params):
        self.integ = None
        if integ_name == 'RK4':
            self.integ = RKInt(mode, integ_params)
            
        # Check to make sure integ is picked
        
    def get(self):
        return self.integ