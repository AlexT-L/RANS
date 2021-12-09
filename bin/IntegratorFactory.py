import DynamicIntegrator

class IntegratorFactory:
    
    def __init__(self, integ_name, model, integ_params):
        self.integ = None
        if integ_name == 'dynam':
            self.integ = DynamicIntegrator(model, integ_params)
            
        # Check to make sure integ is picked
        
    def get(self):
        return self.integ