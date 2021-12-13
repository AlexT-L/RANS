import NavierStokes

class ModelFactory:
    
    def __init__(self, modelName, grid, flo_params, coarse):
        self.model = None
        if modelName == 'Navier':
            self.model = NavierStokes(grid, flo_params, coarse)
        
        # Raise error if model not selected
    
    def get(self):
        return self.model