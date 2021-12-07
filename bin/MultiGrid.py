class MultiGrid:
    
    # Constructor
    def __init__(self, model, workspace, input):
        self.cycle = 0
        self.workspace = workspace
        self.model = model
        self.integrator = 0
        self.expandinstor = 0
        self.contractinator = 0
    
    def loop(self):
        
        for dir in self.cycle.pattern:
            # It don't go down
            if dir < 0: 
                # It do, it to go down