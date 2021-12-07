import numpy as np
import CycleFactory, Workspace, IntegratorFactory, ExpandinatorFactory, ContractinatorFactory


class MultiGrid:
    
    # Constructor
    def __init__(self, model, workspace, input):
        self.cycle = CycleFactory(input)
        self.workspaces = np.array([workspace])
        self.model = model
        
        n_levels = self.cycle.levels
        for l in range(1, n_levels):
            self.workspaces.append(Workspace(self.workspaces[-1]))
        
        self.integrator = IntegratorFactory(input)
        self.expandinator = ExpandinatorFactory(input)
        self.contractinator = ContractinatorFactory(input)
        self.res = 1
    
    def loop(self):
        level = 0
        for dir in self.cycle.pattern:
            level += dir
            # It don't go down
            if dir < 0: 
                # It do, it do go down
                self.contractinator.contract(self.workspace[level+1], self.workspace[level])
                # Find residuals, state field
            elif dir > 0:
                # It don't go down
                # apply corrections, get state, and get residuals
                pass
    
    def res(self):
        dw = self.workspaces[-1].dw
        return np.max(dw)
        
        