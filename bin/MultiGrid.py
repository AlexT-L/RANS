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
    
    def loop(self):
        level = 0
        for dir in self.cycle.pattern:
            level += dir
            # It don't go down
            if dir < 0: 
                # It do, it do go down
                self.contractinator.contract(self.workspaces[level-dir], self.workspaces[level])
                self.workspaces[level].w1 = self.workspaces[level].w
                self.integrator.step(self.workspaces[level])
                wc = self.workspace[level].w - self.workspaces[level].w1
            elif dir > 0:
                # It don't go down
                self.expandinator(self.workspaces[level-dir], self.workspaces[level])
                
                
                # apply corrections, get state, and get residuals
    def res(self):
        dw = self.workspaces[-1].dw
        return np.max(dw)
        
        