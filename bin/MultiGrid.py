import numpy as np
import CycleFactory, Model, Integrator,  ModelFactory, IntegratorFactory, ExpandinatorFactory, ContractinatorFactory, Grid, Field


class MultiGrid:

    # Constructor
    def __init__(self, grid, modelName, integratorName, input):
        # Parameters
        self.f_relax = input.fcoll

        # Objects
        self.cycle = CycleFactory(input)
        self.Models = {}
        self.Integrators = {}
        self.expandinator = ExpandinatorFactory(input)
        self.contractinator = ContractinatorFactory(input)

        # Number of Cycles
        n_levels = self.cycle.levels

        # Direct storage of variables
        self.Grids = {}
        self.W = {}
        self.W1st = {}
        self.WCorrections = {}
        self.Residuals = {}
        self.Fluxes = {}
        self.visits = np.zeros(n_levels, dtype=int)

        n_levels = self.cycle.levels

        # set up objects
        self.Grids[n_levels] = grid
        self.Models[n_levels] = ModelFactory(modelName, grid, input.flo_params, coarse=False)

        for l in range(n_levels-1, 1):
            newGrid = Grid(self.Grids[l+1])
            newModel = ModelFactory(modelName, newGrid, input.flo_params, coarse=True)
            self.Grids[l] = newGrid
            self.Models[l] = newModel
            self.Integrators[l] = IntegratorFactory(integratorName, newModel, input.integrator_params)
        
        
        stateDim = self.Models[1].dim()
        # initialize state variables
        for l in range(1, n_levels):
            grid = self.Grids[l]
            self.W[l] = Field(grid, stateDim)
            self.W1st[l] = Field(grid, stateDim)
            self.WCorrections = Field(grid, stateDim)
            self.Residuals[l] = Field(grid, stateDim)
            self.Fluxes[l] = Field(grid, stateDim)

        
        self.res = 1
    
    def performCycle(self):
        level = self.cycle.levels

        for dir in self.cycle.pattern:
            level += dir
            prev = level-dir
            self.visits[level] += 1

            # set pointers to working variables
            grid = self.Grids[level]
            model = self.Models[level]
            integrator = self.Integrators[level]
            w = self.W[level]
            w1 = self.W1st[level]
            wc = self.WCorrections[level]
            wr = self.Residuals[level]
            Rw = self.Fluxes[level]

            if dir < 0: # go down a level
                # Transfer state and residuals (fluxes) down to coarse mesh
                self.contractinator.sum4way(self.Fluxes[prev], wr)
                wr.scale(self.f_relax)
                self.contractinator.conservative4way(self.W[prev], self.Grids[prev], w)

                # If first time at this grid level, store baseline state into w1
                if self.visits[level] == 1:
                    w.copyTo(w1)
                
                # Find residuals at this grid level and calculate "forcing term"
                model.getFlux(w, Rw)
                Rw.storeDifference(wr, Rw)

                # Perform integration to get new state
                integrator.step(w, Rw)

                # Update Correction
                wc.storeDifference(w, w1)

            elif dir > 0: # go up a level
                # Transer correction(s) from coarser mesh(es)
                self.expandinator.standard4way(self.WCorrections[prev], wc)

                # Update state
                w.add(wc)

                # Update residuals
                model.getFlux(w, Rw)
                
                # Update Correction
                wc.storeDifference(w, w1)
    
    def res(self):
        dw = self.Workspaces[-1].dw
        return np.max(dw)
        
        