import numpy as np
import Cycle, ModelFactory, IntegratorFactory, Expandinator, Contractinator, Grid, Field

class MultiGrid:

    # Constructor
    def __init__(self, grid, modelName, integratorName, input):
        # Parameters
        self.f_relax = input.fcoll

        # Objects
        self.cycle          = Cycle(input)
        self.expandinator   = Expandinator(input)
        self.contractinator = Contractinator(input)

        # Number of Cycles
        n_levels = self.cycle.levels

        # Direct storage of variables
        self.Grids       = [None] * n_levels
        self.Models      = [None] * n_levels
        self.Integrators = [None] * n_levels
        self.W           = [None] * n_levels
        self.W1st        = [None] * n_levels
        self.WCorr       = [None] * n_levels
        self.Res         = [None] * n_levels
        self.Fluxes      = [None] * n_levels
        self.visits      = np.zeros(n_levels, dtype=int)

        n_levels = self.cycle.levels

        # set up objects
        self.Grids[n_levels-1]  = grid
        self.Models[n_levels-1] = ModelFactory(modelName, grid, input.flo_params, False).get()

        for lev in range(n_levels-2, -1, -1):
            newGrid = Grid(self.Grids[lev+1])
            newModel = ModelFactory(modelName, newGrid, input.flo_params, True).get()
            self.Grids[lev] = newGrid
            self.Models[lev] = newModel
            self.Integrators[lev] = IntegratorFactory(integratorName, newModel, input.integrator_params).get()
        
        
        stateDim = self.Models[0].dim()
        # initialize state variables
        for lev in range(n_levels):
            grid             = self.Grids[lev]
            gs = grid.get_size()
            self.W[lev]      = Field(gs, stateDim)
            self.W1st[lev]   = Field(gs, stateDim)
            self.WCorr[lev]  = Field(gs, stateDim)
            self.Res[lev]    = Field(gs, stateDim)
            self.Fluxes[lev] = Field(gs, stateDim)
    
    
    def performCycle(self):
        level = self.cycle.levels

        for dir in self.cycle.pattern:
            level += dir
            prev = level-dir
            self.visits[level] += 1

            # set working variables
            grid = self.Grids[level]
            model = self.Models[level]
            integrator = self.Integrators[level]
            w = self.W[level]
            w1 = self.W1st[level]
            wc = self.WCorr[level]
            wr = self.Res[level]
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
        dw = self.Res[-1]
        return np.max(dw)
        
        