import numpy as np
import Expandinator as expand
import Contractinator as contract
import CycleFactory, Model, Integrator,  ModelFactory, IntegratorFactory, Grid, Field
from bin.Workspace import Workspace


class MultiGrid:

    # Constructor
    def __init__(self, grid, model, integrator, input):
        # Parameters
        self.f_relax = input.fcoll

        # Objects
        self.cycle = CycleFactory(input)
        self.Model = model
        self.Integrator = integrator

        # Number of Cycles
        n_levels = self.cycle.levels

        # Direct storage of variables
        self.Workspaces = {}
        self.W = {}
        self.W1st = {}
        self.WCorrections = {}
        self.Residuals = {}
        self.Fluxes = {}
        self.visits = np.zeros(n_levels, dtype=int)


        n_levels = self.cycle.levels()
        stateDim = self.Model.dim()

        # set up grids
        self.Workspaces[n_levels] = Workspace(grid, True)
        for l in range(n_levels-1, 1):
            newGrid = Grid(self.Workspaces[l+1].getGrid())
            self.Workspaces[l] = Workspace(newGrid, False)
        
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
        n_levels = self.cycle.levels()
        model = self.Model
        integrator = self.Integrator

        level = n_levels
        for dir in self.cycle.pattern:
            level += dir
            prev = level-dir
            self.visits[level] += 1

            # set pointers to working variables
            workspace = self.Workspaces[level]
            w = self.W[level]
            w1 = self.W1st[level]
            wc = self.WCorrections[level]
            wr = self.Residuals[level]
            Rw = self.Fluxes[level]

            if dir < 0: # go down a level
                # Transfer state and residuals (fluxes) down to coarse mesh
                contract.sum4way(self.Fluxes[prev], wr)
                wr.scale(self.f_relax)
                contract.conservative4way(self.W[prev], self.Grids[prev], w)

                # If first time at this grid level, store baseline state into w1
                if self.visits[level] == 1:
                    w.copyTo(w1)
                
                # Find residuals at this grid level and calculate "forcing term"
                model.getFlux(workspace, w, Rw)
                Rw.storeDifference(wr, Rw)

                # Perform integration to get new state
                integrator.step(workspace, w, Rw)

                # Update Correction
                wc.storeDifference(w, w1)

            elif dir > 0: # go up a level
                # Transer correction(s) from coarser mesh(es)
                self.expandinator.standard4way(self.WCorrections[prev], wc)

                # Update state
                w.add(wc)

                # Update residuals
                model.getFlux(workspace, w, Rw)
                
                # Update Correction
                wc.storeDifference(w, w1)
                
    def res(self):
        n_levels = self.cycle.levels()
        dw = self.Workspaces[n_levels].dw
        return np.max(dw)
        
        