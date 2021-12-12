import numpy as np
import math
from math import fmod as mod
import Expandinator as expand
import Contractinator as contract
import Model, Integrator, Grid, Field, Cycle
from bin.Workspace import Workspace


class MultiGrid:

    # Constructor
    def __init__(self, workspace, model, integrator, input):
        # Parameters
        self.f_relax = input.fcoll
        self.num_cycles = 0

        # Objects
        self.cycle = Cycle(input)
        self.Model = model
        self.Integrator = integrator

        # Number of Cycles
        n_levels = self.cycle.levels
        stateDim = self.Model.dim()

        # Direct storage of variables
        self.Workspaces   = dict()
        self.W            = dict()
        self.W1st         = dict()
        self.WCorrections = dict()
        self.Residuals    = dict()
        self.Fluxes       = dict()
        self.visits       = np.zeros(n_levels, dtype=int)
        

        # set up grids
        self.Workspaces[n_levels-1] = workspace
        for l in range(n_levels-2, -1, -1):
            newGrid = Grid(self.Workspaces[l+1].get_grid())
            self.Workspaces[l] = workspace.MakeNew(newGrid, False)
        
        # initialize state variables
        for l in range(n_levels):
            grid = self.Grids[l]
            gridSize = grid.get_size()

            def newStateField():
                return Field(gridSize, stateDim)

            self.W[l]         = newStateField()
            self.W1st[l]      = newStateField()
            self.WCorrections = newStateField()
            self.Residuals[l] = newStateField()
            self.Fluxes[l]    = newStateField()
    

    # perform one iteration of the given cycle
    def performCycle(self):
        n_levels = self.cycle.levels()
        model = self.Model
        integrator = self.Integrator
        cycleIndex = self.num_cycles + 1

        # toggle for updating stability
        if self.num_cycles == 0:
            UPDATE_STABILITY = True
        else:
            UPDATE_STABILITY = (self.num_cycles % self.stabilityUpdateFrequency) == 0

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
                contract.conservative4way(self.W[prev], self.Grids[prev], w)
                contract.sum4way(self.Fluxes[prev], wr)
                wr.scale(self.f_relax)

                # If first time at this grid level, store baseline state into w1
                if self.visits[level] == 1:
                    w.copyTo(w1)

                # Check if stability needs to be updated
                if UPDATE_STABILITY:
                    model.update_stability(workspace, w)
                
                # Perform integration to get new state
                integrator.step(workspace, w, Rw)

                # Update Correction
                wc.storeDifference(w, w1)

            elif dir > 0: # go up a level
                # Allow model to transfer data to next mesh
                model.transfer_down(self.Workspaces[prev], workspace)

                # Transer correction(s) from coarser mesh(es)
                expand.bilinear4way(self.WCorrections[prev], wc)

                # Update state
                w.storeSum(wc, w)

                # Update residuals
                model.get_flux(workspace, w, Rw)
                
                # Update Correction
                wc.storeDifference(w, w1)

        # update number of cycles
        self.num_cycles += 1

                
    def residuals(self):
        dw = self.Workspaces[-1].dw
        return np.max(dw)
        
    def solution(self):
        return self.W[-1]