import numpy as np
import Expandinator as expand
import Contractinator as contract
from Grid import Grid
from Field import Field
from Cycle import Cycle
from Workspace import Workspace
from integrator import Integrator


class MultiGrid:

    # Constructor
    def __init__(self, workspace, model, integrator, input):
        # Parameters
        self.stabilityUpdateFrequency = input['ftim']
        self.wr_relax = input['fcoll']

        # counter variable
        self.num_cycles = 0

        # Objects
        self.cycle = Cycle(input)
        self.Model = model
        self.Integrator = integrator

        # Number of Cycles
        n_levels = self.cycle.levels
        stateDim = model.dim()

        # Direct storage of variables
        self.Workspaces   = [None] * n_levels
        self.W            = [None] * n_levels
        self.W1st         = [None] * n_levels
        self.WCorrections = [None] * n_levels
        self.Residuals    = [None] * n_levels
        self.Fluxes       = [None] * n_levels
        self.Volume       = [None] * n_levels
        self.visits       = np.zeros(n_levels, dtype=int)
        

        # set up grids
        self.Workspaces[-1] = workspace
        for l in range(n_levels-2, -1, -1):
            # newGrid = Grid(self.Workspaces[l+1].get_grid())
            newGrid = workspace.get_grid()
            self.Workspaces[l] = workspace.MakeNew(newGrid, False)
        
        # initialize state variables
        for l in range(n_levels):
            workspace = self.Workspaces[l]
            field_size = workspace.field_size()

            def newStateField():
                return Field(field_size, stateDim)

            self.W[l]            = newStateField()
            self.W1st[l]         = newStateField()
            self.WCorrections[l] = newStateField()
            self.Residuals[l]    = newStateField()
            self.Fluxes[l]       = newStateField()
            vol = Field(field_size, 1)
            self.Volume[l]       = vol

            for i in range(field_size[0]):
                for j in range(field_size[1]):
                    vol[i,j] = 1

        # set initial state values
        model.init_state(self.Workspaces[-1], self.W[-1])
    

    # perform one iteration of the given cycle
    def performCycle(self):
        n_levels = self.cycle.levels
        model = self.Model
        integrator = self.Integrator
        cycleIndex = self.num_cycles + 1

        # toggle for updating stability
        if self.num_cycles == 0:
            UPDATE_STABILITY = True
        else:
            UPDATE_STABILITY = (self.num_cycles % self.stabilityUpdateFrequency) == 0

    ##### first level #####
        # set pointers to working variables
        workspace = self.Workspaces[-1]
        w = self.W[-1]
        Rw = self.Fluxes[-1]
        # Update residuals
        model.get_flux(workspace, w, Rw)
                
        # Check if stability needs to be updated
        if UPDATE_STABILITY:
            model.update_stability(workspace, w)
                
        # Perform integration to get new state
        integrator.step(workspace, w, Rw)
    #####

        # subsequent levels
        level = n_levels-1
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
                grid = workspace.get_grid()
                vol = self.Volume[level]

                # Transfer state and residuals (fluxes) down to coarse mesh
                model.transfer_down(self.Workspaces[prev], workspace)
                contract.conservative4way(self.W[prev], w, vol)
                contract.sum4way(self.Fluxes[prev], wr)

                # relax transferred residuals
                wr.scale(self.wr_relax)

                # If first time at this grid level, store baseline state into w1
                if self.visits[level] == 1:
                    w.copy_to(w1)

                # Check if stability needs to be updated
                if UPDATE_STABILITY:
                    model.update_stability(workspace, w)
                
                # Perform integration to get new state
                integrator.step(workspace, w, Rw)

                # Update Correction
                wc.store_difference(w, w1)

            elif dir > 0: # go up a level
                # Transer correction(s) from coarser mesh(es)
                # expand.bilinear4way(self.WCorrections[prev], wc)

                # Update state
                w.store_sum(wc, w)

                # Update residuals
                model.get_flux(workspace, w, Rw)
                
                # Update Correction
                wc.store_difference(w, w1)

        # update number of cycles
        self.num_cycles += 1

    # copy residuals into output field
    def residuals(self, output):
        residuals = self.Residuals[-1]
        residuals.copy_to(output)
        
    # copy state into output field
    def solution(self, output):
        state = self.W[-1]
        state.copy_to(output)