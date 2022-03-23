import numpy as np
import bin.Expandinator as expand
import bin.Contractinator as contract
from bin.Field import Field
from bin.Cycle import Cycle
from bin.Field import copy


class MultiGrid:
    """ Uses mulitple coarser grids to apply corrections to the values on the current grid.
        Correcting the soluion this way allows for much faster convergence to be achieved
        than if the solution were only updated with the finer grid

    Constructor:
        Args:
            workspace (Workspace): 
                The workspace corresponding to the grid on which the solution will be calculated
            model (Model):
                The physics model to be used
            integrator (Integrator):
                The integration scheme to be used
            input (Dict):
                Dictionary of parameters containing:
                    ftim: the interval at which the stability will be updated
                    fcoll: the relaxation factor on the residuals transferred from the finer mesh

        Returns:
            A new Input object containing five dicts - dims, solv_param, flo_param, geo_param and in_var 

        Notes:
            Check top of Input.py file to see the contents of each of the five dictionanries 
    """

   
    # Constructor
    def __init__(self, workspace, model, integrator, input):
        # Args:
        self.stabilityUpdateFrequency = input['ftim']
        self.wr_relax = float(input['fcoll'])

        # counter variable
        self.num_cycles = 0

        # Objects
        self.cycle = Cycle(input)
        self.Model = model
        self.Integrator = integrator

        # Number of Cycles
        n_levels = self.cycle.depth()
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
            grid = self.Workspaces[l+1].get_grid()
            newGrid = grid.from_grid(grid)
            # newGrid = workspace.get_grid()
            self.Workspaces[l] = workspace.make_new(newGrid)
        
        # initialize state variables
        for l in range(n_levels):
            workspace = self.Workspaces[l]
            [nx, ny] = workspace.field_size()
            shape = (nx, ny, stateDim)

            self.W[l]            = Field.create(shape)
            self.W1st[l]         = Field.create(shape)
            self.WCorrections[l] = Field.create(shape)
            self.Residuals[l]    = Field.create(shape)
            self.Fluxes[l]       = Field.create(shape)

        # set initial state values
        model.init_state(self.Workspaces[-1], self.W[-1])
    

    # perform one iteration of the given cycle
    def performCycle(self):
        """Performs one multi-grid cycle and calculates new state.
        
        """
        n_levels = self.cycle.depth()
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
        integrator.step(workspace, w)
        #####

        # subsequent levels
        level = n_levels-1
        for dir in self.cycle.path():
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
                vol = self.Workspaces[level-dir].get_field("vol")

                # Transfer state and residuals (fluxes) down to coarse mesh
                model.transfer_down(self.Workspaces[prev], workspace)
                contract.conservative4way(self.W[prev], w, vol)
                contract.sum4way(self.Fluxes[prev], wr)

                # relax transferred residuals
                wr *= self.wr_relax

                # If first time at this grid level, store baseline state into w1
                if self.visits[level] == 1:
                    w1[:] = copy(w)
                    self.W1st[level] = w1

                # Check if stability needs to be updated
                if UPDATE_STABILITY:
                    model.update_stability(workspace, w)
                
                # Perform integration to get new state
                integrator.step(workspace, w, wr)

                # Update Correction
                wc[:] = w - w1

            elif dir > 0: # go up a level
                # Transer correction(s) from coarser mesh(es)
                expand.bilinear4way(self.WCorrections[prev], wc)

                # Update state
                w += wc

                # Update residuals
                model.get_flux(workspace, w, Rw)
                
                # Update Correction
                wc[:] = w - w1

            else: # stay on same grid level
                # perform step
                integrator.step(workspace, w, wr)

                # Update Correction
                wc[:] = w - w1
        
        # perform one last step
        integrator.step(workspace, w)

        # update number of cycles
        self.num_cycles += 1


    # copy residuals into output field
    def residuals(self, output):
        """Copies the residual values to the output Field.
        
        Args:
            output (Field): Field that will store the residual values
        """
        residuals = self.Fluxes[-1]
        output[:] = copy(residuals)
        
    # copy state into output field
    def solution(self, output):
        """Copies the state values to the output Field.
        
        Args:
            output (Field): Field that will store the values
        """
        state = self.W[-1]
        output[:] = copy(state)