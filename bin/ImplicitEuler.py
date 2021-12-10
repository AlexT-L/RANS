import Integrator

class ImplicitEuler(Integrator):
    # Constructor
    def __init__(self, model, input):
        
        # set attributes
        self.Model = model
        self.className = "ImplicitEuler"
        self.numStages = 000000000000000


    
    def step(self, workspace, state, forcing):
        model = self.model

        # make sure necessary variables exist in workspace
        self.checkVars(workspace)


        w = state
        wn = workspace.get_field("wn", self.className)
        Rw = workspace.get_field("Rw", self.className)


        # perform implicit euler step
        for stage in range(1, self.numStages):
            # calculate new flux
            model.getFlux(workspace, w, Rw)

            # add forcing
            Rw.add(forcing)

            # set timestep
            dt.storeProduct(kn, dtl)

            # take step
            dw.storeProduct(Rw, dt)

            # update state
            w.storeDifference(wn, dw)


    def updateStability(self):
        pass