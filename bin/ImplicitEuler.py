import Integrator, Field

class ImplicitEuler(Integrator):
    # Constructor
    def __init__(self, model, input):
        
        # set attributes
        self.Model = model
        self.className = "ImplicitEuler"
        self.numStages = 000000000000000
        self.kn = 0000000000000000
        self.cdis = 000000000000000000

    
    def step(self, workspace, state, forcing):
        model = self.model

        # make sure necessary variables exist in workspace
        self.__checkVars(workspace)


        w = state
        wn = workspace.get_field("wn", self.className)
        Rw = workspace.get_field("Rw", self.className)
        dw = workspace.get_field("Rw", self.className)


        # perform implicit euler step
        for stage in range(0, self.numStages-1):
            # calculate new flux
            model.getFlux(workspace, w, Rw, self.cdis[stage])

            # add forcing
            Rw.add(forcing)

            # set timestep
            dt.storeProduct(self.kn[stage], dtl)

            # take step
            dw.storeProduct(Rw, dt)

            # update state
            w.storeDifference(wn, dw)


    def updateStability(self):
        pass

    def getCdis(self):
        return self.cdis[0]

    def __checkVars(self, workspace):
        grid = workspace.getGrid()
        gridSize = grid.getSize()
        stateDim = self.Model.dim()
        className = self.className

        def exist(var):
            return workspace.exist(var, className)

        for stateName in enumerate(["wn", "Rw", "dw"]):
            if ~exist(stateField):
                stateField = Field(gridSize, stateDim)
                workspace.add_field(stateField, stateName, className)