import Integrator, Field

class ImplicitEuler(Integrator):
    # Constructor
    def __init__(self, model, input):
        
        # set attributes
        self.Model = model
        self.className = "ImplicitEuler"
        self.numStages = 000000000000000
        self.kn = 0000000000000000
        self.Flux_update = 000000000000000000 # relaxation/update factor for flux --> 0: no update, 1: full update

    
    def step(self, workspace, state, forcing):
        model = self.model

        # make sure necessary variables exist in workspace
        self.__checkVars(workspace)


        w = state
        wn = workspace.get_field("wn", self.className)
        Rw = workspace.get_field("Rw", self.className)
        dw = workspace.get_field("Rw", self.className)
        dt = workspace.get_field("dt", self.className)

        # subtract baseline residuals from forcing
        model.get_flux(workspace, w, Rw, 1)
        forcing.storeDifference(forcing, Rw)

        # perform implicit euler step
        for stage in range(0, self.numStages-1):
            # calculate new flux
            model.get_flux(workspace, w, Rw, self.Flux_update[stage])

            # add forcing
            Rw.add(forcing)

            # set timestep
            model.get_safe_timestep(workspace, dt)
            dt.scale(self.kn[stage])

            # take step
            dw.storeProduct(Rw, dt)

            # update state
            w.storeDifference(wn, dw)

    # check if dictionary has been initialized
    def __check_vars(self, workspace):
        if ~workspace.has_dict(self.className):
            self.__init_vars(workspace)

    # initialize class workspace fields
    def __init_vars(self, workspace):
        field_size = workspace.get_size()
        stateDim = self.Model.dim()
        className = self.className

        vars = dict()
        vars["dt"] = [field_size, stateDim]
        for stateName in ["wn", "Rw", "dw"]:
            vars[stateName] = [field_size, stateDim]

        workspace.init_vars(className, vars)