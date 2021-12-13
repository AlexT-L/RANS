import Integrator, Field
from bin.Workspace import Workspace

class ImplicitEuler(Integrator):
    # Constructor
    def __init__(self, model, input):
        
        # set attributes
        self.Model = model
        self.className = "ImplicitEuler"
        self.numStages = input.mstage
        self.Flux_update = input.cdis # relaxation/update factor for flux --> 0: no update, 1: full update
        self.c_step = input.cstp    # fraction of timestep to use

    
    def step(self, workspace, state, forcing):
        model = self.model

        # make sure necessary variables exist in workspace
        self.__checkVars(workspace)

        def get(varName):
            return workspace.get_field(varName, self.className)
        w =  get("w")
        wn = get("wn")
        Rw = get("Rw")
        dw = get("dw")
        dt = get("dt")

        # subtract baseline residuals from forcing
        model.get_flux(workspace, w, Rw, 1)
        forcing.store_difference(forcing, Rw)

        # perform implicit euler step
        for stage in range(0, self.numStages-1):
            # calculate new flux
            model.get_flux(workspace, w, Rw, self.Flux_update[stage])

            # add forcing
            Rw.store_sum(Rw, forcing)

            # get local timestep
            model.get_safe_timestep(workspace, dt)

            # get courant number
            cfl = model.get_cfl(workspace)
            
            # scale timestep
            c_dt = self.cfl*self.c_step/2.0
            dt.scale(c_dt)

            # take step
            dw.store_product(Rw, dt)

            # update state
            w.store_difference(wn, dw)

    # initialize class workspace fields
    def __init_vars(self, workspace):
        field_size = workspace.field_size()
        stateDim = self.Model.dim()
        className = self.className

        vars = dict()
        vars["dt"] = [field_size, 1]
        for stateName in ["wn", "Rw", "dw"]:
            vars[stateName] = [field_size, stateDim]

        workspace.init_vars(className, vars)
