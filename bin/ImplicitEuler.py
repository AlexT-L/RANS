from Field import max
from Workspace import Workspace
from Integrator import Integrator

class ImplicitEuler(Integrator):
    # Constructor
    def __init__(self, model, input):
        """Constructor
        
        Parameters
        ----------
        model:
            Model object
        input:
            Dictionary of parameter values

        Returns
        -------
        :
            A new ImplicitEuler integrator object.
        """
        # set attributes
        self.Model = model
        self.className = "ImplicitEuler"
        self.numStages = input['mstage']
        self.Flux_update = input['cdis'] # relaxation/update factor for flux --> 0: no update, 1: full update
        self.c_step = input['cstp']    # fraction of timestep to use
        print("ImplicitEuler")

    
    def step(self, workspace, state, forcing):
        """Returns the local timestep such that stability is maintained.
        
        Parameters
        ----------
        workspace:
            The Workspace object
        state:
            A Field containing the current state
        forcing:
            Field of values on the right hand side of the equation that "force" the ODE
        """
        model = self.Model

        # make sure necessary variables exist in workspace
        self.__check_vars(workspace)

        def get(varName):
            return workspace.get_field(varName, self.className)
        w =  state
        wn = get("wn")
        Rw = get("Rw")
        dw = get("dw")
        dt = get("dt")

        # store initial state
        w.copy_to(wn)

        # subtract baseline residuals from forcing
        model.get_flux(workspace, w, Rw, 1)
        forcing -= Rw

        # perform implicit euler step
        for stage in range(0, self.numStages):
            # calculate new flux
            model.get_flux(workspace, w, Rw, self.Flux_update[stage])

            # add forcing
            Rw += forcing

            # get local timestep
            model.get_safe_timestep(workspace, w, dt)

            # get courant number
            cfl = model.get_cfl(workspace)
            print(cfl)
            
            # scale timestep
            c_dt = cfl*self.c_step[stage]/2.0
            dt *= c_dt

            # take step
            dw.store_product(Rw, dt)

            # update state
            w.store_difference(wn, dw)

    # check if dictionary has been initialized
    def __check_vars(self, workspace):
        if not workspace.has_dict(self.className):
            self.__init_vars(workspace)

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
