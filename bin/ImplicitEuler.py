from operator import mod
from bin.Field import max, copy, mismatch_mul
from bin.Workspace import Workspace
from bin.Integrator import Integrator

class ImplicitEuler(Integrator):
    """Implicit Euler mulistage integration scheme

    Constructor:
        Args:
            model (Model): physics model
            input (Dict): Dictionary with the following items:
                mstage (int):   number of stages in the multistage integration scheme
                cdis:           flux update relaxation factor --> 0: no update, 1: full update
                cstp:           timestep relaxation factor --> 0: no timestep, 1: full step

        Returns:
            A new ImplicitEuler object

    Attributes:
        Model:          physics model
        className (str): name of class
    """


    def __init__(self, model, input):
        # set attributes
        self.Model = model
        self.className = "ImplicitEuler"
        self.numStages = input['mstage']
        self.Flux_update = input['cdis']
        self.c_step = input['cstp']

    
    def step(self, workspace, state, dwf=None):
        """Returns the local timestep such that stability is maintained.
        
        Args:
            workspace:  The Workspace object
            state:      A Field containing the current state
            forcing:    Field of values on the right hand side of the equation that "force" the ODE
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
        rfl = dt*0

        # store initial state
        wn[:] = copy(w)
        
        # account for residual "forcing" from coarser meshes
        forcing = True
        if dwf is None:
            dwf = 0
            forcing = False

        # perform implicit euler step
        for stage in range(0, self.numStages):
            # calculate new flux
            model.get_flux(workspace, w, Rw, self.Flux_update[stage])

            # offset forcing by initial residual
            if forcing and stage==0:
                print("forcing!")
                dwf -= Rw

            # add forcing
            Rw += dwf

            # get local timestep
            model.get_safe_timestep(workspace, w, dt, rfl)

            # get courant number
            cfl = model.get_cfl(workspace)
            
            # scale timestep
            c_dt = cfl*self.c_step[stage]/2.0
            dt *= c_dt*rfl

            # take step
            dw[:] = mismatch_mul(Rw, dt)

            # update state
            w[:] = wn - dw

    # check if dictionary has been initialized
    def __check_vars(self, workspace):
        if not workspace.has_dict(self.className):
            self.__init_vars(workspace)

    # initialize class workspace fields
    def __init_vars(self, workspace):
        [nx, ny] = workspace.field_size()
        stateDim = self.Model.dim()
        className = self.className

        vars = dict()
        vars["dt"] = [(nx, ny)]
        for stateName in ["wn", "Rw", "dw"]:
            vars[stateName] = [(nx, ny, stateDim)]

        workspace.init_vars(className, vars)
