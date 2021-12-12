from Model import Model
from Workspace import Workspace
from Input import Input

class NavierStokes(Model):
    
    def __init__(self, bcmodel, input):
        self.className = "NavierStokes"
        self.BCmodel = bcmodel
        self.padding = 2 # size of halo
        self.flo_params = input.flo_params # grab physical parameters
        

    # flux calculations
    from .model_funcs import eflux_wrap,nsflux_wrap, dflux_wrap, dfluxc_wrap

    def get_flux(self, workspace, state, output, update_factor=1):
        
        # initialize the variables we want in the workspace 
        self.__check_vars(self, workspace)

        # set rfil value
        rfil = update_factor

        # retrieve necessary workspace fields
        def get(varName):
            return workspace.get_field(varName, self.className)
        w = get("w")
        fw = get("fw")
        vw = get("vw")
        dw = get("dw")


        # copy state into padded array
        self.__copy_in(state, w)

        # calculate residuals
        self.update_viscocity(self,workspace,state)

        # copy residuals into output array
        self.__copy_out(dw, output)



    def get_safe_timestep(self, workspace, state, timestep):
        # retrieve necessary workspace fields
        def get(varName):
            return workspace.get_field(varName, self.className)
        dtl = get("dtl")
        rfl = get("rfl")

        # default is to set dt to dtl
        dt = dtl

        # export dt
        self.__copy_out(dt, timestep)


    # update rev and rlv
    def update_physics(self, workspace, state):
        self.BCmodel.update_physics(self, workspace, state)


    # calls 'step.f' to update stability conditions
    def update_stability(self, workspace, state):
        # retrieve necessary workspace fields
        def get(varName):
            return workspace.get_field(varName, self.className)
        radi = get("radi")
        radj = get("radj")
        rfl = get("rfl")
        dtl = get("dtl")
        rfli = get("rfli")
        rflj = get("rflj")


        ##### NEED TO FINISH #####
        

        # set boundary values
        self.BCmodel.update_stability(workspace, fields)

    
    def transfer_down(self, workspace1, workspace2):
        # This needs to be filled in

        self.BCmodel.transfer_down(workspace1, workspace2, fields1, fields2)

    # copy non-padded fields into padded fields
    def __copy_in(self, field, paddedField):
        # get field size
        [leni, lenj] = field.size()
        p = self.padding

        # perform copy operation
        for i in range(0, leni):
            for j in range(0, lenj):
                paddedField[i+p][j+p] = field[i][j]

    # extract data from a padded field
    def __copy_out(self, paddedField, field):
        # get field size
        [leni, lenj] = field.size()
        p = self.padding

        # perform copy operation
        for i in range(0, leni):
            for j in range(0, lenj):
                field[i][j] = paddedField[i+p][j+p]

    # check if dictionary has been initialized
    def __check_vars(self, workspace):
        if not workspace.has_dict(self.className):
            self.__init_vars(workspace)

    # initialize class workspace fields
    def __init_vars(self, workspace):
        field_size = workspace.get_size()
        stateDim = self.Model.dim()
        className = self.className

        vars = dict()
        vars["dt"] = [field_size, stateDim]
        for stateName in ["P", "vw", "fw","porI","porJ","radi","radj","rfl","dtl","rfli","rflj"]:
            vars[stateName] = [field_size, stateDim]

        workspace.init_vars(className, vars)