from Model import Model
from Workspace import Workspace
from model_funcs import eflux_wrap, nsflux_wrap

class NavierStokes(Model):
    
    def __init__(self, bcmodel, input):
        self.className = "NavierStokes"
        self.BCmodel = bcmodel
        self.padding = 2 # size of halo

    # flux calculations
    from .model_funcs import eflux_wrap,nsflux_wrap



    def get_flux(self, workspace, state, output, update_factor=1):
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