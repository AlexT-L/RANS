from numpy.core.numeric import Infinity
from Model import Model
from Workspace import Workspace
from Field import Field

class NavierStokes(Model):
    
    def __init__(self, bcmodel, input):
        self.className = "NavierStokes"
        self.BCmodel = bcmodel
        self.padding = bcmodel.padding
        self.params = input # grab physical parameters
        self.dimensions = 4

        # courant number
        self.cfl_fine = input['cflf']
        self.cfl_coarse = input['cflc']
        self.cfl_lim = Infinity
        self.cfl = min([self.cfl_fine, self.cfl_lim])
        

    # initialize state
    def init_state(self, workspace, state):
        self.__check_vars(workspace)

        # copy into padded state field
        w = workspace.get_field("w", self.className)
        self.__copy_in(state, w)

        # pass off to boundary condition model to initialize
        self.BCmodel.init_state(self, workspace, w)
        
        # copy out to non-padded field
        self.__copy_out(w, state)


    # flux calculations
    # from .model_funcs import eflux_wrap,nsflux_wrap, dflux_wrap, dfluxc_wrap

    def get_flux(self, workspace, state, output, update_factor=1):
        self.__check_vars(workspace)
        
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

        # update boundary conditions
        bcmodel = self.BCmodel
        bcmodel.bc_all(self, workspace, state)

        # calculate residuals
        # eflux_wrap.eflux(self, workspace, w, dw)
        # dflux_wrap.dflux(self, workspace, w, dw, fw, rfil)
        # if self.params.kvis > 0 and False:
            # nsflux_wrap.nsflux(self, workspace, w, dw, vw, rfil)

        # copy residuals into output array
        self.__copy_out(dw, output)



    def get_safe_timestep(self, workspace, state, timestep):
        self.__check_vars(workspace)
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
        self.__check_vars(workspace)

        # copy state into padded field
        w = workspace.get_field("w", self.className)
        self.__copy_in(state, w)

        self.BCmodel.update_physics(self, workspace, state)


    # calls 'step.f' to update stability conditions
    def update_stability(self, workspace, state):
        self.__check_vars(workspace)
        
        # copy state into padded field
        w = workspace.get_field("w", self.className)
        self.__copy_in(state, w)

        self.BCmodel.update_stability(self, workspace, state)

    
    # specify upper bound on cfl
    def update_cfl_limit(self, cfl_lim):
        self.cfl_lim = cfl_lim


    # get courant number
    def get_cfl(self, workspace):

        # set courant number
        cfl = self.cfl_coarse
        if workspace.isFinest():
            cfl = self.cfl_fine
        self.cfl = min(cfl, self.cfl_lim)

        # return courant number
        return self.cfl
            

    def transfer_down(self, workspace1, workspace2):
        self.__check_vars(workspace1)
        self.__check_vars(workspace2)
        self.BCmodel.transfer_down(self, workspace1, workspace2)

    # return state dimensions
    def dim(self):
        return self.dimensions

    # copy non-padded fields into padded fields
    def __copy_in(self, field, paddedField):
        # get field size
        [leni, lenj] = field.size()
        p = self.padding

        # perform copy operation
        for i in range(0, leni):
            for j in range(0, lenj):
                paddedField[i+p,j+p] = field[i,j]

    # extract data from a padded field
    def __copy_out(self, paddedField, field):
        # get field size
        [leni, lenj] = field.size()
        p = self.padding

        # perform copy operation
        for i in range(0, leni):
            for j in range(0, lenj):
                field[i,j] = paddedField[i+p,j+p]

    # check if dictionary has been initialized
    def __check_vars(self, workspace):
        if not workspace.has_dict(self.className):
            self.__init_vars(workspace)

    # initialize class workspace fields
    def __init_vars(self, workspace):
        pad = self.padding
        [nx, ny] = workspace.field_size()
        grid_size = workspace.grid_size()
        field_size = [pad+nx+pad, pad+ny+pad]
        stateDim = self.dimensions
        className = self.className

        # initialize list of variables to add
        vars = dict()

        # add state variables stored at cell center with padding
        for stateName in ["w", "dw", "vw", "fw"]:
            vars[stateName] = [field_size, stateDim]

        # add scalar variables stored at cell center with padding
        for stateName in ["p","radI","radJ","rfl","dtl","rfli","rflj","vol","xc","rev","rlv"]:
            vars[stateName] = [field_size, 1]

        # add scalar variables stored at edges
        for stateName in ["porI","porJ"]:
            vars[stateName] = [grid_size, 1]

        workspace.init_vars(className, vars)

        # set porosity values
        bcmodel = self.BCmodel
        porI = workspace.get_field("porI", self.className)
        porJ = workspace.get_field("porI", self.className)
        pori = bcmodel.get_pori(workspace)
        porj = bcmodel.get_porj(workspace)

        # copy over porosity values
        pori.copy_to(porI)
        porj.copy_to(porJ)

        # copy over volume and centers
        VOL = workspace.get_field("vol")
        vol = workspace.get_field("vol", self.className)
        self.__copy_in(VOL, vol)

        
        XC = workspace.get_field("xc")
        xc = workspace.get_field("xc", self.className)
        self.__copy_in(XC, xc)
