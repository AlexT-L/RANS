from numpy.core.numeric import Infinity
from Model import Model
from Workspace import Workspace
from Field import Field, max
from model_funcs.eflux import eflux

class NavierStokes(Model):
    
    def __init__(self, bcmodel, input):
        """Constructor
        
        Parameters
        ----------
        bcmodel:
            A BoundaryConditioner object
        input:
            Dictionary with parameter values

        Returns
        -------
        :
            A new NavierStokes object.
        """
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
        """Initializes state field.
        
        Parameters
        ----------
        workspace:
            The Workspace object
        state:
            A field containing the current state
        """
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
        """Calculates the spatial flux given the current state.
        
        Parameters
        ----------
        workspace:
            The Workspace object
        state:
            A Field containing the current state
        output:
            A Field where the flux values will be stored

        Returns
        -------
        :
            A new AirfoilMap object.
        """
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
        # eflux(self, workspace, w, dw)
        # eflux_wrap.eflux(self, workspace, w, dw)
        # dflux_wrap.dflux(self, workspace, w, dw, fw, rfil)
        # if self.params.kvis > 0 and False:
            # nsflux_wrap.nsflux(self, workspace, w, dw, vw, rfil)

        # copy residuals into output array
        self.__copy_out(dw, output)



    def get_safe_timestep(self, workspace, state, timestep):
        """Returns the local timestep such that stability is maintained.
        
        Parameters
        ----------
        workspace:
            The Workspace object
        state:
            A Field containing the current state
        timestep:
            A Field where the time steps will be stored

        """
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
        """Updates physical properties of system based on state
        
        Parameters
        ----------
        workspace:
            The Workspace object
        state:
            A Field containing the current state
        """
        self.__check_vars(workspace)

        # copy state into padded field
        w = workspace.get_field("w", self.className)
        self.__copy_in(state, w)

        self.BCmodel.update_physics(self, workspace, state)


    # calls 'step.f' to update stability conditions
    def update_stability(self, workspace, state):
        """Updates the stability parameters given the current state.
        
        Parameters
        ----------
        workspace:
            The Workspace object
        state:
            A Field containing the current state
        """
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
        """Calculates the spatial flux given the current state.
        
        Parameters
        ----------
        workspace1:
            The Workspace object for the finer level
        workspace2:
            The Workspace object for the coarser level
        """
        self.__check_vars(workspace1)
        self.__check_vars(workspace2)
        self.BCmodel.transfer_down(self, workspace1, workspace2)

    # return state dimensions
    def dim(self):
        return self.dimensions

    # copy non-padded fields into padded fields
    def __copy_in(self, field, paddedField):
        # get field size
        [nx, ny] = field.size()
        pad = self.padding

        # perform copy operation
        paddedField[pad:nx+pad, pad:ny+pad, :].copy_from(field)
        # for i in range(0, leni):
        #     for j in range(0, lenj):
        #         paddedField[i+pad,j+pad] = field[i,j]

    # extract data from a padded field
    def __copy_out(self, paddedField, field):
        # get field size
        [nx, ny] = field.size()
        pad = self.padding

        # perform copy operation
        paddedField[pad:nx+pad, pad:ny+pad, :].copy_to(field)
        # for i in range(0, nx):
        #     for j in range(0, ny):
        #         field[i,j] = paddedField[i+pad,j+pad]

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
        for stateName in ["p","radI","radJ","rfl","dtl","rfli","rflj","vol","rev","rlv"]:
            vars[stateName] = [field_size, 1]

        # xc has 2 dimensions
        vars['xc'] = [field_size, 2]

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

        # set geometric values in the halo
        bcmodel.halo_geom(self, workspace)