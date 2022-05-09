from numpy.core.numeric import Infinity
import numpy as np
from bin.Model import Model
from bin.Workspace import Workspace
from bin.Field import Field, max, min, isfinite, pos_diff, copy, minimum
from bin.model_funcs.eflux import eflux
from bin.model_funcs.dflux import dflux
from bin.model_funcs.dfluxc import dfluxc
from bin.model_funcs.nsflux import nsflux

UPDATE_FORTRAN_DATA = True

# Validation
from bin.model_funcs.BaldwinLomax import turbulent_viscosity
from bin.model_funcs.BoundaryThickness import boundary_thickness
if UPDATE_FORTRAN_DATA:
    from bin.model_funcs.fortran_versions.eflux_wrap import eflux as eflux_fortran
    from bin.model_funcs.fortran_versions.dflux_wrap import dflux as dflux_fortran
    from bin.model_funcs.fortran_versions.dfluxc_wrap import dfluxc as dfluxc_fortran
    from bin.model_funcs.fortran_versions.nsflux_wrap import nsflux as nsflux_fortran
    from bin.model_funcs.fortran_versions.turb2_wrap import turb_BL as turb2
    from bin.model_funcs.fortran_versions.delt_wrap import thickness as thickness_fortran
    from bin.model_funcs.fortran_versions.viscf_wrap import viscosity as viscosity_fortran
    from bin.model_funcs.fortran_versions.bcwall_wrap import bc_wall as bcwall_fortran
    from bin.model_funcs.fortran_versions.bcfar_wrap import bc_far as bcfar_fortran
    from bin.model_funcs.fortran_versions.halo_wrap import halo as halo_fortran
    from bin.model_funcs.fortran_versions.stability_wrap import stability as stability_fortran

from NS_Airfoil import IGNORE_NAN

class NavierStokes(Model):
    """Physics model for fluid flow based on the Reynolds Averaged Navier Stokes (RANS) equations 
    for use in a multigrid scheme. The state w is composed of Fields with density, x-momentum, y-momentum and energy. 
    Based on a finite volume formulation with ability to compute fluxes, update scheme stable timestep, 
    and update eddy viscocities. Contains terms for convective, artificial dissipative, and viscous fluxes. 

    Constructor:
        Args:
            bcmodel (BoundaryConditioner): instance of BoundaryConditioner class
            input (dictionary): dictionary with parameter values

        Returns:
            A new NavierStokes object


    Attributes:
      className (str): name of class for acessing it's dictionaries in the workspace
      BCmodel (BoundaryConditioner): boundary condition model, instance of BoundaryConditioner
      padding (int): outter padding for boundary condition implementation
      params (dictionary): physics parameters from the input
      dimensions (int): number of states (4)
      cfl_fine (np.array): courant number on fine mesh
      cfl_coarse (np.array): courant number on coarse mesh 
      cfl_lim (float): upper limit on courant number
      cfl (float): minimum cfl between fine and coarse grids

    Note:
        See report for more details on the physics"""
    
    def __init__(self, bcmodel, input):
        self.className = "NavierStokes"
        self.BCmodel = bcmodel
        self.padding = bcmodel.padding
        self.params = input # grab physical parameters
        self.dimensions = 4

        # courant number
        self.cfl_fine = abs(input['cflf'])
        self.cfl_coarse = abs(input['cflc'])
        self.cfl_lim = Infinity
        self.cfl = minimum(self.cfl_fine, self.cfl_lim)
        

    # initialize state
    def init_state(self, workspace, state):
        """Finds max number of columns in a row in the input file.
        
        Args:
            workspace (Workspace): containins relevant fields to initialize
            state (Field): has the state variables
        """
        self.__check_vars(workspace)

        # copy into padded state field
        w = workspace.get_field("w", self.className)
        self.__copy_in(state, w)

        # pass off to private method
        self.__init_state(workspace, w)
        
        # copy out to non-padded field
        self.__copy_out(w, state)
        assert(isfinite(state))


    # flux calculations
    def get_flux(self, workspace, state, output, update_factor=1):
        """Calculates the spatial flux given the current state.
        
        Args:
            workspace (Workspace): contains the relevant fields
            state (Field): the current state
            output (Field): where the flux values will be stored
            """
        assert(isfinite(state)) or IGNORE_NAN
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

        # update pressure
        self.__update_pressure(workspace, w)
        assert isfinite(get("p")) or IGNORE_NAN
        
        # update boundary conditions
        bcmodel = self.BCmodel
        bcmodel.bc_all(self, workspace, w)
        
        assert isfinite(w) or IGNORE_NAN

        # calculate residuals
        eflux(self, workspace, w, dw)
        assert isfinite(dw) or IGNORE_NAN

        # if workspace.is_finest():
        #     dflux_fortran(self, workspace, w, dw, rfil)
        #     assert isfinite(dw) or IGNORE_NAN
        # else:
        #     dfluxc_fortran(self, workspace, w, dw, rfil)
        #     assert isfinite(dw) or IGNORE_NAN
            
        # if UPDATE_FORTRAN_DATA:
        #     nsflux_fortran(self, workspace, w, dw, rfil)
        # else: nsflux(self, workspace, w, dw, rfil)

        if IGNORE_NAN:
            np.nan_to_num(dw,False,0.5)

        # copy residuals into output array
        self.__copy_out(dw, output)
        assert(isfinite(output)) or IGNORE_NAN



    def get_safe_timestep(self, workspace, state, timestep, rfl_out=None):
        """Returns the local timestep such that stability is maintained.
        
        Args:
            workspace (Workspace): contains the relevant fields
            state (Field): the current state
            timestep (Field): where the time steps will be stored

        """
        assert(isfinite(state)) or IGNORE_NAN
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
        assert(isfinite(timestep)) or IGNORE_NAN

        # export rfl
        if rfl is not None:
            self.__copy_out(rfl, rfl_out)
        

    # update ev and lv
    def update_physics(self, workspace, state):
        """Updates physical properties of system based on state
        
        Args:
            workspace (Workspace): contains the relevant fields
            state (Field): the current state
        """
        assert(isfinite(state)) or IGNORE_NAN
        self.__check_vars(workspace)

        # copy state into padded field
        w = workspace.get_field("w", self.className)
        self.__copy_in(state, w)

        # update pressure
        self.__update_pressure(workspace, w)

        self.BCmodel.update_physics(self, workspace, w)


    # calls python implementation 'step.f' to update stability conditions
    def update_stability(self, workspace, state):
        """Updates the stability parameters given the current state.
        
        Args:
            workspace (Workspace): contains the relevant fields
            state (Field): the current state
        """
        assert(isfinite(state)) or IGNORE_NAN
        self.__check_vars(workspace)

        # copy state into padded field
        w = workspace.get_field("w", self.className)
        self.__copy_in(state, w)

        # update pressure
        self.__update_pressure(workspace, w)
        
        self.BCmodel.update_stability(self, workspace, w)

    # get courant number
    def get_cfl(self, workspace):

        # set courant number
        cfl = self.cfl_coarse
        if workspace.is_finest():
            cfl = self.cfl_fine
        self.cfl = minimum(cfl, self.cfl_lim)

        # return courant number
        return self.cfl
            

    def transfer_down(self, workspace1, workspace2):
        """ Move workspace1 on fine mesh to workspace2 on coarse mesh
        
        Args:
            workspace1 (Workspace): The Workspace object for the finer level
            workspace2 (Workspace): The Workspace object for the coarser level
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
        [nx, ny] = Field.size(field)
        pad = self.padding

        # perform copy operation
        paddedField[pad:nx+pad, pad:ny+pad] = copy(field)
        
    # extract data from a padded field
    def __copy_out(self, paddedField, field):
        # get field size
        [nx, ny] = Field.size(field)
        pad = self.padding

        # perform copy operation
        if ny <= 1:
            field[:] = copy(paddedField[pad:nx+pad])
        elif nx <= 1:
            field[:] = copy(paddedField[:, pad:ny+pad])
        else:
            field[:] = copy(paddedField[pad:nx+pad, pad:ny+pad])
        
    # check if dictionary has been initialized
    def __check_vars(self, workspace):
        if not workspace.has_dict(self.className):
            self.__init_vars(workspace)

    # initialize class workspace fields
    def __init_vars(self, workspace):
        pad = self.padding
        [nx, ny] = workspace.field_size()
        [nxp, nyp] = [pad+nx+pad, pad+ny+pad]
        stateDim = self.dimensions
        className = self.className

        # initialize list of variables to add
        vars = dict()

        # add state variables stored at cell center with padding
        for stateName in ["w", "dw", "vw", "fw"]:
            shape = [nxp, nyp, stateDim]
            vars[stateName] = [shape]

        # add scalar variables stored at cell center with padding
        for stateName in ["p","radI","radJ","rfl","dtl","rfli","rflj","vol",'ev','lv']:
            shape = (nxp, nyp)
            vars[stateName] = [shape]

        # xc has 2 dimensions
        vars['xc'] = [(nxp, nyp, 2)]

        # add scalar variables stored at edges
        vars["porI"] = [(nx+1, ny)]
        vars["porJ"] = [(nx, ny+1)]
        
        # viscous flux
        vars["q"] = [(nxp,nyp,3,2)]
        vars["u"] = [(nxp,nyp,3)]
        vars["fs"] = [(nxp,nyp,3)]

        workspace.init_vars(className, vars)

        # set porosity values
        bcmodel = self.BCmodel
        porI = workspace.get_field("porI", self.className)
        porJ = workspace.get_field("porI", self.className)
        pori = bcmodel.get_pori(workspace)
        porj = bcmodel.get_porj(workspace)

        # copy over porosity values
        porI[:] = copy(pori)
        porJ[:] = copy(porj)

        # copy over volume and centers
        VOL = workspace.get_field("vol")
        vol = workspace.get_field("vol", self.className)
        self.__copy_in(VOL, vol)

        assert min(VOL) > 0
        assert min(vol[2:nx+2, 2:ny+2]) > 0

        
        XC = workspace.get_field("xc")
        xc = workspace.get_field("xc", self.className)
        self.__copy_in(XC, xc)

        # set geometric values in the halo
        bcmodel.halo_geom(self, workspace)

    # initialize state
    def __init_state(self, workspace, state):
        # get pressure
        p = workspace.get_field("p", self.className)

        # set initial values
        rho0 = self.params['rho0']
        u0 = self.params['u0']
        v0 = self.params['v0']
        h0 = self.params['h0']
        p0 = self.params['p0']

        state[:,:,0] = rho0
        state[:,:,1] = rho0*u0
        state[:,:,2] = rho0*v0
        state[:,:,3] = rho0*h0 - p0
        p[:,:] = p0

    # calculate pressure
    def __update_pressure(self, workspace, state):
        # retrieve variables
        p = workspace.get_field('p', self.className)
        w = state
        gamma = self.params['gamma']
        pad = self.padding
        nx, ny = workspace.field_size()
        ip, jp = 2, 2
        ie, je = nx+pad, ny+pad

        rqq = 0.5*(w[ip:ie, jp:je, 1]**2 + w[ip:ie, jp:je, 2])/w[ip:ie, jp:je, 0]
        p[ip:ie, jp:je] = (gamma-1)*pos_diff(w[ip:ie, jp:je, 3], rqq)


### TESTING

    def test(self, workspace, state, method, code=''):
        """Calculates the spatial eflux given the current state.
        
        Args:
            workspace (Workspace): contains the relevant fields
            state (Field): the current state
            output (Field): where the flux values will be stored
            """
        assert(isfinite(state))
        self.__check_vars(workspace)
        
        # retrieve necessary workspace fields
        def get(varName):
            return workspace.get_field(varName, self.className)
        w = get("w")
        dw = get("dw")
        bcmodel = self.BCmodel

        # copy state into padded array
        self.__copy_in(state, w)
        dw *= 0

        # update pressure
        self.__update_pressure(workspace, w)
        
        # viscosity tests
        if method=='turb':
            if code=='fortran':
                turb2(self, workspace, w)
            else:
                turbulent_viscosity(self, workspace, w)
            # return ev
            return copy(get('ev'))
        
        if method=='ynot':
            if code=='fortran':
                [ynot,dsti] = thickness_fortran(self, workspace, w)
            else:
                [ynot,dsti] = boundary_thickness(self, workspace, w)
            # return ynot
            return copy(ynot)
        
        if method=='dsti':
            if code=='fortran':
                [ynot,dsti] = thickness_fortran(self, workspace, w)
            else:
                [ynot,dsti] = boundary_thickness(self, workspace, w)
            # return dsti
            return copy(dsti)

        if method=='ev':
            if code=='fortran':
                viscosity_fortran(self, workspace, w)
            else:
                bcmodel.update_physics(self, workspace, w)
            # return ev
            return copy(get('ev'))

        if method=='lv':
            if code=='fortran':
                viscosity_fortran(self, workspace, w)
            else:
                bcmodel.update_physics(self, workspace, w) 
            # return lv
            return copy(get('lv'))
        
        # update viscosity
        bcmodel.update_physics(self, workspace, w)
        
        
        # test stability
        if method=='radi':
            if code=='fortran':
                stability_fortran(bcmodel, self, workspace, w)
            else:
                self.BCmodel.update_stability(self, workspace, w)
            return copy(get('radI'))
        if method=='radj':
            if code=='fortran':
                stability_fortran(bcmodel, self, workspace, w)
            else:
                self.BCmodel.update_stability(self, workspace, w)
            return copy(get('radJ'))
        if method=='rfl':
            if code=='fortran':
                stability_fortran(bcmodel, self, workspace, w)
            else:
                self.BCmodel.update_stability(self, workspace, w)
            return copy(get('rfl'))
        if method=='rfli':
            if code=='fortran':
                stability_fortran(bcmodel, self, workspace, w)
            else:
                self.BCmodel.update_stability(self, workspace, w)
            return copy(get('rfli'))
        if method=='rflj':
            if code=='fortran':
                stability_fortran(bcmodel, self, workspace, w)
            else:
                self.BCmodel.update_stability(self, workspace, w)
            return copy(get('rflj'))
        if method=='dtl':
            if code=='fortran':
                stability_fortran(bcmodel, self, workspace, w)
            else:
                self.BCmodel.update_stability(self, workspace, w)
            return copy(get('dtl'))
        if method=='dtlc':
            if code=='fortran':
                dtlc = stability_fortran(bcmodel, self, workspace, w)
            else:
                dtlc = self.BCmodel.update_stability(self, workspace, w)
            return copy(dtlc)
        
        
        # test boundary conditions
        if method=='bcwall':
            if code=='fortran':
                bcwall_fortran(bcmodel, self, workspace, w)
            else:
                bcmodel.bc_wall(self, workspace, w)
        if method=='bcfar':
            bcmodel.bc_wall(self, workspace, w)
            if code=='fortran':
                bcfar_fortran(bcmodel, self, workspace, w)
            else:
                bcmodel.bc_far(self, workspace, w)
        if method=='halo':
            bcmodel.bc_wall(self, workspace, w)
            # bcmodel.bc_far(self, workspace, w)
            if code=='fortran':
                halo_fortran(bcmodel, self, workspace, w)
            else:
                bcmodel.halo(self, workspace, w)
                
        if method=='bcwall' or method=='bcfar' or method=='halo':
            return copy(w)
        
        # update boundary conditions
        bcmodel.bc_wall(self, workspace, w)
        bcfar_fortran(bcmodel, self, workspace, w)
        bcmodel.halo(self, workspace, w)
        # bcmodel.bc_all(self, workspace, w)
        # w += 1
        # dw += 1

        # calculate residuals
        if method=='eflux':
            if code=='fortran':
                eflux_fortran(self, workspace, w, dw)
            else:
                eflux(self, workspace, w, dw)
        if method=='dflux':
            if code=='fortran':
                dflux_fortran(self, workspace, w, dw, 1)
            else:
                dflux(self, workspace, w, dw, 1)
            dw += 1
        if method=='dfluxc':
            if code=='fortran':
                dfluxc_fortran(self, workspace, w, dw, 1)
            else:
                dfluxc(self, workspace, w, dw, 1)
            dw += 1
        if method=='nsflux':
            if code=='fortran':
                nsflux_fortran(self, workspace, w, dw, 1)
            else:
                nsflux(self, workspace, w, dw, 1)
                
        # return copy of residuals
        return copy(dw)