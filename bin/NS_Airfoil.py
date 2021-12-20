from numpy.core.defchararray import mod
from bin.BoundaryConditioner import BoundaryConditioner
import bin.model_funcs.bc_transfer as tf
from bin.model_funcs.bc_metric import halo_geom
from bin.model_funcs.bcwall import wall
from bin.model_funcs.bcfar import far_field
from bin.model_funcs.halo import halo
from bin.model_funcs.stability_fast import stability

class NS_Airfoil(BoundaryConditioner):
    
    
    def __init__(self, input):
        """Initializes NS_Airfoil object.
        
        Parameters
        ----------
        input:
            Dictionary of input parameters that must contain values for:
                - vt
                - bc
        """
        self.className = "NS_Airfoil"
        self.padding = 2
        self.local_timestepping = not bool(input['vt'])
        self.bc = input['bc']

    # update ev and lv
    def update_physics(self, model, workspace, state):
        """Updates laminar viscosity and eddy viscosity
        
        Parameters
        ----------
        model:
            The physics model
        workspace:
            The Workspace
        state:
            A field containing the current state
        """
        self.__check_vars(workspace)
        ### This method should call Baldwin Lomax ###
        # turbulent_viscosity(params, dims)
    
    # update stability
    def update_stability(self, model, workspace, state):
        """Updates local time step bounds to ensure stability
        
        Parameters
        ----------
        model:
            The physics model
        workspace:
            The Workspace
        state:
            A field containing the current state
        """
        self.__check_vars(workspace)
        stability(self, model, workspace, state)
    
    # apply far-field boundary conditions
    def bc_far(self, model, workspace, state):
        """Sets the far-field boundary conditions
        
        Parameters
        ----------
        model:
            The physics model
        workspace:
            The Workspace
        state:
            A field containing the current state
        """
        self.__check_vars(workspace)
        far_field(self, model, workspace, state)

    # apply wall boundary conditions
    def bc_wall(self, model, workspace, state):
        """Sets the wall boundary conditions
        
        Parameters
        ----------
        model:
            The physics model
        workspace:
            The Workspace
        state:
            A field containing the current state
        """
        self.__check_vars(workspace)
        wall(self, model, workspace, state)

    # apply halo boundary conditions
    def halo(self, model, workspace, state):
        """Sets the values in the outer halo
        
        Parameters
        ----------
        model:
            The physics model
        workspace:
            The Workspace
        state:
            A field containing the current state
        """
        self.__check_vars(workspace)
        halo(self, model, workspace, state)

    # apply all boundary conditions
    def bc_all(self, model, workspace, state):
        """Sets all boundary conditions (wall, far-field, outer halo)
        
        Parameters
        ----------
        model:
            The physics model
        workspace:
            The Workspace
        state:
            A field containing the current state
        """
        self.__check_vars(workspace)
        self.bc_wall(model, workspace, state)
        self.bc_far(model, workspace, state)
        self.halo(model, workspace, state)

    # transfer data between workspaces
    def transfer_down(self, model, workspace1, workspace2):
        """returns the porosity in the i diretion
        
        Parameters
        ----------
        model:
            The physics model
        workspace1:
            The Workspace corresponding to the finer grid
        workspace2:
            The Workspace corresponding to the coarser grid
        """
        self.__check_vars(workspace1)
        self.__check_vars(workspace2)
        tf.transfer_down(self, model, workspace1, workspace2)

    # Get porosity
    def get_pori(self, workspace):
        """returns the porosity in the i diretion
        
        Parameters
        ----------
        workspace:
            The Workspace
        """
        self.__check_vars(workspace)
        return workspace.get_field("pori", self.className)

    def get_porj(self, workspace):
        """returns the porosity in the j diretion
        
        Parameters
        ----------
        workspace:
            The Workspace
        """
        self.__check_vars(workspace)
        return workspace.get_field("porj", self.className)

    # set geometry values in the halo
    def halo_geom(self, model, workspace):
        """Sets geometry values in the halo
        
        Parameters
        ----------
        model:
            The physics model
        workspace:
            The Workspace
        """
        self.__check_vars(workspace)
        halo_geom(self, model, workspace)

    # check if dictionary has been initialized
    def __check_vars(self, workspace):
        if not workspace.has_dict(self.className):
            self.__init_vars(workspace)

    # initialize class workspace fields
    def __init_vars(self, workspace):
        [nx, ny] = workspace.field_size()
        p = self.padding
        field_size = [p+nx+p, p+ny+p]

        vars = dict()

        # edge porosities
        vars["pori"] = [(nx+1, ny)]
        vars["porj"] = [(nx, ny+1)]

        # stability
        vars["s"] = [field_size]
        vars["dtlc"] = [field_size]

        # Cp and Cf
        vars["cp"] = [p+nx+p]
        vars["cf"] = [p+nx+p]

        workspace.init_vars(self.className, vars)

        self.__set_porosity(workspace)

    # set the porosity values
    def __set_porosity(self, workspace):
        # get relevant geometry parameters
        dims = workspace.get_dims()
        itl = dims['itl']
        itu = dims['itu']

        # get porosity
        pori = workspace.get_field("pori", self.className)
        porj = workspace.get_field("porj", self.className)

        # set the porosity to unity
        pori[:] = 1.0
        porj[:] = 1.0

        # flag the wall at the j boundaries
        porj[itl:itu,0]   = 0.0
