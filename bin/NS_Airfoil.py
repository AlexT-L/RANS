from numpy.core.defchararray import mod
from bin.BoundaryConditioner import BoundaryConditioner
import bin.model_funcs.bc_transfer as tf
from bin.model_funcs.bc_metric import halo_geom
from bin.model_funcs.bcwall import wall
from bin.model_funcs.bcfar import far_field
from bin.model_funcs.halo import halo
from bin.model_funcs.stability_fast import stability

class NS_Airfoil(BoundaryConditioner):
    """
    Description
    
    Implements boundary conditions for Navier Stokes based model of flow over an airfoil.
    A halo is formed around the mesh containing ghost nodes. Two types of boundary conditions
    are implemented: wall boundaries and far field boundaries

    Attributes
    
    class_name: 
        name of class for accessing Fields in the workspace

    padding: 
        number of ghost nodes on boundary

    local_timestepping: 
        parameter for use in calculating stable time step
        
    bc: 
        some parameter

    Libraries/Modules
    
    numpy
    BoundaryConditioner
    NS_Airfoil_imp
    bcfar
    bcwall
    halo
    stability

    Notes
    
    Based on bcfar.f and bcwall.f """
    
    
    def __init__(self, input):
        """Constructor
        
        Args:
        
        input:
            dictionary of values containing vt and bc

        Returns
        
        A new NS_Arifoil object 

        """
        self.className = "NS_Airfoil"
        self.padding = 2
        self.local_timestepping = not bool(input['vt'])
        self.bc = input['bc']

    # Methods for applying boundary conditions

    # update rev and rlv
    def update_physics(self, model, workspace, state):
        """
        updates the turbulent viscocity for calculation of boundary conditions
        
        Args:
        
        model:
            instance of NavierStokes model class

        workspace:
            instance of workspace class with the relevant fields

        state:
            current state of the system (density, momentum, energy)
    
        """
        ### This method should call Baldwin Lomax ###
        # turbulent_viscosity(params, dims)
    
    # update stability
    def update_stability(self, model, workspace, state):
        """
        updates stability parameters for time step calculations
        
        Args:
        
        model:
            instance of NavierStokes model class

        workspace:
            instance of workspace class with the relevant fields

        state:
            current state of the system (density, momentum, energy)
    
        """
        self.__check_vars(workspace)
        stability(self, model, workspace, state)
    
    # apply far-field boundary conditions
    def bc_far(self, model, workspace, state):
        """
        apply boundary condition in the far field
        
        Args:
        
        model:
            instance of NavierStokes model class

        workspace:
            instance of workspace class with the relevant fields

        state:
            current state of the system (density, momentum, energy)
    
        """
        self.__check_vars(workspace)
        far_field(self, model, workspace, state)
        # implementation.bc_far(self, model, workspace, state)


    # apply wall boundary conditions
    def bc_wall(self, model, workspace, state):
        """
        apply boundary condition along the wall
        
         Args:
         
         model:
            instance of NavierStokes model class

         workspace:
            instance of workspace class with the relevant fields

         state:
            current state of the system (density, momentum, energy)
        
        
        """
        self.__check_vars(workspace)
        wall(self, model, workspace, state)
        # implementation.bc_wall(self, model, workspace, state)


    # apply halo boundary conditions
    def halo(self, model, workspace, state):
        """
        set the values in the ghost cells
        
         Args:
         
         model:
            instance of NavierStokes model class

         workspace:
            instance of workspace class with the relevant fields

         state:
            current state of the system (density, momentum, energy)
        
        
        """
        self.__check_vars(workspace)
        halo(self, model, workspace, state)
        # implementation.halo(self, model, workspace, state)


    # apply all boundary conditions
    def bc_all(self, model, workspace, state):
        """
        do wall boundaries, far field and set halo values at once
        
         Args:
         
         model:
            instance of NavierStokes model class

         workspace:
            instance of workspace class with the relevant fields

         state:
            current state of the system (density, momentum, energy)
        
        """
        self.__check_vars(workspace)
        self.bc_wall(model, workspace, state)
        self.bc_far(model, workspace, state)
        self.halo(model, workspace, state)

    # transfer data between workspaces
    def transfer_down(self, model, workspace1, workspace2):
        self.__check_vars(workspace1)
        self.__check_vars(workspace2)
        tf.transfer_down(self, model, workspace1, workspace2)

    # Get porosity
    def get_pori(self, workspace):
        """
        grab porosity in i direction from the workspace
        
         Args:
         
         workspace:
            instance of workspace class with the relevant field
        
        """
        self.__check_vars(workspace)
        return workspace.get_field("pori", self.className)

    def get_porj(self, workspace):
        """
        grab porosity in j direction from the workspace
        
         Args:
         
         workspace:
            instance of workspace class with the relevant field
        
        """
        self.__check_vars(workspace)
        return workspace.get_field("pori", self.className)

    # set geometry values in the halo
    def halo_geom(self, model, workspace):
        """
        grab porosity in i direction from the workspace
        
         Args:
         
         workspace:
            instance of workspace class with the relevant field
        
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
