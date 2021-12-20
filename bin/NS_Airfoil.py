from numpy.core.defchararray import mod
from bin.BoundaryConditioner import BoundaryConditioner
import bin.utils.NS_Airfoil_imp as implementation
from bin.model_funcs.bcfar import far_field
from bin.model_funcs.halo import halo
from bin.model_funcs.bcwall import wall
# from bin.model_funcs.stability import stability
from bin.model_funcs.stability_fast import stability

class NS_Airfoil(BoundaryConditioner):
    
    
    def __init__(self, input):
        self.className = "NS_Airfoil"
        self.padding = 2
        self.local_timestepping = not bool(input['vt'])
        self.bc = input['bc']

    # Methods for applying boundary conditions

    # update rev and rlv
    def update_physics(self, model, workspace, state):
        self.__check_vars(workspace)
        implementation.update_physics(self, model, workspace, state)
    
    # update stability
    def update_stability(self, model, workspace, state):
        self.__check_vars(workspace)
        stability(self, model, workspace, state)
    
    # apply far-field boundary conditions
    def bc_far(self, model, workspace, state):
        self.__check_vars(workspace)
        far_field(self, model, workspace, state)
        # implementation.bc_far(self, model, workspace, state)


    # apply wall boundary conditions
    def bc_wall(self, model, workspace, state):
        self.__check_vars(workspace)
        wall(self, model, workspace, state)
        # implementation.bc_wall(self, model, workspace, state)


    # apply halo boundary conditions
    def halo(self, model, workspace, state):
        self.__check_vars(workspace)
        halo(self, model, workspace, state)
        # implementation.halo(self, model, workspace, state)


    # apply all boundary conditions
    def bc_all(self, model, workspace, state):
        self.__check_vars(workspace)
        self.bc_wall(model, workspace, state)
        self.bc_far(model, workspace, state)
        self.halo(model, workspace, state)

    # transfer data between workspaces
    def transfer_down(self, model, workspace1, workspace2):
        self.__check_vars(workspace1)
        self.__check_vars(workspace2)
        implementation.transfer_down(self, model, workspace1, workspace2)

    # Get porosity
    def get_pori(self, workspace):
        self.__check_vars(workspace)
        return workspace.get_field("pori", self.className)

    def get_porj(self, workspace):
        self.__check_vars(workspace)
        return workspace.get_field("pori", self.className)

    # set geometry values in the halo
    def halo_geom(self, model, workspace):
        self.__check_vars(workspace)
        implementation.halo_geom(self, model, workspace)

    # check if dictionary has been initialized
    def __check_vars(self, workspace):
        if not workspace.has_dict(self.className):
            self.__init_vars(workspace)

    # initialize class workspace fields
    def __init_vars(self, workspace):
        [nx, ny] = workspace.field_size()
        p = self.padding
        field_size = [p+nx+p, p+ny+p]
        grid_size = workspace.grid_size()
        className = self.className

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

        workspace.init_vars(className, vars)

        self.__set_porosity(workspace)

    # set the porosity values
    def __set_porosity(self, workspace):
        implementation.set_porosity(self, workspace)
