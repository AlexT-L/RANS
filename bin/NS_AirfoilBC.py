from bin.BoundaryConditioner import BoundaryConditioner
from bin.NS_AirfoilBC_imp import NS_AirfoilBC_imp as implementation

class NS_AirfoilBC(BoundaryConditioner):
    
    
    def __init__(self, input):
        self.padding = 2
        self.local_timestepping = not bool(input.vt)
        self.bc = input.bc

    # initialize state
    def init_state(self, model, workspace, state):
        self.__check_vars(workspace)
        implementation.init_state(self, model, workspace, state)

    # Methods for applying boundary conditions

    # update rev and rlv
    def update_physics(self, model, workspace, state):
        self.__check_vars(workspace)
        implementation.update_physics(self, model, workspace, state)
    
    # update stability
    def update_stability(self, model, workspace, state):
        self.__check_vars(workspace)
        implementation.update_physics(self, model, workspace, state)
    
    # apply far-field boundary conditions
    def bc_far(self, model, workspace, state):
        self.__check_vars(workspace)
        implementation.bc_far(self, model, workspace, state)


    # apply wall boundary conditions
    def bc_wall(self, model, workspace, state):
        self.__check_vars(workspace)
        implementation.bc_wall(self, model, workspace, state)


    # apply halo boundary conditions
    def halo(self, model, workspace, state):
        self.__check_vars(workspace)
        implementation.halo(self, model, workspace, state)


    # apply all boundary conditions
    def bc_all(self, model, workspace, state):
        self.__check_vars(workspace)
        self.bc_wall(self, model, workspace, state)
        self.bc_far(self, model, workspace, state)
        self.halo(self, model, workspace, state)

    # transfer data between workspaces
    def transfer_down(self, model, workspace1, workspace2):
        self.__check_vars(workspace1)
        self.__check_vars(workspace2)
        implementation.transfer_down(self, model, workspace1, workspace2)


    # initialize class workspace fields
    def __init_vars(self, workspace):
        [nx, ny] = workspace.field_size()
        p = self.padding
        field_size = [p+nx+p, p+ny+p]
        grid_size = workspace.grid_size()
        className = self.className

        vars = dict()

        # edge porosities
        vars["pori"] = [grid_size, 1]
        vars["porj"] = [grid_size, 1]

        # stability
        vars["s"] = [field_size, 1]
        vars["dtlc"] = [field_size, 1]

        # Cp and Cf
        vars["cp"] = [[p+nx+p,1], 1]
        vars["cf"] = [[p+nx+p,1], 1]

        workspace.init_vars(className, vars)

        self.__set_porosity(workspace)

    # set the porosity values
    def __set_porosity(self, workspace):
        implementation.set_porosity(self, workspace)
