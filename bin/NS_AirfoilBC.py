from bin.BoundaryConditioner import BoundaryConditioner
from bin.NS_AirfoilBC_imp import NS_AirfoilBC_imp as implementation

class NS_AirfoilBC(BoundaryConditioner):
    
    
    def __init__(self, input):
        pass

# Methods for applying boundary conditions
    
    # apply far-field boundary conditions
    def bc_far(self, workspace, state, fields):

        # extract fields
        p = fields.p

        implementation.bc_far(self, workspace, state, p)


    # apply wall boundary conditions
    def bc_wall(self, workspace, state, fields):

        # extract fields
        p = fields.p
        vol = fields.vol
        rev = fields.rev

        implementation.bc_wall(self, workspace, state, p, vol, rev)


    # apply halo boundary conditions
    def halo(self, workspace, state, fields):

        # extract fields
        p = fields.p
        vol = fields.vol

        implementation.halo(self, workspace, state, p, vol)


    # apply all boundary conditions
    def bc_all(self, workspace, state, fields):
        self.bc_wall(workspace, state, fields)
        self.bc_far(workspace, state, fields)
        self.halo(workspace, state, fields)

    # Get porosity
    def get_pori(self, i, j):
        return self.pori(i,j)

    def get_porj(self, i, j):
        return self.porj(i,j)

    # transfer data between workspaces
    def transfer_down(self, workspace1, workspace2, fields1, fields2):

        # extract fields
        rev1 = fields1.rev
        rlv1 = fields1.rlv
        rev2 = fields2.rev
        rlv2 = fields2.rlv

        implementation.transfer_down(self, workspace1, workspace2, rev1, rlv1, rev2, rlv2)