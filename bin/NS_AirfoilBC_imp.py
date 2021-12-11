

class NS_AirfoilBC_imp():

    def bc_far(self, this, workspace, state, p):
        # get geometry dictionary
        geom = workspace.get_geom()


    def bc_wall(self, this, workspace, state, p, vol, rev):
        # get geometry dictionary
        geom = workspace.get_geom()


    def halo(self, this, workspace, state, p, vol):
        # get geometry dictionary
        geom = workspace.get_geom()

    
    def transfer_down(self, this, workspace1, workspace2, rev1, rlv1, rev2, rlv2):
        # get geometry dictionary
        geom1 = workspace1.get_geom()
        geom2 = workspace2.get_geom()