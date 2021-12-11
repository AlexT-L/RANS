from Model import Model
from Workspace import Workspace
from model_funcs import eflux_wrap, nsflux_wrap

class NavierStokes(Model):
    
    def __init__(self, bcmodel, input):
        self.BCmodel = bcmodel

    # flux calculations
    from .model_funcs import eflux_wrap,nsflux_wrap



    def get_flux(self, workspace, state, output, update_factor=1):
        pass


    def get_safe_timestep(self, workspace, state, dt):
        pass


    def update_stability(self, workspace, state):
        pass

    
    def transfer_down(self, workspace1, workspace2):
        # This needs to be filled in

        self.BCmodel.transfer_down(workspace1, workspace2, fields1, fields2)
