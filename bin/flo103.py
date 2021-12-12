import Input, AirfoilMap, NavierStokes, Workspace, MultiGrid, CellCenterWS
from bin.ImplicitEuler import ImplicitEuler

if __name__ == '__main__':
    # Comment later
    filename = 'data.dat'
    TOLERANCE = 0.001
    physicsUpdateFrequency = 1
    
    # Command line inputs: Cycle type, Integrator type
    input = Input(filename) # Will actually take all command line inputs
    grid = AirfoilMap(input.grid)
    workspace = CellCenterWS(grid)
    model = NavierStokes(input.model)
    state = model.init_state(workspace)
    integrator = ImplicitEuler(model, input.integrator)
    mg = MultiGrid(workspace, model, integrator, state, input.multigrid)

    # Update model viscosity with current state
    
    
    num_iterations = 0
    while mg.res() < TOLERANCE:
        # update rev and rlv at specified interval
        updatePhysics = True
        if num_iterations != 0:
            if (num_iterations % physicsUpdateFrequency) != 0:
                updatePhysics = False
        if updatePhysics: 
            model.update_physics(workspace, state)
        
        mg.performCycle()
    
    sol = mg.solution()
    
    # Take solution and plot and save info