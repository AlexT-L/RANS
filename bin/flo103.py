import Input, SqrtGrid, NavierStokes, Workspace, MultiGrid
from bin.ImplicitEuler import ImplicitEuler

if __name__ == '__main__':
    # Comment later
    filename = 'data.dat'
    TOLERANCE = 0.001
    
    # Command line inputs: Cycle type, Integrator type
    input = Input(filename) # Will actually take all command line inputs
    grid = SqrtGrid(input.grid)
    workspace = Workspace(grid)
    model = NavierStokes(input.model)
    integrator = ImplicitEuler(model, input.integrator)
    mg = MultiGrid(workspace, model, integrator, input.multigrid)
    
    while mg.res() < TOLERANCE:
        mg.loop()
    
    sol = mg.solution()
    
    # Take solution and plot and save info